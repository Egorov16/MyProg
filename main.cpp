#include <algorithm>
#include <array>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math.h>
#include <set>
#include <string>
#include <vector>

// Структура для узла
struct Node {
  double x;
  double y;
  double z;
};

// Точки интегрирования (по умолчанию, для куба, центр в (0, 0, 0))
const double commonCrd1 = 0.57735026918962573;
const std::array<Node, 8> integrationPoints{
    Node{-commonCrd1, -commonCrd1, -commonCrd1}, // Точка 0
    Node{-commonCrd1, -commonCrd1, commonCrd1},  // Точка 1
    Node{-commonCrd1, commonCrd1, -commonCrd1},  // Точка 2
    Node{-commonCrd1, commonCrd1, commonCrd1},   // Точка 3
    Node{commonCrd1, -commonCrd1, -commonCrd1},  // Точка 4
    Node{commonCrd1, -commonCrd1, commonCrd1},   // Точка 5
    Node{commonCrd1, commonCrd1, -commonCrd1},   // Точка 6
    Node{commonCrd1, commonCrd1, commonCrd1}     // Точка 7
};

using indexes = std::vector<size_t>;
using Crd = std::vector<Node>;
enum class Type { hexa8 = 0, hexa20 = 1 };
// Структура для конечного элемента
struct Element {
  Type type;
  indexes inds;
};
using Elems = std::vector<Element>;
// Структура для хранения векторов
struct SolutionData {
  std::vector<float> Temperature;
  std::vector<float> SigmaXX;
  std::vector<float> SigmaYY;
  std::vector<float> SigmaZZ;
  std::vector<float> SigmaXY;
  std::vector<float> SigmaYZ;
  std::vector<float> SigmaZX;
  std::vector<float> DefXX;
  std::vector<float> DefYY;
  std::vector<float> DefZZ;
  std::vector<float> DefXY;
  std::vector<float> DefYZ;
  std::vector<float> DefZX;
};

// Зачем две похожие структуры?
//// Структура для хранения данных в точках интегрирования
// struct IntegrationPointData {
//   int elementNumber;
//   int gaussPointNumber;
//   double localCoords[3]; // Массив для хранения локальных координат
//   double globalCoords[3]; // Массив для хранения локальных координат
//   double stress[6]; // Массив для хранения напряжений
// };

enum class ValueType { CrdX, CrdY, CrdZ, SXX, SYY, SZZ, SXY, SYZ, SZX };

struct GaussPointData {
  size_t elementNumber;
  size_t gaussPointNumber;
  Node localCoords; // Массив для хранения локальных координат
  Node globalCoords; // Массив для хранения локальных координат
  std::array<double, 6> stress; // Массив для хранения напряжений
  std::array<double, 6> strain; // Массив для хранения деформаций

  double getValue(ValueType type) {
    switch (type) {
      // TODO Дописать остальное
    case ValueType::SXX:
      return stress[0];
      break;
    default:
      return std::numeric_limits<double>::quiet_NaN();
      break;
    }
  }
};

struct Mesh {
  Crd crd;
  Elems elems;
  SolutionData solutionData;
  // Внешний вектор - по элементам, внутренний по точкам
  std::vector<std::vector<GaussPointData>> intData;
  // Карта принадлежности элементов узлам
  std::map<size_t, std::vector<size_t>> nodeElements;
};

// Чтение сетки и данных
Mesh readMesh(const std::string &prefix);
// Проверка интерполяций
void checkData(const Mesh &mesh);

// Функция формы
// double shapeFunction(double localCoord, double xi) {
//  // Линейная функция формы
//  return 0.5 * (1.0 + xi * localCoord);
//}

// Получить локальный индекс точки интегрирования, ближайшей к узлу с глобальным
// номером nodeID
size_t getNearestGPID(const Mesh &mesh, size_t elemID, size_t nodeID) {
  assert(mesh.elems[elemID].type == Type::hexa8);

  const auto &ids = mesh.elems[elemID].inds;
  auto pos = std::find(ids.cbegin(), ids.cend(), nodeID);
  assert(pos != ids.cend());
  auto localID = pos - ids.cbegin();
  switch (localID) {
  case 0:
    return 0;
    break;
  case 1:
    return 4;
    break;
  case 2:
    return 6;
    break;
  case 3:
    return 2;
    break;
  case 4:
    return 1;
    break;
  case 5:
    return 5;
    break;
  case 6:
    return 7;
    break;
  case 7:
    return 3;
    break;
  default:
    throw std::runtime_error("Should't happen");
    break;
  }
}

// Поиск ближайших точек интегрировани
auto findNearestGPs(const Mesh &mesh, size_t nodeID) {
  std::vector<std::pair<size_t /*elemID*/, std::vector<size_t> /*gPointIDs*/>>
      nearest;

  auto elems = mesh.nodeElements.at(nodeID);
  assert(elems.size());

  nearest.resize(elems.size());

  for (size_t ct = 0; ct < elems.size(); ++ct) {
    nearest[ct] = {elems[ct], std::vector<size_t>(
                                  {getNearestGPID(mesh, elems[ct], nodeID)})};
  }

  return nearest;
}

double calculateValue(Node localCrd, const std::vector<double> &nodalValues,
                      const Element &element) {

  assert(element.inds.size() == nodalValues.size());
  // TODO написать другие типы
  assert(element.type == Type::hexa8);

  std::vector<double> h(8);

  double onePlsX = (1. + localCrd.x);
  double oneMnsX = (1. - localCrd.x);
  double onePlsY = (1. + localCrd.y);
  double oneMnsY = (1. - localCrd.y);
  double onePlsZ = (1. + localCrd.z);
  double oneMnsZ = (1. - localCrd.z);

  // Вычисляем значения функций формы
  h[0] = oneMnsX * oneMnsY * oneMnsZ / 8.0;
  h[1] = onePlsX * oneMnsY * oneMnsZ / 8.0;
  h[2] = onePlsX * onePlsY * oneMnsZ / 8.0;
  h[3] = oneMnsX * onePlsY * oneMnsZ / 8.0;
  h[4] = oneMnsX * oneMnsY * onePlsZ / 8.0;
  h[5] = onePlsX * oneMnsY * onePlsZ / 8.0;
  h[6] = onePlsX * onePlsY * onePlsZ / 8.0;
  h[7] = oneMnsX * onePlsY * onePlsZ / 8.0;

  // Вычисляем значение в точке
  double value = 0;
  for (size_t ct = 0, size = h.size(); ct < size; ++ct)
    value += h[ct] * nodalValues[ct];
  return value;
};

// Функция для расчета расстояния между двумя точками
double distance(const Node &point1, const Node &point2) {
  return std::sqrt(std::pow(point1.x - point2.x, 2) +
                   std::pow(point1.y - point2.y, 2) +
                   std::pow(point1.z - point2.z, 2));
}

//// Функция для поиска ближайшей точки интегрирования
// int findClosestIntegrationPoint(const Node &targetNode,
//                                 const Element &element) {
//   double minDistance = std::numeric_limits<double>::max();
//   int closestPointIndex = 0;

//  for (size_t i = 0; i < element.integrationPoints.size(); ++i) {
//    double currentDistance = distance(targetNode,
//    element.integrationPoints[i]); if (currentDistance < minDistance) {
//      minDistance = currentDistance;
//      closestPointIndex = i;
//    }
//  }

//  return closestPointIndex;
//}

// Собрать вектор узловых значений
std::vector<double> constructNodalData(const Mesh &mesh, size_t elemID,
                                       ValueType vType) {
  std::vector<double> data;
  assert(mesh.elems[elemID].type == Type::hexa8);

  data.resize(8);
  const auto &ids = mesh.elems[elemID].inds;
  for (size_t ct = 0; ct < ids.size(); ++ct) {
    if (vType == ValueType::CrdX)
      data[ct] = mesh.crd[ids[ct]].x;
    else if (vType == ValueType::CrdY)
      data[ct] = mesh.crd[ids[ct]].y;
    else if (vType == ValueType::CrdZ)
      data[ct] = mesh.crd[ids[ct]].z;
    else if (vType == ValueType::SXX)
      data[ct] = mesh.solutionData.SigmaXX[ids[ct]];
    else // TODO дописать остальные поля
      assert(false);
  }
  return data;
}

// Вычисление невязки в отдельной точке
// Возвращает абсолютную и относительную погрешность
std::pair<double, double> calculateError(const Mesh &mesh, size_t elemId,
                                         size_t pointId, double targetValue,
                                         std::vector<double> nodeValues) {
  auto value = calculateValue(integrationPoints[pointId], nodeValues,
                              mesh.elems.at(elemId));

  auto err = value - targetValue;
  double relErr;
  // TODO Надо провирить на пограничных случаях
  if (std::fabs(targetValue) > std::numeric_limits<double>::epsilon())
    relErr = err / targetValue;
  else {
    relErr = 0;
  }
  return {err, relErr};
}

// Функция для вычисления невязки
std::vector<std::vector<double>>
calculateDiscrepancy(const std::vector<GaussPointData> &interpolatedData,
                     const std::vector<GaussPointData> &fileData) {

  // Внешний вектор для хранения невязок по элементам
  // Внутренний вектор для хранения невязок по компонентам напряжения в элементе
  std::vector<std::vector<double>> discrepancy;

  // Проверка, что вектора одинаковой длины
  if (interpolatedData.size() != fileData.size()) {
    std::cerr << "Вектора данных имеют разную длину!\n";
    return discrepancy; // Возвращаем пустой вектор в случае ошибки
  }

  // Итерация по данным и вычисление невязки
  int currentElement = -1;
  // Текущий номер элемента
  std::vector<double> elementDiscrepancy;
  // Вектор для хранения невязок для текущего элемента

  for (size_t i = 0; i < interpolatedData.size(); ++i) {
    const auto &interpPoint = interpolatedData[i];
    const auto &filePoint = fileData[i];

    // Проверка на совпадение номеров элемента и точки Гаусса
    if (interpPoint.elementNumber != filePoint.elementNumber ||
        interpPoint.gaussPointNumber != filePoint.gaussPointNumber) {
      std::cerr << "Несовпадение номеров элемента или точки Гаусса!\n";
      continue; // Переходим к следующей точке
    }

    // Если это новый элемент, добавляем предыдущий в общий вектор
    if (currentElement != interpPoint.elementNumber) {
      if (!elementDiscrepancy.empty()) {
        discrepancy.push_back(elementDiscrepancy);
      }
      elementDiscrepancy.clear(); // Очищаем вектор невязок для нового элемента
      currentElement = interpPoint.elementNumber;
    }

    // Вычисление невязки для каждого компонента напряжения
    for (int j = 0; j < 6; ++j) {

      double k = 100 * std::abs(interpPoint.stress[j] - filePoint.stress[j]) /
                 std::max(std::abs(interpPoint.stress[j]),
                          std::abs(filePoint.stress[j]));

      elementDiscrepancy.push_back(k);
    }
  }

  return discrepancy;
}

// Функция пробует улучшить узловое значение в узле nodeID
void tryToImproove(const Mesh &mesh, size_t nodeId, ValueType type) {
  // Найти все узлы, соседние с целевым

  // Найти все точки интегрирования, ближайшие к целевому узлу

  // Создать функцию невязки

  // Варьирование
}

int main(int numArgs, char **args) {
  // Не понимаю зачем это
  setlocale(LC_ALL, "rus");

  // Путь по умолчанию
  std::string workDir = "";
  // В аргументах передан каталог
  if (numArgs == 2) {
    // Передан рабочий каталог
    workDir = std::string(args[1]) + "/";
  }
  // В аргументах явно ошибка
  if (numArgs > 2)
    throw std::runtime_error("Wrong number of params. Should be 0 or 1 (path)");

  const auto mesh = readMesh(workDir);
  auto &crd = mesh.crd;
  auto elems = mesh.elems;
  auto &solutionData = mesh.solutionData;
  const size_t numNodes = crd.size();

  checkData(mesh);

  // Проверяем напряжения SXX
  double maxErr = 0, maxRelErr = 0;
  for (size_t elemCt = 0, size = mesh.elems.size(); elemCt < size; ++elemCt) {
    auto data = constructNodalData(mesh, elemCt, ValueType::SXX);
    for (auto &gPoint : mesh.intData[elemCt]) {
      auto val = calculateValue(gPoint.localCoords, data, mesh.elems[elemCt]);
      auto target = gPoint.stress[0];
      auto err = target - val;
      auto relErr = err / target;
      if (std::fabs(err) > std::fabs(maxErr))
        maxErr = err;
      if (std::fabs(relErr) > std::fabs(maxRelErr))
        maxRelErr = relErr;
    }
  }
  std::cout << "SXX:\t" << maxErr << "\t" << maxRelErr << std::endl;

  //   Вектор для хранения данных в точках интегрирования
  //  std::vector<IntegrationPointData> integrationPointData;

  //  for (int targetNode = 1; targetNode <= crd.size(); ++targetNode) {

  //    for (int elementIndex : nodeElements[targetNode]) {

  //      // Получаем номера узлов из elems:
  //      auto nodeIndexes = elems[elementIndex - 1].second;

  //      // Создаем объект Element:
  //      Element elem = {Type::hexa8, nodeIndexes};

  //      // Индекс ближайшей точки интегрирования
  //      int locNumNode = -1; // Индекс числа, по умолчанию -1 (не найдено)
  //      for (int i = 0; i < 8; ++i) {
  //        if (nodeIndexes[i] == targetNode) {
  //          locNumNode = i;
  //          break; // Прерываем цикл, как только найдем число
  //        }
  //      }
  //      int closestPointIndex;
  //      switch (locNumNode) {
  //      case 0:
  //        closestPointIndex = 0;
  //        break;
  //      case 1:
  //        closestPointIndex = 4;
  //        break;
  //      case 2:
  //        closestPointIndex = 6;
  //        break;
  //      case 3:
  //        closestPointIndex = 2;
  //        break;
  //      case 4:
  //        closestPointIndex = 1;
  //        break;
  //      case 5:
  //        closestPointIndex = 5;
  //        break;
  //      case 6:
  //        closestPointIndex = 7;
  //        break;
  //      case 7:
  //        closestPointIndex = 3;
  //        break;
  //      }

  //      // Получаем координаты точки интегрирования:
  //      double localX = elem.integrationPoints[closestPointIndex].x;
  //      double localY = elem.integrationPoints[closestPointIndex].y;
  //      double localZ = elem.integrationPoints[closestPointIndex].z;

  //      // Создаём вектора узловых значений координат
  //      std::vector<double> coordsNodeX;
  //      std::vector<double> coordsNodeY;
  //      std::vector<double> coordsNodeZ;

  //      for (int i = 0; i < 8; i++) {
  //        int nodeIndex = nodeIndexes[i];
  //        coordsNodeX.push_back(crd[nodeIndex - 1].x);
  //        coordsNodeY.push_back(crd[nodeIndex - 1].y);
  //        coordsNodeZ.push_back(crd[nodeIndex - 1].z);
  //      }

  //      // Глобальные координаты точки интегрирования
  //      double globalX =
  //          calculateValue(localX, localY, localZ, coordsNodeX, elem);
  //      double globalY =
  //          calculateValue(localX, localY, localZ, coordsNodeY, elem);
  //      double globalZ =
  //          calculateValue(localX, localY, localZ, coordsNodeZ, elem);

  //      // Создаем массив для хранения интерполированных значений
  //      double stress[6];

  //      // Итерация по компонентам напряжения и интерполяция
  //      for (int i = 0; i < 6; ++i) {
  //        // Получаем значения в узлах элемента для текущего компонента
  //        std::vector<double> nodalValues;
  //        for (int node : nodeIndexes) {
  //          // Используем правильный доступ к компонентам solutionData
  //          switch (i) {
  //          case 0:
  //            nodalValues.push_back(solutionData.SigmaXX[node - 1]);
  //            break;
  //          case 1:
  //            nodalValues.push_back(solutionData.SigmaYY[node - 1]);
  //            break;
  //          case 2:
  //            nodalValues.push_back(solutionData.SigmaZZ[node - 1]);
  //            break;
  //          case 3:
  //            nodalValues.push_back(solutionData.SigmaXY[node - 1]);
  //            break;
  //          case 4:
  //            nodalValues.push_back(solutionData.SigmaYZ[node - 1]);
  //            break;
  //          case 5:
  //            nodalValues.push_back(solutionData.SigmaZX[node - 1]);
  //            break;
  //          }
  //        }

  //        // Интерполяция значения в ближайшую точку интегрирования
  //        stress[i] = calculateValue(localX, localY, localZ, nodalValues,
  //        elem);
  //      }
  //      // Добавляем данные в вектор integrationPointData
  //      IntegrationPointData data;
  //      data.elementNumber = elementIndex;
  //      data.gaussPointNumber = closestPointIndex;
  //      data.localCoords[0] = localX;
  //      data.localCoords[1] = localY;
  //      data.localCoords[2] = localZ;
  //      data.globalCoords[0] = globalX;
  //      data.globalCoords[1] = globalY;
  //      data.globalCoords[2] = globalZ;

  //      // Копируем интерполированные значения в data.stress
  //      std::copy(std::begin(stress), std::end(stress),
  //      std::begin(data.stress));

  //      integrationPointData.push_back(data);
  //    }
  //  }

  //  std::sort(integrationPointData.begin(), integrationPointData.end(),
  //            [](const IntegrationPointData &a, const IntegrationPointData &b)
  //            {
  //              // Сначала сортируем по номеру элемента
  //              if (a.elementNumber != b.elementNumber) {
  //                return a.elementNumber < b.elementNumber;
  //              } else {
  //                // Если элементы одинаковые, сортируем по номеру точки
  //                // интегрирования
  //                return a.gaussPointNumber < b.gaussPointNumber;
  //              }
  //            });

  //  // Вычисление невязки
  //  std::vector<std::vector<double>> discrepancy =
  //      calculateDiscrepancy(integrationPointData, data);

  //  // Поиск максимального значения невязки
  //  double maxDiscrepancy = 0.0;
  //  for (const auto &elementDiscrepancy : discrepancy) {
  //    for (const auto &value : elementDiscrepancy) {
  //      maxDiscrepancy = std::max(maxDiscrepancy, value);
  //    }
  //  }

  //  std::cout << "\nМаксимальная невязка: " << maxDiscrepancy << "%" <<
  //  std::endl;

  return 0;
}

Mesh readMesh(const std::string &prefix) {
  Mesh mesh;
  auto &crd = mesh.crd;
  // crd
  {
    std::ifstream in(prefix + "crd.sba", std::ifstream::binary);
    if (!in.good())
      throw std::runtime_error("Cant open Crd file: " + prefix + "crd.sba");
    int numNodes, kort;
    in.read((char *)&numNodes, sizeof(int));
    in.read((char *)&kort, sizeof(int));
    numNodes /= kort;
    crd.resize(numNodes);
    std::cout << numNodes << "\t" << kort << std::endl;
    std::vector<float> tmp(numNodes * kort);
    in.read((char *)tmp.data(), sizeof(float) * tmp.size());

    for (size_t ct = 0; ct < numNodes; ++ct) {
      crd[ct].x = tmp[0 * numNodes + ct];
      crd[ct].y = tmp[1 * numNodes + ct];
      crd[ct].z = tmp[2 * numNodes + ct];
    }
  }
  auto &elems = mesh.elems;
  // ind
  {
    std::ifstream in(prefix + "ind.sba", std::ifstream::binary);
    if (!in.good())
      throw std::runtime_error("Cant open Ind file: " + prefix + "ind.sba");
    int numTypes;
    in.read((char *)&numTypes, sizeof(int));
    std::vector<std::pair<int, int>> types(numTypes);
    for (auto &t : types) {
      in.read((char *)&(t.first), sizeof(int));
      in.read((char *)&(t.second), sizeof(int));
    }

    // Общая длина читается один раз
    int numIndexes;
    in.read((char *)&numIndexes, sizeof(int));

    int indexElem = 0;
    for (auto t : types) {
      // TODO добавить остальные типы
      Type type = Type::hexa8;
      size_t indLength;
      if (t.first == 0) {
        type = Type::hexa8;
        indLength = 8;
      } else if (t.first == 1) {
        type = Type::hexa20;
        indLength = 20;
      }

      std::vector<int> ind;
      ind.resize(indLength);

      for (int i = 0; i < t.second; ++i) {
        in.read((char *)ind.data(), sizeof(int) * ind.size());
        for (auto &i : ind)
          --i;
        Element elem;
        elem.type = type;
        elem.inds.insert(elem.inds.begin(), ind.begin(), ind.end());
        elems.push_back(elem);
      }
    }

    // Заполняем вектор nodeElements
    for (size_t i = 0; i < elems.size(); ++i) {
      for (const auto &node : elems[i].inds) {
        mesh.nodeElements[node].push_back(i); // Добавляем элемент в вектор узла
      }
    }
  }
  auto numNodes = crd.size();
  auto &solutionData = mesh.solutionData;
  {
    std::ifstream in(prefix + "solution0001.sba", std::ifstream::binary);

    if (!in.good())
      std::runtime_error("Ошибка открытия файла solution0001.sba");

    // Читаем флаги
    std::vector<int> flags(20);
    in.read((char *)(flags.data()), sizeof(int) * flags.size());

    // Читаем массивы (если они присутствуют)
    for (size_t i = 0; i < 20; ++i) {
      if (flags[i] == 1) {
        // Читаем массив
        std::vector<float> data(numNodes);
        in.read((char *)(data.data()), sizeof(float) * numNodes);
        // Заполняем соответствующий вектор в SolutionData
        switch (i) {
        case 0:
          solutionData.Temperature = data;
          break;
        case 1:
          solutionData.SigmaXX = data;
          break;
        case 2:
          solutionData.SigmaYY = data;
          break;
        case 3:
          solutionData.SigmaZZ = data;
          break;
        case 4:
          solutionData.SigmaXY = data;
          break;
        case 5:
          solutionData.SigmaYZ = data;
          break;
        case 6:
          solutionData.SigmaZX = data;
          break;
        }
      }
    }
  }
  {
    auto &data = mesh.intData;
    data.reserve(mesh.elems.size());

    std::ifstream in(prefix + "int.dat", std::ifstream::binary);

    if (!in.good())
      throw std::runtime_error("Cant open data file: " + prefix + "int.dat");

    size_t lastReadElemID = std::numeric_limits<size_t>::max();
    while (!in.eof()) {
      GaussPointData point;
      in.read(reinterpret_cast<char *>(&point.elementNumber), sizeof(size_t));
      in.read(reinterpret_cast<char *>(&point.gaussPointNumber),
              sizeof(size_t));
      // Считываем локальные координаты в структуру
      in.read((char *)(&point.localCoords.x), sizeof(double));
      in.read((char *)(&point.localCoords.y), sizeof(double));
      in.read((char *)(&point.localCoords.z), sizeof(double));
      // Считываем глобальные координаты в структуру
      in.read((char *)(&point.globalCoords.x), sizeof(double));
      in.read((char *)(&point.globalCoords.y), sizeof(double));
      in.read((char *)(&point.globalCoords.z), sizeof(double));
      // Считываем напряжения в массив
      in.read((char *)(point.stress.data()), sizeof(double) * 6);
      // Считываем деформации в массив
      in.read((char *)(point.strain.data()), sizeof(double) * 6);

      if (lastReadElemID != point.elementNumber) {
        data.push_back(std::vector<GaussPointData>());
        lastReadElemID = point.elementNumber;
      }
      data.back().push_back(point);
    }

    in.close();
  }

  assert(mesh.elems.size() == mesh.intData.size());

  return mesh;
}

void checkData(const Mesh &mesh) {
  // Проверяем консистентность точек интегрирования интерполяцией координат
  for (size_t elemCt = 0, sizeElems = mesh.elems.size(); elemCt < sizeElems;
       ++elemCt) {
    const auto &elem = mesh.elems.at(elemCt);
    assert(elem.type == Type::hexa8);

    // проход по точкам интегрирования внутри элемента
    std::vector<double> nodesData;
    nodesData.resize(8);
    for (size_t pointCt = 0; pointCt < 8; ++pointCt) {
      const auto &gPoint = mesh.intData.at(elemCt).at(pointCt);
      // X
      assert(std::fabs(gPoint.localCoords.x - integrationPoints[pointCt].x) <
             std::numeric_limits<double>::epsilon());
      nodesData = constructNodalData(mesh, elemCt, ValueType::CrdX);
      auto x = calculateValue(integrationPoints[pointCt], nodesData, elem);
      assert(std::fabs(x - gPoint.globalCoords.x) <
             std::numeric_limits<double>::epsilon());
      // Y
      assert(std::fabs(gPoint.localCoords.y - integrationPoints[pointCt].y) <
             std::numeric_limits<double>::epsilon());
      nodesData = constructNodalData(mesh, elemCt, ValueType::CrdY);
      auto y = calculateValue(integrationPoints[pointCt], nodesData, elem);
      assert(std::fabs(y - gPoint.globalCoords.y) <
             std::numeric_limits<double>::epsilon());
      // Z
      assert(std::fabs(gPoint.localCoords.z - integrationPoints[pointCt].z) <
             std::numeric_limits<double>::epsilon());
      nodesData = constructNodalData(mesh, elemCt, ValueType::CrdZ);
      auto z = calculateValue(integrationPoints[pointCt], nodesData, elem);
      assert(std::fabs(z - gPoint.globalCoords.z) <
             std::numeric_limits<double>::epsilon());
    }
  }
}
