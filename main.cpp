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

#include <dlib/global_optimization.h>
#include <dlib/optimization.h>

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
  std::vector<float> SigmaI;
  std::vector<float> Sigma1;
  std::vector<float> Sigma2;
  std::vector<float> Sigma3;
  std::vector<float> DefXX;
  std::vector<float> DefYY;
  std::vector<float> DefZZ;
  std::vector<float> DefXY;
  std::vector<float> DefYZ;
  std::vector<float> DefZX;
  std::vector<float> DefI;
  std::vector<float> DefPlast;
  std::vector<float> None;
};

enum class ValueType {
  CrdX,
  CrdY,
  CrdZ,
  SXX,
  SYY,
  SZZ,
  SXY,
  SYZ,
  SZX,
  EXX,
  EYY,
  EZZ,
  EXY,
  EYZ,
  EZX,
};

struct GaussPointData {
  size_t elementNumber;
  size_t gaussPointNumber;
  Node localCoords; // Массив для хранения локальных координат
  Node globalCoords; // Массив для хранения локальных координат
  std::array<double, 6> stress; // Массив для хранения напряжений
  std::array<double, 6> strain; // Массив для хранения деформаций

  double getValue(ValueType type) const {
    switch (type) {
      // TODO Дописать остальное
    case ValueType::SXX:
      return stress[0];
      break;
    case ValueType::SYY:
      return stress[1];
      break;
    case ValueType::SZZ:
      return stress[2];
      break;
    case ValueType::SXY:
      return stress[3];
      break;
    case ValueType::SYZ:
      return stress[4];
      break;
    case ValueType::SZX:
      return stress[5];
      break;
    case ValueType::EXX:
      return strain[0];
      break;
    case ValueType::EYY:
      return strain[1];
      break;
    case ValueType::EZZ:
      return strain[2];
      break;
    case ValueType::EXY:
      return strain[3];
      break;
    case ValueType::EYZ:
      return strain[4];
      break;
    case ValueType::EZX:
      return strain[5];
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

  // Получить полный вектор узловых значений
  std::vector<double> getNodesData(ValueType type) const {
    std::vector<double> data;
    data.reserve(crd.size());

    // TODO написать остальные
    if (type == ValueType::SXX)
      data.insert(data.end(), solutionData.SigmaXX.cbegin(),
                  solutionData.SigmaXX.cend());
    else if (type == ValueType::SYY)
      data.insert(data.end(), solutionData.SigmaYY.cbegin(),
                  solutionData.SigmaYY.cend());
    else if (type == ValueType::SZZ)
      data.insert(data.end(), solutionData.SigmaZZ.cbegin(),
                  solutionData.SigmaZZ.cend());
    else if (type == ValueType::SXY)
      data.insert(data.end(), solutionData.SigmaXY.cbegin(),
                  solutionData.SigmaXY.cend());
    else if (type == ValueType::SYZ)
      data.insert(data.end(), solutionData.SigmaYZ.cbegin(),
                  solutionData.SigmaYZ.cend());
    else if (type == ValueType::SZX)
      data.insert(data.end(), solutionData.SigmaZX.cbegin(),
                  solutionData.SigmaZX.cend());
    else if (type == ValueType::EXX)
      data.insert(data.end(), solutionData.DefXX.cbegin(),
                  solutionData.DefXX.cend());
    else if (type == ValueType::EYY)
      data.insert(data.end(), solutionData.DefYY.cbegin(),
                  solutionData.DefYY.cend());
    else if (type == ValueType::EZZ)
      data.insert(data.end(), solutionData.DefZZ.cbegin(),
                  solutionData.DefZZ.cend());
    else if (type == ValueType::EXY)
      data.insert(data.end(), solutionData.DefXY.cbegin(),
                  solutionData.DefXY.cend());
    else if (type == ValueType::EYZ)
      data.insert(data.end(), solutionData.DefYZ.cbegin(),
                  solutionData.DefYZ.cend());
    else if (type == ValueType::EZX)
      data.insert(data.end(), solutionData.DefZX.cbegin(),
                  solutionData.DefZX.cend());

    else {
      throw std::runtime_error("Unimplemented");
    }

    return data;
  }

  // Получить локальный вектор узловых значений
  std::vector<double> getLocalNodesData(size_t elemID, ValueType type) const {
    std::vector<double> data;
    const auto &inds = elems[elemID].inds;
    data.reserve(inds.size());
    auto globData = getNodesData(type);

    for (auto id : inds) {
      data.push_back(globData[id]);
    }

    return data;
  }
};

// Чтение сетки и данных
Mesh readMesh(const std::string &prefix);
// Проверка интерполяций
void checkData(const Mesh &mesh);

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
    else if (vType == ValueType::SYY)
      data[ct] = mesh.solutionData.SigmaYY[ids[ct]];
    else if (vType == ValueType::SZZ)
      data[ct] = mesh.solutionData.SigmaZZ[ids[ct]];
    else if (vType == ValueType::SXY)
      data[ct] = mesh.solutionData.SigmaXY[ids[ct]];
    else if (vType == ValueType::SYZ)
      data[ct] = mesh.solutionData.SigmaYZ[ids[ct]];
    else if (vType == ValueType::SZX)
      data[ct] = mesh.solutionData.SigmaZX[ids[ct]];
    else if (vType == ValueType::EXX)
      data[ct] = mesh.solutionData.DefXX[ids[ct]];
    else if (vType == ValueType::EYY)
      data[ct] = mesh.solutionData.DefYY[ids[ct]];
    else if (vType == ValueType::EZZ)
      data[ct] = mesh.solutionData.DefZZ[ids[ct]];
    else if (vType == ValueType::EXY)
      data[ct] = mesh.solutionData.DefXY[ids[ct]];
    else if (vType == ValueType::EYZ)
      data[ct] = mesh.solutionData.DefYZ[ids[ct]];
    else if (vType == ValueType::EZX)
      data[ct] = mesh.solutionData.DefZX[ids[ct]];

    else
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

typedef dlib::matrix<double, 0, 1> column_vector;
// Функция пробует улучшить узловое значение в узле nodeID
double tryToImproove(const Mesh &mesh, size_t nodeId, ValueType type) {
  // Найти все узлы, соседние с целевым
  indexes indsToVariate;

  // Найти все точки интегрирования, ближайшие к целевому узлу
  auto gpIDs = findNearestGPs(mesh, nodeId);
  // Первый подход - добавить в список вариации все узлы всех сосудних элементов
  for (const auto &d : gpIDs) {
    auto elID = d.first;
    const auto &ids = mesh.elems[elID].inds;
    indsToVariate.insert(indsToVariate.end(), ids.cbegin(), ids.cend());
  }
  auto pos = std::unique(indsToVariate.begin(), indsToVariate.end());
  indsToVariate.erase(pos, indsToVariate.end());

  // Создать функцию невязки

  auto func = [&](const column_vector &deltas) {
    // Вычислить среднеквадратическое отклонение
    // Копия узлового вектора
    auto nodesData = mesh.getNodesData(type);
    for (size_t ct = 0, size = indsToVariate.size(); ct < size; ++ct) {
      nodesData[indsToVariate[ct]] += deltas(ct);
    }

    double err = 0;

    for (auto &data : gpIDs) {
      auto elID = data.first;
      for (auto pID : data.second) {
        const auto &ids = mesh.elems[elID].inds;
        std::vector<double> localData;
        localData.reserve(ids.size());
        for (auto id : ids) {
          localData.push_back(nodesData[id]);
        }
        auto val = calculateValue(mesh.intData[elID][pID].localCoords,
                                  localData, mesh.elems[elID]);
        err += std::pow(val - mesh.intData[elID][pID].getValue(type), 2);
      }
    }

    return err;
  };

  // Варьирование

  // Приращения к узловым значениям
  column_vector deltas;
  deltas.set_size(indsToVariate.size());
  for (auto &d : deltas)
    d = 0;

  std::cout << func(deltas) << std::endl;
  find_min_using_approximate_derivatives(
      dlib::bfgs_search_strategy(), dlib::objective_delta_stop_strategy(1e-7),
      func, deltas, -1);
  for (size_t ct = 0; ct < indsToVariate.size(); ++ct) {
    std::cout << indsToVariate[ct] << "\t" << deltas(ct) << std::endl;
  }
  std::cout << func(deltas) << std::endl;

  auto pas = std::find(indsToVariate.cbegin(), indsToVariate.cend(), nodeId);

  assert(pas != indsToVariate.cend());

  return deltas(pas - indsToVariate.cbegin());
}

std::vector<int> flags(20);

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

  auto mesh = readMesh(workDir);
  auto &crd = mesh.crd;
  auto elems = mesh.elems;
  auto &solutionData = mesh.solutionData;
  const size_t numNodes = crd.size();

  checkData(mesh);

  //Понимаю, что алгоритм ниже неэффективен для plate и будет долго компилироваться

  // Цикл по каждому узлу и минимизация ошибки
  for (size_t nodeID = 0; nodeID < numNodes; ++nodeID) {
    // Вычисляем  приращения  для  текущего  узла
    float deltaSXX = tryToImproove(mesh, nodeID, ValueType::SXX);
    float deltaSYY = tryToImproove(mesh, nodeID, ValueType::SYY);
    float deltaSZZ = tryToImproove(mesh, nodeID, ValueType::SZZ);
    float deltaSXY = tryToImproove(mesh, nodeID, ValueType::SXY);
    float deltaSYZ = tryToImproove(mesh, nodeID, ValueType::SYZ);
    float deltaSZX = tryToImproove(mesh, nodeID, ValueType::SZX);
    float deltaEXX = tryToImproove(mesh, nodeID, ValueType::EXX);
    float deltaEYY = tryToImproove(mesh, nodeID, ValueType::EYY);
    float deltaEZZ = tryToImproove(mesh, nodeID, ValueType::EZZ);
    float deltaEXY = tryToImproove(mesh, nodeID, ValueType::EXY);
    float deltaEYZ = tryToImproove(mesh, nodeID, ValueType::EYZ);
    float deltaEZX = tryToImproove(mesh, nodeID, ValueType::EZX);

    // Обновляем  значения  напрямую
    solutionData.SigmaXX[nodeID] += deltaSXX; 
    solutionData.SigmaYY[nodeID] += deltaSYY;
    solutionData.SigmaZZ[nodeID] += deltaSZZ;
    solutionData.SigmaXY[nodeID] += deltaSXY;
    solutionData.SigmaYZ[nodeID] += deltaSYZ;
    solutionData.SigmaZX[nodeID] += deltaSZX;
    solutionData.DefXX[nodeID] += deltaEXX;
    solutionData.DefYY[nodeID] += deltaEYY;
    solutionData.DefZZ[nodeID] += deltaEZZ;
    solutionData.DefXY[nodeID] += deltaEXY;
    solutionData.DefYZ[nodeID] += deltaEYZ;
    solutionData.DefZX[nodeID] += deltaEZX;
  }

  // Записываем обновленные данные в новый файл
  std::ofstream out(workDir + "solution0002.sba", std::ofstream::binary);
  if (!out.good()) {
    std::cerr << "Ошибка открытия файла solution0002.sba.sba\n";
    return 1;
  }

  // Записываем флаги
  out.write(reinterpret_cast<const char *>(flags.data()),
            sizeof(int) * flags.size());

  // Записываем массивы (если они присутствуют)
  for (size_t i = 0; i < 20; ++i) {
    if (flags[i] == 1) {
      switch (i) {
      case 0:
        out.write(
            reinterpret_cast<const char *>(solutionData.Temperature.data()),
            sizeof(float) * numNodes);
        break;
      case 1:
        out.write(reinterpret_cast<const char *>(solutionData.SigmaXX.data()),
                  sizeof(float) * numNodes);
        break;
      case 2:
        out.write(reinterpret_cast<const char *>(solutionData.SigmaYY.data()),
                  sizeof(float) * numNodes);
        break;
      case 3:
        out.write(reinterpret_cast<const char *>(solutionData.SigmaZZ.data()),
                  sizeof(float) * numNodes);
        break;
      case 4:
        out.write(reinterpret_cast<const char *>(solutionData.SigmaXY.data()),
                  sizeof(float) * numNodes);
        break;
      case 5:
        out.write(reinterpret_cast<const char *>(solutionData.SigmaYZ.data()),
                  sizeof(float) * numNodes);
        break;
      case 6:
        out.write(reinterpret_cast<const char *>(solutionData.SigmaZX.data()),
                  sizeof(float) * numNodes);
        break;
      case 7:
        out.write(reinterpret_cast<const char *>(solutionData.SigmaI.data()),
                  sizeof(float) * numNodes);
        break;
      case 8:
        out.write(reinterpret_cast<const char *>(solutionData.Sigma1.data()),
                  sizeof(float) * numNodes);
        break;
      case 9:
        out.write(reinterpret_cast<const char *>(solutionData.Sigma2.data()),
                  sizeof(float) * numNodes);
        break;
      case 10:
        out.write(reinterpret_cast<const char *>(solutionData.Sigma3.data()),
                  sizeof(float) * numNodes);
        break;
      case 11:
        out.write(reinterpret_cast<const char *>(solutionData.DefXX.data()),
                  sizeof(float) * numNodes);
        break;
      case 12:
        out.write(reinterpret_cast<const char *>(solutionData.DefYY.data()),
                  sizeof(float) * numNodes);
        break;
      case 13:
        out.write(reinterpret_cast<const char *>(solutionData.DefZZ.data()),
                  sizeof(float) * numNodes);
        break;
      case 14:
        out.write(reinterpret_cast<const char *>(solutionData.DefXY.data()),
                  sizeof(float) * numNodes);
        break;
      case 15:
        out.write(reinterpret_cast<const char *>(solutionData.DefYZ.data()),
                  sizeof(float) * numNodes);
        break;
      case 16:
        out.write(reinterpret_cast<const char *>(solutionData.DefZX.data()),
                  sizeof(float) * numNodes);
        break;
      case 17:
        out.write(reinterpret_cast<const char *>(solutionData.DefI.data()),
                  sizeof(float) * numNodes);
        break;
      case 18:
        out.write(reinterpret_cast<const char *>(solutionData.DefPlast.data()),
                  sizeof(float) * numNodes);
        break;
      case 19:
        out.write(reinterpret_cast<const char *>(solutionData.None.data()),
                  sizeof(float) * numNodes);
        break;
      }
    }
  }

  out.close();

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
        case 7:
          solutionData.SigmaI = data;
          break;
        case 8:
          solutionData.Sigma1 = data;
          break;
        case 9:
          solutionData.Sigma2 = data;
          break;
        case 10:
          solutionData.Sigma3 = data;
          break;
        case 11:
          solutionData.DefXX = data;
          break;
        case 12:
          solutionData.DefYY = data;
          break;
        case 13:
          solutionData.DefZZ = data;
          break;
        case 14:
          solutionData.DefXY = data;
          break;
        case 15:
          solutionData.DefYZ = data;
          break;
        case 16:
          solutionData.DefZX = data;
          break;
        case 17:
          solutionData.DefI = data;
          break;
        case 18:
          solutionData.DefPlast = data;
          break;
        case 19:
          solutionData.None = data;
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
