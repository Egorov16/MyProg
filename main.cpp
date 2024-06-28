#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

// Структура для узла
struct Node {
  double x;
  double y;
  double z;
};

// Структура для конечного элемента
struct Element {
  int node1;
  int node2;
  int node3;
  int node4;
  int node5;
  int node6;
  int node7;
  int node8;

  // Точки интегрирования (по умолчанию, для куба, центр в (0, 0, 0))
  std::vector<Node> integrationPoints = {
      {-0.57735026918962573, -0.57735026918962573,
       -0.57735026918962573}, // Точка 0
      {-0.57735026918962573, -0.57735026918962573,
       0.57735026918962573}, // Точка 1
      {-0.57735026918962573, 0.57735026918962573,
       -0.57735026918962573}, // Точка 2
      {-0.57735026918962573, 0.57735026918962573,
       0.57735026918962573}, // Точка 3
      {0.57735026918962573, -0.57735026918962573,
       -0.57735026918962573}, // Точка 4
      {0.57735026918962573, -0.57735026918962573,
       0.57735026918962573}, // Точка 5
      {0.57735026918962573, 0.57735026918962573,
       -0.57735026918962573}, // Точка 6
      {0.57735026918962573, 0.57735026918962573,
       0.57735026918962573}, // Точка 7
  };
};

// Функция формы
double shapeFunction(double localCoord, double xi) {
  // Линейная функция формы
  return 0.5 * (1.0 + xi * localCoord);
}

double calculateValue(double localX, double localY, double localZ,
                      const std::vector<double> &nodalValues,
                      const Element &element) {

  // Вычисляем значения функций формы
  float N1 = shapeFunction(localX, -1.0) * shapeFunction(localY, -1.0) *
             shapeFunction(localZ, -1.0);
  float N2 = shapeFunction(localX, 1.0) * shapeFunction(localY, -1.0) *
             shapeFunction(localZ, -1.0);
  float N3 = shapeFunction(localX, 1.0) * shapeFunction(localY, 1.0) *
             shapeFunction(localZ, -1.0);
  float N4 = shapeFunction(localX, -1.0) * shapeFunction(localY, 1.0) *
             shapeFunction(localZ, -1.0);
  float N5 = shapeFunction(localX, -1.0) * shapeFunction(localY, -1.0) *
             shapeFunction(localZ, 1.0);
  float N6 = shapeFunction(localX, 1.0) * shapeFunction(localY, -1.0) *
             shapeFunction(localZ, 1.0);
  float N7 = shapeFunction(localX, 1.0) * shapeFunction(localY, 1.0) *
             shapeFunction(localZ, 1.0);
  float N8 = shapeFunction(localX, -1.0) * shapeFunction(localY, 1.0) *
             shapeFunction(localZ, 1.0);

  // Вычисляем значение в точке
  float value = N1 * nodalValues[0] + // Индекс 0 соответствует node1 в Element
                N2 * nodalValues[1] + // Индекс 1 соответствует node2 в Element
                N3 * nodalValues[2] + // Индекс 2 соответствует node3 в Element
                N4 * nodalValues[3] + // Индекс 3 соответствует node4 в Element
                N5 * nodalValues[4] + // Индекс 4 соответствует node5 в Element
                N6 * nodalValues[5] + // Индекс 5 соответствует node6 в Element
                N7 * nodalValues[6] + // Индекс 6 соответствует node7 в Element
                N8 * nodalValues[7]; // Индекс 7 соответствует node8 в Element

  return value;
};

// Функция для расчета расстояния между двумя точками
double distance(const Node &point1, const Node &point2) {
  return std::sqrt(std::pow(point1.x - point2.x, 2) +
                   std::pow(point1.y - point2.y, 2) +
                   std::pow(point1.z - point2.z, 2));
}

// Функция для поиска ближайшей точки интегрирования
int findClosestIntegrationPoint(const Node &targetNode,
                                const Element &element) {
  double minDistance = std::numeric_limits<double>::max();
  int closestPointIndex = 0;

  for (size_t i = 0; i < element.integrationPoints.size(); ++i) {
    double currentDistance = distance(targetNode, element.integrationPoints[i]);
    if (currentDistance < minDistance) {
      minDistance = currentDistance;
      closestPointIndex = i;
    }
  }

  return closestPointIndex;
}

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

// Структура для хранения данных в точках интегрирования
struct IntegrationPointData {
  int elementNumber;
  int gaussPointNumber;
  double localCoords[3]; // Массив для хранения локальных координат
  double globalCoords[3]; // Массив для хранения локальных координат
  double stress[6]; // Массив для хранения напряжений
};

struct GaussPointData {
  int elementNumber;
  int gaussPointNumber;
  double localCoords[3]; // Массив для хранения локальных координат
  double globalCoords[3]; // Массив для хранения локальных координат
  double stress[6]; // Массив для хранения напряжений
  double strain[6]; // Массив для хранения деформаций
};

// Функция для вычисления невязки
std::vector<std::vector<double>>
calculateDiscrepancy(const std::vector<IntegrationPointData> &interpolatedData,
                     const std::vector<GaussPointData> &fileData) {

  std::vector<std::vector<double>> discrepancy;
  // Внешний вектор для хранения невязок по элементам
  // Внутренний вектор для хранения невязок по компонентам напряжения в элементе

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

int main() {
  setlocale(LC_ALL, "rus");
  // crd

  std::vector<Node> crd;
  {
    std::ifstream in("C:/git/MyProg/plate/crd.sba", std::ifstream::binary);
    if (!in.good())
      throw std::runtime_error("Crd");
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

  const size_t numNodes = crd.size();

  using indexes = std::vector<int>;

  std::vector<std::pair<int, indexes>> elems;

  // ind
  {
    std::ifstream in("C:/git/MyProg/plate/ind.sba", std::ifstream::binary);
    if (!in.good())
      throw std::runtime_error("Ind");
    int numTypes;
    in.read((char *)&numTypes, sizeof(int));
    std::vector<std::pair<int, int>> types(numTypes);
    for (auto &t : types) {
      in.read((char *)&(t.first), sizeof(int));
      in.read((char *)&(t.second), sizeof(int));
    }

    int indexElem = 0;
    for (auto t : types) {
      int numIndexes;
      in.read((char *)&numIndexes, sizeof(int));

      for (int i = 0; i < t.second; ++i) {
        indexes ind;
        ind.resize(numIndexes / t.second);
        in.read((char *)ind.data(), sizeof(int) * ind.size());
        elems.push_back(std::make_pair(
            t.first,
            ind)); // Заполняем вектор elems парами Тип элемента - список узлов
      }
    }
  }

  SolutionData solutionData; // Создаем объект SolutionData

  {
    std::ifstream in("C:/git/MyProg/plate/solution0001.sba",
                     std::ifstream::binary);

    if (!in.good()) {
      std::cerr << "Ошибка открытия файла solution0001.sba\n";
      return 1;
    }

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

  // Создаем мапу для хранения информации о узлах и принадлежащих им конечных
  // элементах
  std::map<int, std::vector<int>> nodeElements;

  // Заполняем вектор nodeElements
  for (size_t i = 0; i < elems.size(); ++i) {
    for (const auto &node : elems[i].second) {
      nodeElements[node].push_back(i + 1); // Добавляем элемент в вектор узла
    }
  }

  // Вектор для хранения данных в точках интегрирования
  std::vector<IntegrationPointData> integrationPointData;

  for (int targetNode = 1; targetNode <= crd.size(); ++targetNode) {

    for (int elementIndex : nodeElements[targetNode]) {

      // Получаем номера узлов из elems:
      auto nodeIndexes = elems[elementIndex - 1].second;

      // Создаем объект Element:
      Element elem = {nodeIndexes[0], nodeIndexes[1], nodeIndexes[2],
                      nodeIndexes[3], nodeIndexes[4], nodeIndexes[5],
                      nodeIndexes[6], nodeIndexes[7]};

      // Индекс ближайшей точки интегрирования
      int locNumNode = -1; // Индекс числа, по умолчанию -1 (не найдено)
      for (int i = 0; i < 8; ++i) {
        if (nodeIndexes[i] == targetNode) {
          locNumNode = i;
          break; // Прерываем цикл, как только найдем число
        }
      }
      int closestPointIndex;
      switch (locNumNode) {
      case 0:
        closestPointIndex = 0;
        break;
      case 1:
        closestPointIndex = 4;
        break;
      case 2:
        closestPointIndex = 6;
        break;
      case 3:
        closestPointIndex = 2;
        break;
      case 4:
        closestPointIndex = 1;
        break;
      case 5:
        closestPointIndex = 5;
        break;
      case 6:
        closestPointIndex = 7;
        break;
      case 7:
        closestPointIndex = 3;
        break;
      }

      // Получаем координаты точки интегрирования:
      double localX = elem.integrationPoints[closestPointIndex].x;
      double localY = elem.integrationPoints[closestPointIndex].y;
      double localZ = elem.integrationPoints[closestPointIndex].z;

      // Создаём вектора узловых значений координат
      std::vector<double> coordsNodeX;
      std::vector<double> coordsNodeY;
      std::vector<double> coordsNodeZ;

      for (int i = 0; i < 8; i++) {
        int nodeIndex = nodeIndexes[i];
        coordsNodeX.push_back(crd[nodeIndex - 1].x);
        coordsNodeY.push_back(crd[nodeIndex - 1].y);
        coordsNodeZ.push_back(crd[nodeIndex - 1].z);
      }

      // Глобальные координаты точки интегрирования
      double globalX =
          calculateValue(localX, localY, localZ, coordsNodeX, elem);
      double globalY =
          calculateValue(localX, localY, localZ, coordsNodeY, elem);
      double globalZ =
          calculateValue(localX, localY, localZ, coordsNodeZ, elem);

      // Создаем массив для хранения интерполированных значений
      double stress[6];

      // Итерация по компонентам напряжения и интерполяция
      for (int i = 0; i < 6; ++i) {
        // Получаем значения в узлах элемента для текущего компонента
        std::vector<double> nodalValues;
        for (int node : nodeIndexes) {
          // Используем правильный доступ к компонентам solutionData
          switch (i) {
          case 0:
            nodalValues.push_back(solutionData.SigmaXX[node - 1]);
            break;
          case 1:
            nodalValues.push_back(solutionData.SigmaYY[node - 1]);
            break;
          case 2:
            nodalValues.push_back(solutionData.SigmaZZ[node - 1]);
            break;
          case 3:
            nodalValues.push_back(solutionData.SigmaXY[node - 1]);
            break;
          case 4:
            nodalValues.push_back(solutionData.SigmaYZ[node - 1]);
            break;
          case 5:
            nodalValues.push_back(solutionData.SigmaZX[node - 1]);
            break;
          }
        }

        // Интерполяция значения в ближайшую точку интегрирования
        stress[i] = calculateValue(localX, localY, localZ, nodalValues, elem);
      }
      // Добавляем данные в вектор integrationPointData
      IntegrationPointData data;
      data.elementNumber = elementIndex;
      data.gaussPointNumber = closestPointIndex;
      data.localCoords[0] = localX;
      data.localCoords[1] = localY;
      data.localCoords[2] = localZ;
      data.globalCoords[0] = globalX;
      data.globalCoords[1] = globalY;
      data.globalCoords[2] = globalZ;

      // Копируем интерполированные значения в data.stress
      std::copy(std::begin(stress), std::end(stress), std::begin(data.stress));

      integrationPointData.push_back(data);
    }
  }

  std::sort(integrationPointData.begin(), integrationPointData.end(),
            [](const IntegrationPointData &a, const IntegrationPointData &b) {
              // Сначала сортируем по номеру элемента
              if (a.elementNumber != b.elementNumber) {
                return a.elementNumber < b.elementNumber;
              } else {
                // Если элементы одинаковые, сортируем по номеру точки
                // интегрирования
                return a.gaussPointNumber < b.gaussPointNumber;
              }
            });

  // int

  std::vector<GaussPointData> data;

  std::ifstream file("C:/git/MyProg/plate/int.dat", std::ifstream::binary);

  if (!file.is_open()) {
    std::cerr << "Ошибка открытия файла int.dat\n";
    return 1;
  }

  while (!file.eof()) {
    GaussPointData point;
    file.read(reinterpret_cast<char *>(&point.elementNumber), sizeof(double));
    file.read(reinterpret_cast<char *>(&point.gaussPointNumber),
              sizeof(double));
    // Считываем локальные координаты в массив
    file.read((char *)(point.localCoords), sizeof(double) * 3);
    // Считываем глобальные координаты в массив
    file.read((char *)(point.globalCoords), sizeof(double) * 3);
    // Считываем напряжения в массив
    file.read((char *)(point.stress), sizeof(double) * 6);
    // Считываем деформации в массив
    file.read((char *)(point.strain), sizeof(double) * 6);

    // Увеличиваем номер КЭ на 1
    point.elementNumber++;

    data.push_back(point);
  }

  file.close();

  data.pop_back();

  // Вычисление невязки
  std::vector<std::vector<double>> discrepancy =
      calculateDiscrepancy(integrationPointData, data);

  // Поиск максимального значения невязки
  double maxDiscrepancy = 0.0;
  for (const auto &elementDiscrepancy : discrepancy) {
    for (const auto &value : elementDiscrepancy) {
      maxDiscrepancy = std::max(maxDiscrepancy, value);
    }
  }

  std::cout << "\nМаксимальная невязка: " << maxDiscrepancy << "%" << std::endl;

  return 0;
}