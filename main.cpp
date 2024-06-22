#include <array>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>


// Функция формы
float shapeFunction(float localCoord, float xi) {
  // Линейная функция формы
  return 0.5 * (1.0 + xi * localCoord);
}

// Функция, которая вычисляет значение в заданных локальных координатах
float calculateValue(float localX, float localY, float localZ,
                     const std::vector<float> &nodalValues) {
  // Проверяем размер массива узловых значений
  if (nodalValues.size() != 8) {
    std::cerr << "Ошибка: Массив узловых значений должен иметь 8 элементов!"
              << std::endl;
    return 0.0;
  }

  // Вычисляем значения функций формы
  float N1 = shapeFunction(localX, -1.0) * shapeFunction(localY, -1.0) *
             shapeFunction(localZ, -1.0);
  std::cout << "N1 = " << N1 << std::endl;
  float N2 = shapeFunction(localX, 1.0) * shapeFunction(localY, -1.0) *
             shapeFunction(localZ, -1.0);
  std::cout << "N2 = " << N2 << std::endl;
  float N3 = shapeFunction(localX, 1.0) * shapeFunction(localY, 1.0) *
             shapeFunction(localZ, -1.0);
  std::cout << "N3 = " << N3 << std::endl;
  float N4 = shapeFunction(localX, -1.0) * shapeFunction(localY, 1.0) *
             shapeFunction(localZ, -1.0);
  std::cout << "N4 = " << N4 << std::endl;
  float N5 = shapeFunction(localX, -1.0) * shapeFunction(localY, -1.0) *
             shapeFunction(localZ, 1.0);
  std::cout << "N5 = " << N5 << std::endl;
  float N6 = shapeFunction(localX, 1.0) * shapeFunction(localY, -1.0) *
             shapeFunction(localZ, 1.0);
  std::cout << "N6 = " << N6 << std::endl;
  float N7 = shapeFunction(localX, 1.0) * shapeFunction(localY, 1.0) *
             shapeFunction(localZ, 1.0);
  std::cout << "N7 = " << N7 << std::endl;
  float N8 = shapeFunction(localX, -1.0) * shapeFunction(localY, 1.0) *
             shapeFunction(localZ, 1.0);
  std::cout << "N8 = " << N8 << std::endl;

  // Вычисляем значение в точке
  float value = N1 * nodalValues[0] + N2 * nodalValues[1] +
                N3 * nodalValues[2] + N4 * nodalValues[3] +
                N5 * nodalValues[4] + N6 * nodalValues[5] +
                N7 * nodalValues[6] + N8 * nodalValues[7];

  return value;
}

int main() {
  setlocale(LC_ALL, "rus");

  // crd
  struct Node {
    float x, y, z;
  };

  std::vector<Node> crd;
  {
    std::ifstream in("crd.sba", std::ifstream::binary);
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

  // Находим максимальную и минимальную координаты по каждой оси
  float maxX = std::numeric_limits<float>::min();
  float minX = std::numeric_limits<float>::max();
  float maxY = std::numeric_limits<float>::min();
  float minY = std::numeric_limits<float>::max();
  float maxZ = std::numeric_limits<float>::min();
  float minZ = std::numeric_limits<float>::max();

  for (const auto &node : crd) {
    maxX = std::max(maxX, node.x);
    minX = std::min(minX, node.x);
    maxY = std::max(maxY, node.y);
    minY = std::min(minY, node.y);
    maxZ = std::max(maxZ, node.z);
    minZ = std::min(minZ, node.z);
  }

  // Вывод результата
  std::cout << "Максимальная координата X: " << maxX << std::endl;
  std::cout << "Минимальная координата X: " << minX << std::endl;
  std::cout << "Максимальная координата Y: " << maxY << std::endl;
  std::cout << "Минимальная координата Y: " << minY << std::endl;
  std::cout << "Максимальная координата Z: " << maxZ << std::endl;
  std::cout << "Минимальная координата Z: " << minZ << std::endl;

  using indexes = std::vector<int>;

  std::vector<std::pair<int, indexes>> elems;

  // ind
  {
    std::ifstream in("ind.sba", std::ifstream::binary);
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
            t.first, ind)); // Заполняем вектор elems парами номер КЭ - список
                            // узлов, которые в него включены
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

  // Выводим полученный вектор
  std::cout << "Вектор nodeElements:" << std::endl;
  for (const auto &item : nodeElements) { // Перебираем пары в nodeElements
    std::cout << "  " << item.first << "  -  "; // Выводим ключ (узел)
    for (const auto &element : item.second) { // Перебираем элементы узла
      std::cout << "  " << element << "  ";
    }
    std::cout << std::endl;
  }

  // Пример использования функции calculateValue
  std::vector<float> nodalValues = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  float localX = 0.5;
  float localY = 0.25;
  float localZ = -0.75;

  float value = calculateValue(localX, localY, localZ, nodalValues);

  std::cout << "Значение в точке (" << localX << ", " << localY << ", "
            << localZ << "): " << value << std::endl;

  return 0;
}