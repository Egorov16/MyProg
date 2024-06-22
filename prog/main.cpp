#include <array>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

int main() {
  setlocale(LC_ALL, "rus");

  // crd
  struct Node {
    float x, y, z;
  };
  std::vector<Node> crd;
  {
    std::ifstream in("crd.sba",
                     std::ifstream::binary);
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
    std::ifstream in("ind.sba",
                     std::ifstream::binary);
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

  return 0;
}