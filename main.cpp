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

const std::array<Node, 4> TetrahedronNodes{
    Node{0, 0, 0}, // Узел 1
    Node{1, 0, 0}, // Узел 2
    Node{0, 1, 0}, // Узел 3
    Node{0, 0, 1}, // Узел 4
};

const std::array<Node, 10> quadTetrahedronNodes{
    Node{0, 0, 0},     // Узел 1
    Node{1, 0, 0},     // Узел 2
    Node{0, 1, 0},     // Узел 3
    Node{0, 0, 1},     // Узел 4
    Node{0.5, 0, 0},   // Узел 5
    Node{0, 0.5, 0},   // Узел 6
    Node{0, 0, 0.5},   // Узел 7
    Node{0.5, 0.5, 0}, // Узел 8
    Node{0.5, 0, 0.5}, // Узел 9
    Node{0, 0.5, 0.5}, // Узел 10
};

using indexes = std::vector<size_t>;
using Crd = std::vector<Node>;
enum class Type { hexa8 = 0, hexa20 = 1, tet4 = 2, tet10 = 3 };
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
};

enum class ValueType { CrdX, CrdY, CrdZ, SXX };

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

double quadcalculateValue(Node localCrd, const std::vector<double> &nodalValues,
                          const Element &element) {

  assert(element.inds.size() == nodalValues.size());
  // TODO написать другие типы
  assert(element.type == Type::hexa20);

  std::vector<double> h(20);

  // Вершины
  h[0] = 0.125 * (1.0 - localCrd.x) * (1.0 - localCrd.y) * (1.0 - localCrd.z) *
         (-localCrd.x - localCrd.y - localCrd.z - 2);
  h[1] = 0.125 * (1.0 + localCrd.x) * (1.0 - localCrd.y) * (1.0 - localCrd.z) *
         (localCrd.x - localCrd.y - localCrd.z - 2);
  h[2] = 0.125 * (1.0 + localCrd.x) * (1.0 + localCrd.y) * (1.0 - localCrd.z) *
         (localCrd.x + localCrd.y - localCrd.z - 2);
  h[3] = 0.125 * (1.0 - localCrd.x) * (1.0 + localCrd.y) * (1.0 - localCrd.z) *
         (-localCrd.x + localCrd.y - localCrd.z - 2);
  h[4] = 0.125 * (1.0 - localCrd.x) * (1.0 - localCrd.y) * (1.0 + localCrd.z) *
         (-localCrd.x - localCrd.y + localCrd.z - 2);
  h[5] = 0.125 * (1.0 + localCrd.x) * (1.0 - localCrd.y) * (1.0 + localCrd.z) *
         (localCrd.x - localCrd.y + localCrd.z - 2);
  h[6] = 0.125 * (1.0 + localCrd.x) * (1.0 + localCrd.y) * (1.0 + localCrd.z) *
         (localCrd.x + localCrd.y + localCrd.z - 2);
  h[7] = 0.125 * (1.0 - localCrd.x) * (1.0 + localCrd.y) * (1.0 + localCrd.z) *
         (-localCrd.x + localCrd.y + localCrd.z - 2);

  // Ребра
  h[8] = 0.25 * (1 - localCrd.x * localCrd.x) * (1 - localCrd.y) *
         (1 + localCrd.z);
  h[9] = 0.25 * (1 + localCrd.x) * (1 - localCrd.y * localCrd.y) *
         (1 + localCrd.z);
  h[10] = 0.25 * (1 - localCrd.x * localCrd.x) * (1 + localCrd.y) *
          (1 + localCrd.z);
  h[11] = 0.25 * (1 - localCrd.x) * (1 - localCrd.y * localCrd.y) *
          (1 + localCrd.z);
  h[12] = 0.25 * (1 - localCrd.x * localCrd.x) * (1 - localCrd.y) *
          (1 - localCrd.z);
  h[13] = 0.25 * (1 + localCrd.x) * (1 - localCrd.y * localCrd.y) *
          (1 - localCrd.z);
  h[14] = 0.25 * (1 - localCrd.x * localCrd.x) * (1 + localCrd.y) *
          (1 - localCrd.z);
  h[15] = 0.25 * (1 - localCrd.x) * (1 - localCrd.y * localCrd.y) *
          (1 - localCrd.z);
  h[16] = 0.25 * (1 + localCrd.x) * (1 + localCrd.y) *
          (1 - localCrd.z * localCrd.z);
  h[17] = 0.25 * (1 - localCrd.x) * (1 + localCrd.y) *
          (1 - localCrd.z * localCrd.z);
  h[18] = 0.25 * (1 - localCrd.x) * (1 - localCrd.y) *
          (1 - localCrd.z * localCrd.z);
  h[19] = 0.25 * (1 + localCrd.x) * (1 - localCrd.y) *
          (1 - localCrd.z * localCrd.z);

  // Вычисляем значение в точке
  double value = 0;
  for (size_t ct = 0, size = h.size(); ct < size; ++ct)
    value += h[ct] * nodalValues[ct];
  return value;
};

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

double TetrahedroncalculateValue(Node localCrd,
                                 const std::vector<double> &nodalValues,
                                 const Element &element) {

  assert(element.inds.size() == nodalValues.size());
  // TODO написать другие типы
  assert(element.type == Type::tet4);

  std::vector<double> h(8);
  {
    // Координаты узлов тетраэдра
    double x1 = TetrahedronNodes[0].x;
    double y1 = TetrahedronNodes[0].y;
    double z1 = TetrahedronNodes[0].z;
    double x2 = TetrahedronNodes[1].x;
    double y2 = TetrahedronNodes[1].y;
    double z2 = TetrahedronNodes[1].z;
    double x3 = TetrahedronNodes[2].x;
    double y3 = TetrahedronNodes[2].y;
    double z3 = TetrahedronNodes[2].z;
    double x4 = TetrahedronNodes[3].x;
    double y4 = TetrahedronNodes[3].y;
    double z4 = TetrahedronNodes[3].z;

    // Координаты точки внутри тетраэдра
    double x = localCrd.x;
    double y = localCrd.y;
    double z = localCrd.z;

    // Вычисление функций формы
    std::vector<double> N(4);

    h[0] = ((x2 * y3 * z4 - x2 * z3 * y4 + x3 * y4 * z2 - x3 * y2 * z4 +
             x4 * y2 * z3 - x4 * y3 * z2) +
            (y2 * z4 - y3 * z4 + z3 * y4 - y4 * z2 - y2 * z3 + y3 * z2) * x +
            (x3 * z4 - x2 * z4 + x4 * z2 + x2 * z3 - x3 * z2 - x4 * z3) * y +
            (-x4 * y2 + x4 * y3 - x2 * y3 - x3 * y4 + x2 * y4 + x3 * y2) * z) /
           (x1 * (y2 * z4 - y3 * z4 + z3 * y4 - y4 * z2 - y2 * z3 + y3 * z2) +
            x2 * (-z3 * y4 + y3 * z4 - y3 * z1 + y1 * z3 + y4 * z1 - y1 * z4) +
            x3 * (y1 * z4 + y4 * z2 - y2 * z4 - y4 * z1 - y1 * z2 + y2 * z1) +
            x4 * (y2 * z3 + y3 * z1 + y1 * z2 - y2 * z1 - y3 * z2 - y1 * z3));

    h[1] = ((-x1 * y3 * z4 + x1 * z3 * y4 - x3 * y4 * z1 + x3 * y1 * z4 -
             x4 * y1 * z3 + x4 * y3 * z1) +
            (-z3 * y4 + y3 * z4 - y3 * z1 + y1 * z3 + y4 * z1 - y1 * z4) * x +
            (x4 * z3 + x1 * z4 - x3 * z4 - x4 * z1 - x1 * z3 + x3 * z1) * y +
            (-x4 * y3 + x3 * y4 + x4 * y1 - x1 * y4 - x3 * y1 + x1 * y3) * z) /
           (x1 * (y2 * z4 - y3 * z4 + z3 * y4 - y4 * z2 - y2 * z3 + y3 * z2) +
            x2 * (-z3 * y4 + y3 * z4 - y3 * z1 + y1 * z3 + y4 * z1 - y1 * z4) +
            x3 * (y1 * z4 + y4 * z2 - y2 * z4 - y4 * z1 - y1 * z2 + y2 * z1) +
            x4 * (y2 * z3 + y3 * z1 + y1 * z2 - y2 * z1 - y3 * z2 - y1 * z3));

    h[2] = ((x2 * y4 * z1 - x4 * y2 * z1 + x1 * y2 * z4 - x1 * y4 * z2 -
             x2 * y1 * z4 + x4 * y1 * z2) +
            (y1 * z4 + y4 * z2 - y2 * z4 - y4 * z1 - y1 * z2 + y2 * z1) * x +
            (-x1 * z4 + x4 * z1 - x2 * z1 + x2 * z4 - x4 * z2 + x1 * z2) * y +
            (x4 * y2 - x2 * y4 + x1 * y4 + x2 * y1 - x1 * y2 - x4 * y1) * z) /
           (x1 * (y2 * z4 - y3 * z4 + z3 * y4 - y4 * z2 - y2 * z3 + y3 * z2) +
            x2 * (-z3 * y4 + y3 * z4 - y3 * z1 + y1 * z3 + y4 * z1 - y1 * z4) +
            x3 * (y1 * z4 + y4 * z2 - y2 * z4 - y4 * z1 - y1 * z2 + y2 * z1) +
            x4 * (y2 * z3 + y3 * z1 + y1 * z2 - y2 * z1 - y3 * z2 - y1 * z3));

    h[3] = ((x1 * y3 * z2 - x3 * y1 * z2 - x2 * y3 * z1 + x3 * y2 * z1 +
             x2 * y1 * z3 - x1 * y2 * z3) +
            (y2 * z3 + y3 * z1 + y1 * z2 - y2 * z1 - y3 * z2 - y1 * z3) * x +
            (-x2 * z3 - x1 * z2 + x1 * z3 - x3 * z1 + x2 * z1 + x3 * z2) * y +
            (x1 * y2 - x2 * y1 - x3 * y2 + x2 * y3 + x3 * y1 - x1 * y3) * z) /
           (x1 * (y2 * z4 - y3 * z4 + z3 * y4 - y4 * z2 - y2 * z3 + y3 * z2) +
            x2 * (-z3 * y4 + y3 * z4 - y3 * z1 + y1 * z3 + y4 * z1 - y1 * z4) +
            x3 * (y1 * z4 + y4 * z2 - y2 * z4 - y4 * z1 - y1 * z2 + y2 * z1) +
            x4 * (y2 * z3 + y3 * z1 + y1 * z2 - y2 * z1 - y3 * z2 - y1 * z3));
  }

  // Вычисляем значение в точке
  double value = 0;
  for (size_t ct = 0, size = h.size(); ct < size; ++ct)
    value += h[ct] * nodalValues[ct];
  return value;
};

double quadTetrahedroncalculateValue(Node localCrd,
                                 const std::vector<double> &nodalValues,
                                 const Element &element) {

  assert(element.inds.size() == nodalValues.size());
  // TODO написать другие типы
  assert(element.type == Type::tet10);

  std::vector<double> h(10);
  {
    // Координаты узлов тетраэдра
    double x1 = TetrahedronNodes[0].x;
    double y1 = TetrahedronNodes[0].y;
    double z1 = TetrahedronNodes[0].z;
    double x2 = TetrahedronNodes[1].x;
    double y2 = TetrahedronNodes[1].y;
    double z2 = TetrahedronNodes[1].z;
    double x3 = TetrahedronNodes[2].x;
    double y3 = TetrahedronNodes[2].y;
    double z3 = TetrahedronNodes[2].z;
    double x4 = TetrahedronNodes[3].x;
    double y4 = TetrahedronNodes[3].y;
    double z4 = TetrahedronNodes[3].z;

    // Координаты точки внутри тетраэдра
    double x = localCrd.x;
    double y = localCrd.y;
    double z = localCrd.z;

    double L1 =
        ((x2 * y3 * z4 - x2 * z3 * y4 + x3 * y4 * z2 - x3 * y2 * z4 +
          x4 * y2 * z3 - x4 * y3 * z2) +
         (y2 * z4 - y3 * z4 + z3 * y4 - y4 * z2 - y2 * z3 + y3 * z2) * x +
         (x3 * z4 - x2 * z4 + x4 * z2 + x2 * z3 - x3 * z2 - x4 * z3) * y +
         (-x4 * y2 + x4 * y3 - x2 * y3 - x3 * y4 + x2 * y4 + x3 * y2) * z) /
        (x1 * (y2 * z4 - y3 * z4 + z3 * y4 - y4 * z2 - y2 * z3 + y3 * z2) +
         x2 * (-z3 * y4 + y3 * z4 - y3 * z1 + y1 * z3 + y4 * z1 - y1 * z4) +
         x3 * (y1 * z4 + y4 * z2 - y2 * z4 - y4 * z1 - y1 * z2 + y2 * z1) +
         x4 * (y2 * z3 + y3 * z1 + y1 * z2 - y2 * z1 - y3 * z2 - y1 * z3));

    double L2 =
        ((-x1 * y3 * z4 + x1 * z3 * y4 - x3 * y4 * z1 + x3 * y1 * z4 -
          x4 * y1 * z3 + x4 * y3 * z1) +
         (-z3 * y4 + y3 * z4 - y3 * z1 + y1 * z3 + y4 * z1 - y1 * z4) * x +
         (x4 * z3 + x1 * z4 - x3 * z4 - x4 * z1 - x1 * z3 + x3 * z1) * y +
         (-x4 * y3 + x3 * y4 + x4 * y1 - x1 * y4 - x3 * y1 + x1 * y3) * z) /
        (x1 * (y2 * z4 - y3 * z4 + z3 * y4 - y4 * z2 - y2 * z3 + y3 * z2) +
         x2 * (-z3 * y4 + y3 * z4 - y3 * z1 + y1 * z3 + y4 * z1 - y1 * z4) +
         x3 * (y1 * z4 + y4 * z2 - y2 * z4 - y4 * z1 - y1 * z2 + y2 * z1) +
         x4 * (y2 * z3 + y3 * z1 + y1 * z2 - y2 * z1 - y3 * z2 - y1 * z3));

    double L3 =
        ((x2 * y4 * z1 - x4 * y2 * z1 + x1 * y2 * z4 - x1 * y4 * z2 -
          x2 * y1 * z4 + x4 * y1 * z2) +
         (y1 * z4 + y4 * z2 - y2 * z4 - y4 * z1 - y1 * z2 + y2 * z1) * x +
         (-x1 * z4 + x4 * z1 - x2 * z1 + x2 * z4 - x4 * z2 + x1 * z2) * y +
         (x4 * y2 - x2 * y4 + x1 * y4 + x2 * y1 - x1 * y2 - x4 * y1) * z) /
        (x1 * (y2 * z4 - y3 * z4 + z3 * y4 - y4 * z2 - y2 * z3 + y3 * z2) +
         x2 * (-z3 * y4 + y3 * z4 - y3 * z1 + y1 * z3 + y4 * z1 - y1 * z4) +
         x3 * (y1 * z4 + y4 * z2 - y2 * z4 - y4 * z1 - y1 * z2 + y2 * z1) +
         x4 * (y2 * z3 + y3 * z1 + y1 * z2 - y2 * z1 - y3 * z2 - y1 * z3));

    double L4 =
        ((x1 * y3 * z2 - x3 * y1 * z2 - x2 * y3 * z1 + x3 * y2 * z1 +
          x2 * y1 * z3 - x1 * y2 * z3) +
         (y2 * z3 + y3 * z1 + y1 * z2 - y2 * z1 - y3 * z2 - y1 * z3) * x +
         (-x2 * z3 - x1 * z2 + x1 * z3 - x3 * z1 + x2 * z1 + x3 * z2) * y +
         (x1 * y2 - x2 * y1 - x3 * y2 + x2 * y3 + x3 * y1 - x1 * y3) * z) /
        (x1 * (y2 * z4 - y3 * z4 + z3 * y4 - y4 * z2 - y2 * z3 + y3 * z2) +
         x2 * (-z3 * y4 + y3 * z4 - y3 * z1 + y1 * z3 + y4 * z1 - y1 * z4) +
         x3 * (y1 * z4 + y4 * z2 - y2 * z4 - y4 * z1 - y1 * z2 + y2 * z1) +
         x4 * (y2 * z3 + y3 * z1 + y1 * z2 - y2 * z1 - y3 * z2 - y1 * z3));

    // Формулы для функций формы квадратичного тетраэдра

    h[0] = (2 * L1 - 1) * L1; // Вершина 0
    h[1] = (2 * L2 - 1) * L2; // Вершина 1
    h[2] = (2 * L3 - 1) * L3; // Вершина 2
    h[3] = (2 * L4 - 1) * L4; // Вершина 3
    h[4] = 4 * L1 * L2;       // Средняя точка на ребре 0-1
    h[5] = 4 * L1 * L3;       // Средняя точка на ребре 0-2
    h[6] = 4 * L1 * L4;       // Средняя точка на ребре 0-3
    h[7] = 4 * L2 * L3;       // Средняя точка на ребре 1-2
    h[8] = 4 * L2 * L4;       // Средняя точка на ребре 1-3
    h[9] = 4 * L3 * L4;       // Средняя точка на ребре 2-3
  }

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

typedef dlib::matrix<double, 0, 1> column_vector;

double norm(const column_vector &v) { return sqrt(dot(v, v)); }

typedef dlib::matrix<double, 0, 1> column_vector;

double dot(const column_vector &a, const column_vector &b) {
  double sum = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    sum += a(i) * b(i);
  }
  return sum;
}

typedef dlib::matrix<double, 0, 1> column_vector;

// Функция оптимизации методом сопряжённого градиента

column_vector conjugateGradientDescent(
    const std::function<double(const column_vector &)> &func,
    const column_vector &initial_deltas, double tolerance, int max_iterations) {
  column_vector x = initial_deltas;
  column_vector grad(x.size());
  column_vector d(x.size());
  column_vector grad_prev(
      x.size()); // Добавили переменную для предыдущего градиента

  // Вычисление начального градиента
  double h = 1e-6;
  for (size_t i = 0; i < x.size(); ++i) {
    column_vector x_plus_h = x;
    x_plus_h(i) += h;
    grad(i) = (func(x_plus_h) - func(x)) / h;
  }

  d = -grad;
  grad_prev = grad; // Инициализировали grad_prev

  int iteration = 0;
  while (iteration < max_iterations) {
    // Вычисление шага alpha
    double alpha = 1.0;
    double f_x = func(x);

    while (func(x + alpha * d) > f_x && alpha > 1e-10) {
      alpha /= 2.0;
    }

    // Обновление значения x
    x = x + alpha * d;

    // Вычисление нового градиента
    for (size_t i = 0; i < x.size(); ++i) {
      column_vector x_plus_h = x;
      x_plus_h(i) += h;
      grad(i) = (func(x_plus_h) - func(x)) / h;
    }

    // Обновление направления спуска
    double beta =
        dot(grad, grad) / dot(grad_prev, grad_prev); // Исправленная формула
    d = -grad + beta * d;

    // Сохранение предыдущего градиента
    grad_prev = grad;

    // Проверка условия останова
    if (norm(grad) <= tolerance) {
      break;
    }

    iteration++;
  }

  return x;
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

  std::cout << std::endl;
  std::cout << "Невязка до оптимизации: " << func(deltas) << std::endl;
  std::cout << std::endl;

  // Вызов градиентного спуска
  deltas = conjugateGradientDescent(func, deltas, 1e-7, 100);

  std::cout << "Невязка после оптимизации: " << func(deltas) << std::endl;

  auto pas = std::find(indsToVariate.cbegin(), indsToVariate.cend(), nodeId);

  assert(pas != indsToVariate.cend());

  return deltas(pas - indsToVariate.cbegin());
}

std::vector<int> flags(20);

int main(int numArgs, char **args) {

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
  auto solutionData = mesh.solutionData;
  const size_t numNodes = crd.size();

  checkData(mesh);

  int i = 0;

  // Цикл по каждому узлу и минимизация ошибки
  for (size_t nodeID = 0; nodeID < numNodes; ++nodeID) {
    if (i == 10) {
      break;
    }
    std::cout << std::endl;
    std::cout << "До оптимизации: " << solutionData.SigmaXX[nodeID]
              << std::endl;
    auto delta =
        static_cast<float>(tryToImproove(mesh, nodeID, ValueType::SXX));
    std::cout << std::endl;
    solutionData.SigmaXX[nodeID] += delta;
    std::cout << "Напряжение после оптимизации: "
              << solutionData.SigmaXX[nodeID] << std::endl;
    std::cout << std::endl;
    i++;
  }

  std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << std::endl;

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
    std::ifstream in(prefix + "solution0002.sba", std::ifstream::binary);

    if (!in.good())
      std::runtime_error("Ошибка открытия файла solution0002.sba");

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