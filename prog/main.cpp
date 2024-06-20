//
// КОДИРОВКА ТЕКСТА СТРОГО UTF-8
//

#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>

using namespace std;

int main() {
  // setlocale(LC_ALL, "rus");
  // crd
  struct Node {
    float x, y, z;
  };
  std::vector<Node> crd;
  {
    // Не должно быть жестко прописанных путей в коде
    std::ifstream in("crd.sba", std::ifstream::binary);
    if (!in.good())
      throw std::runtime_error("Crd");
    int numNodes, kort;
    in.read((char *)&numNodes, sizeof(int));
    in.read((char *)&kort, sizeof(int));
    numNodes /= kort;
    crd.resize(numNodes);
    std::vector<float> tmp;
    tmp.resize(numNodes * kort);
    in.read((char *)tmp.data(), sizeof(float) * tmp.size());

    for (size_t ct = 0; ct < numNodes; ++ct) {
      crd[ct].x = tmp[0 * numNodes + ct];
      crd[ct].y = tmp[1 * numNodes + ct];
      crd[ct].z = tmp[2 * numNodes + ct];
    }
  }
  const size_t numNodes = crd.size();

  // mesh crd info

  using indexes = std::vector<int>;

  std::vector<std::pair<int, indexes>> elems;
  // ind
  {
    // Не должно быть жестко прописанных путей в коде
    std::ifstream in("ind.sba", std::ifstream::binary);
    if (!in.good())
      throw std::runtime_error("Ind");
    int numTypes;
    in.read((char *)&numTypes, sizeof(int));
    std::vector<std::pair<int, int>> types;
    types.resize(numTypes);
    for (auto &t : types) {
      in.read((char *)&(t.first), sizeof(int));
      std::cout << t.first << std::endl;
      in.read((char *)&(t.second), sizeof(int));
      std::cout << t.second << std::endl;
    }
    // read

    for (auto t : types) {
      int numIndexes;
      in.read((char *)&numIndexes,
              sizeof(int)); // Ñ÷èòûâàåì äëèíó ìàññèâà èíäåêñîâ
      std::cout << numIndexes << std::endl;
      for (int indexElement = 0; indexElement < t.second; ++indexElement) {
        indexes ind;
        ind.resize(numIndexes / t.second);
        in.read((char *)ind.data(), sizeof(int) * ind.size());
        // Читерство. У нас элементы лежат в векторе - контейнере с произвольным
        // доступом. Нет никакого смысла класть туда еще и номер элемента. Тут
        // должен был быть признак типа для отличения элементов по типу.
        elems.push_back(std::make_pair(
            indexElement, ind)); // Çàïîëíÿåì âåêòîð elems ïàðàìè íîìåð ÊÝ -
                                 // ñïèñîê óçëîâ, êîòîðûå â íåãî âêëþ÷åíû
      }
    }
  }

  // Ñîçäàåì ìàïó äëÿ õðàíåíèÿ èíäåêñàöèè óçëîâ
  map<int, set<int>> elemsOfNodeID;

  // Çàïîëíÿåì ìàïó
  for (const auto &el : elems) {
    for (const auto &node : el.second) {
      // См. выше - в первом поле должен лежать признак элемента (его тим), а не
      // номер элемента
      elemsOfNodeID[node].insert(el.first);
    }
  }

  // Âûâîäèì ðåçóëüòàò
  for (const auto &it : elemsOfNodeID) {
    cout << it.first << " - "; // it.first - êëþ÷ (íîìåð óçëà)
    for (const auto &el : it.second) { // it.second - íàáîð èíäåêñîâ ÊÝ
      cout << el << " ";
    }
    cout << endl;
  }

  // Это зачем?
  // system("pause");

  return 0;
}