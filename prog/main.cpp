#include <fstream>
#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <set>

using namespace std;

int main()
{
	setlocale(LC_ALL, "rus");
	// crd
	struct Node {
		float x, y, z;
	};
	std::vector<Node> crd;
	{
		std::ifstream in("C:/Users/Дима/Desktop/b/toCopy/brick/crd.sba", std::ifstream::binary);
		if (!in.good())
			throw std::runtime_error("Crd");
		int numNodes, kort;
		in.read((char*)&numNodes, sizeof(int));
		in.read((char*)&kort, sizeof(int));
		numNodes /= kort;
		crd.resize(numNodes);
		std::vector<float> tmp;
		tmp.resize(numNodes * kort);
		in.read((char*)tmp.data(), sizeof(float) * tmp.size());

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
		std::ifstream in("C:/Users/Дима/Desktop/b/toCopy/brick/ind.sba", std::ifstream::binary);
		if (!in.good())
			throw std::runtime_error("Ind");
		int numTypes;
		in.read((char*)&numTypes, sizeof(int));
		std::vector<std::pair<int, int>> types;
		types.resize(numTypes);
		for (auto& t : types) {
			in.read((char*)&(t.first), sizeof(int));
			std::cout << t.first << std::endl;
			in.read((char*)&(t.second), sizeof(int));
			std::cout << t.second << std::endl;
		}
		// read

		for (auto t : types) {
			int numIndexes;
			in.read((char*)&numIndexes, sizeof(int)); // Считываем длину массива индексов
			std::cout << numIndexes << std::endl;
			for (int indexElement = 0; indexElement < t.second; ++indexElement) {
				indexes ind;
				ind.resize(numIndexes / t.second);
				in.read((char*)ind.data(), sizeof(int) * ind.size());
				elems.push_back(std::make_pair(indexElement, ind)); // Заполняем вектор elems парами номер КЭ - список узлов, которые в него включены
			}
		}
	}

	// Создаем мапу для хранения индексации узлов
	map<int, set<int>> elemsOfNodeID;

	// Заполняем мапу
	for (const auto& el : elems) {
		for (const auto& node : el.second) {
			elemsOfNodeID[node].insert(el.first);
		}
	}

	// Выводим результат
	for (const auto& it : elemsOfNodeID) {
		cout << it.first << " - "; // it.first - ключ (номер узла)
		for (const auto& el : it.second) { // it.second - набор индексов КЭ
			cout << el << " ";
		}
		cout << endl;
	}

	system("pause");

	return 0;
}