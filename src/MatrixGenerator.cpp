#include "MatrixCollection.h"

namespace matrixInit {
	//const long int limit = 24137569;

	template<class T>
	T expn(T x, unsigned int y) {
		T r = 1;
		for (unsigned int i = 0; i < y; i++) r *= x;
		return r;
	}

	bool getNextVector(std::vector<Ush> &input) {
		for (auto rit = input.rbegin(); rit != input.rend(); rit++) {
			if (std::distance(rit, input.rend()) == 1) std::cout << "Last element touched\n";

			if (*rit == Matrix<1>::elem.size() - 1) {
				*rit = 0;
			}
			else {
				*rit = *rit + 1;
				return true;
			}
		}
		return false;
	}

	template<size_t N>
	MatrixCollection<N> getAllUnfilteredMatrices() {
		if (N <= 2) throw "Unfiltered matrices can't be generated for N <= 2 dimensions.\n";

		MatrixCollection<N> mc;
		std::vector<Ush> v;
		for (size_t i = 0; i < N * (N - 1) / 2; i++) v.push_back(0);

		unsigned long long int i = 0;
		do {
			i++;
			if (i % 1000000 == 0) std::cout << i / 1000000 << " million.\n";
			Matrix<N> tmp(v);
			mc.add(tmp);
		} while (getNextVector(v));

		return mc;
	}

	template<size_t N>
	void takeOutPermutationOfIt(const Matrix<N> &m, std::vector<bool> &check) {
		Ush perm[N];
		for (Ush j = 0; j < N; j++) perm[j] = j;

		while ( std::next_permutation(perm, perm + N) ) {
			Matrix<N> tmp = m.permutateBy(perm);
			size_t matrixi = tmp.getIndexOfMatrix();
			if (expn(Matrix<0>::elem.size(), N * (N - 1) / 2) + 1 < matrixi) throw "Out of bounds index was about to be used!\n";
			check[matrixi] = true;
		}
	}

	template<size_t N>
	MatrixCollection<N> takeOutAllPermutations(MatrixCollection<N> &vecold) {
		std::vector<bool> checkedMatrices(expn(Matrix<0>::elem.size(), N * (N - 1) / 2) + 1);
		for (std::vector<bool>::iterator it = checkedMatrices.begin();
			it != checkedMatrices.end(); it++) {
			*it = false;
		}

		MatrixCollection<N> vecnew;
		for (auto it = vecold.begin(); it != vecold.end(); it++) {
			size_t tmpi = it - vecold.begin();
			if (!checkedMatrices[tmpi]) {
				vecnew.add(vecold[tmpi]);
				takeOutPermutationOfIt(vecold[tmpi], checkedMatrices);
			}
			//k++;
		}
		return vecnew;
	}

	template<size_t N>
	bool generateAllToFile();

	template<>
	bool generateAllToFile<5>(){
		size_t sizeN = expn(Matrix<0>::elem.size(), 5 * (5 - 1) / 2);

		std::cout << "Generating vector...\n";

		std::cout << sizeN << std::endl;
		std::vector<bool> checkedMatrices(sizeN);
		for (std::vector<bool>::iterator it = checkedMatrices.begin();
			it != checkedMatrices.end(); it++) {
			std::cout << std::distance(checkedMatrices.begin(), it) << "\n";
			*it = false;
		}

		std::cout << "Iterating through indexes...\n";
		MatrixCollection<5> matrices;
		for (size_t i = 0; i < sizeN; i++) {
			if (i % (sizeN / 100) == 0) std::cout << i / (sizeN / 100) << " percent!\n";
			if (!checkedMatrices[i]) {
				Matrix<5> m = Matrix<5>::getMatrixOfIndex(i);
				matrices.add(m);

				takeOutPermutationOfIt(m, checkedMatrices);
			}
		}

		std::cout << "Saving to file " << matrices.size() << " elements.\n";

		matrices.saveToFile("../res/all" + std::to_string(5) + "matrices.mt");
		return true;
	}

	template<>
	bool generateAllToFile<4>() {
		std::cout << "Procedure started\n";

		std::cout << "Getting all matrices...\n";
		MatrixCollection<4> all = getAllUnfilteredMatrices<4>();

		std::cout << "The size is " << ((all.size() == expn(Matrix<0>::elem.size(), 4 * (4 - 1) / 2)) ? "correct." : "not correct.") << std::endl;
		if (all.size() != expn(Matrix<0>::elem.size(), 4 * (4 - 1) / 2)) throw "The sizes are different - error in the method please check.\n";

		std::cout << "Taking out duplicated elements...\n";
		MatrixCollection<4> vecfiltered = takeOutAllPermutations(all);
		std::cout << "Saving to file " << vecfiltered.size() << " elements.\n";

		vecfiltered.saveToFile("../res/all4matrices.mt");
		return true;
	}
}
