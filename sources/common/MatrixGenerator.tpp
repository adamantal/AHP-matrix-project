#include <vector>
#include <bitset>

#include "Matrix.hpp"
#include "MatrixCollection.hpp"
#include "Paths.hpp"

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

	/*void takeOutPermutationOfItSpec5(const Matrix5Spec &m, std::bitset<1> &s, const MatrixCollection<4> mc, const size_t& outerIndex) {
		//template specify:
		Ush N = 5;
		static unsigned long long int COUNT = 0;

		Ush perm[N];
		for (Ush j = 0; j < N; j++) perm[j] = j;

		Matrix<5> M5 = m.getMatrix();

		while ( std::next_permutation(perm, perm + N) ) { //5!
			Matrix<5> tmp = M5.permutateBy(perm);
			unsigned long long int i = 0;
			if (Matrix5Spec::getIndexOutOfMatrix(tmp, i, mc)) {
				if (outerIndex > i) {
					std::cout << "outerIndex: " << outerIndex << std::endl;
					std::cout << "i: " << i << std::endl;
					std::cout << "Original matrix:\n" << M5 << std::endl << std::endl;
					std::cout << "Permutated matrix:\n" << tmp << std::endl << std::endl;
					throw "Error - scrolling through previously watched elements.\n";
				}
				if (!s[i]) COUNT++;
				s[i] = true;
			}
		}
		if (COUNT % 1000 < 100) std::cout << COUNT << std::endl;
	}*/

	/*template<>
	bool generateAllToFile<5>() {
		std::cout << "Procedure started.\n";

		std::cout << "Reading 4x4 files...\n";
		MatrixCollPtr<4> mc = MatrixCollection<4>::readFromFile(PATH_4_ALL);

		std::cout << "Collecting memory for vector...\n";
		//std::vector<bool>

		//std::bitset<84113665016> s; //UNDO THIS
		std::bitset<1> s;

		//std::vector<bool> s(84113665016); // 17^4 * 1007096

		std::cout << "Vector initialized...\n";

		//checking validity:
		{
			std::cout << "\tChecking validity... \n";
			std::vector<Ush> v = {0, 0, 0, 0, 0, 0};
			Matrix5Spec m5(v);
			unsigned long long int i = 99;
			Matrix5Spec::getIndexOutOfMatrix(m5.getMatrix(), i, mc);
			if (i != 0) throw "Error: test failed - please check getIndexOutOfMatrix function.\n";
			std::cout << "\tTEST OK\n";
		}

		std::cout << "Filtering...\n";
		for (unsigned long long int tmpi = 0; tmpi < 84113665016; tmpi++) {
			//tmpi = it - s.begin();
			if (!s[tmpi]) {
				takeOutPermutationOfItSpec5(Matrix5Spec::getMatrixOutOfIndex(tmpi, mc), s, mc, tmpi);
			}
		}
		return true;
	}*/

	template<>
	bool generateAllToFile<4>() {
		std::cout << "Procedure started\n";

		std::cout << "Getting all matrices...\n";
		MatrixCollection<4> all = getAllUnfilteredMatrices<4>();

		std::cout << "The size is " << ((all.size() == expn(Matrix<0>::elem.size(), 4 * (4 - 1) / 2)) ? "correct." : "not correct.") << std::endl;
		if (all.size() != expn(Matrix<0>::elem.size(), 4 * (4 - 1) / 2)) throw "The sizes are different - error in the method please check.\n";

		std::cout << "Taking out duplicated elements...\n";
		MatrixCollection<4> vecfiltered = takeOutAllPermutations(all);

		std::cout << "Regularizing...\n";
		vecfiltered.regularize ();

		std::cout << "Sorting with custom sort...\n";
		vecfiltered.sort();

		std::cout << "Saving to file " << vecfiltered.size() << " elements.\n";

		vecfiltered.saveToFile(PATH_4_ALL);
		return true;
	}
}
