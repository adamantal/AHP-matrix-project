namespace matrixInit {
	const long int limit = 24137569;

	MatrixCollection<4> getAllMatrices() {
		//std::cout << "Loading all the matrices... " << std::endl;
		MatrixCollection<4> vec;
		long int i = 0;
		for (unsigned short i1 = 0; i1 < Matrix<4>::elem.size(); i1++) {
			//std::cout << "Loading datas: " << i1 << "/" << Matrix::elem.size() << std::endl;
			for (unsigned short i2 = 0; i2 < Matrix<4>::elem.size(); i2++) {
				for (unsigned short i3 = 0; i3 < Matrix<4>::elem.size(); i3++) {
					for (unsigned short i4 = 0; i4 < Matrix<4>::elem.size(); i4++) {
						for (unsigned short i5 = 0; i5 < Matrix<4>::elem.size(); i5++) {
							for (unsigned short i6 = 0; i6 < Matrix<4>::elem.size(); i6++) {
								std::vector<Ush> k = {i1,i2,i3,i4,i5,i6};
								Matrix<4> tmp(k);
								vec.add(tmp);
								i++;
							}
						}
					}
				}
			}
		}
		return vec;
	}

	template<size_t N>
	void takeOutPermutationOfIt(MatrixCollection<N> &vec, size_t &i, std::vector<bool> &check) {
		Ush perm[] = {0,1,2,3};
		while ( std::next_permutation(perm,perm+4) ) {
			Matrix<N> tmp = vec[i].permutateBy(perm);
			size_t matrixi = tmp.getIndexOfMatrix();
			check[matrixi] = true;
		}
	}

	template<size_t N>
	MatrixCollection<N> takeOutAllPermutations(MatrixCollection<N> &vecold, std::vector<bool> &check) {
		size_t k = 0;
		MatrixCollection<N> vecnew;
		for (auto it = vecold.begin(); it != vecold.end(); it++) {
			/*if (k % 241376 == 0) {
				std::cout << k/241376.0 << " percent has so far" << std::endl;
			}*/
			size_t tmpi = it - vecold.begin();
			if (!check[tmpi]) {
				vecnew.add(vecold[tmpi]);
				takeOutPermutationOfIt(vecold, tmpi, check);
			}
			k++;
		}
		return vecnew;
	}

	template<size_t N>
	bool generateAllToFile() {
		std::cout << "Procedure started.\n";
		std::vector<bool> checkedMatrices(limit);
		for (std::vector<bool>::iterator it = checkedMatrices.begin();
			it != checkedMatrices.end(); it++) {
			*it = false;
		}

		std::cout << "Getting all matrices.\n";
		MatrixCollection<N> vec = getAllMatrices();

		std::cout << "The size is " << ((vec.size() == limit)? "correct." : "not correct.") << std::endl;
		std::cout << "Taking out duplicated elements.\n";
		MatrixCollection<N> vecfiltered = takeOutAllPermutations(vec, checkedMatrices);
		std::cout << "Saving to file " << vecfiltered.size() << " elements.\n";

		vecfiltered.saveToFile("../res/allMatrices.mt");
		return true;
	}
}
