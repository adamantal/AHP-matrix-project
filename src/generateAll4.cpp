namespace matrixInit {
	const long int limit = 24137569;

	MatrixCollection getAllMatrices() {
		//std::cout << "Loading all the matrices... " << std::endl;
		MatrixCollection vec;
		long int i = 0;
		for (int i1 = 0; i1 < Matrix::elem.size(); i1++) {
			//std::cout << "Loading datas: " << i1 << "/" << Matrix::elem.size() << std::endl;
			for (int i2 = 0; i2 < Matrix::elem.size(); i2++) {
				for (int i3 = 0; i3 < Matrix::elem.size(); i3++) {
					for (int i4 = 0; i4 < Matrix::elem.size(); i4++) {
						for (int i5 = 0; i5 < Matrix::elem.size(); i5++) {
							for (int i6 = 0; i6 < Matrix::elem.size(); i6++) {
								std::vector<int> k = {i1,i2,i3,i4,i5,i6};
								Matrix tmp(k);
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

	void takeOutPermutationOfIt(MatrixCollection &vec, size_t &i, std::vector<bool> &check) {
		int perm[] = {0,1,2,3};
		while ( std::next_permutation(perm,perm+4) ) {
			Matrix tmp = vec[i].permutateBy(perm);
			size_t matrixi = tmp.getIndexOfMatrix();
			check[matrixi] = true;
		}
	}

	MatrixCollection takeOutAllPermutations(MatrixCollection &vecold, std::vector<bool> &check) {
		size_t k = 0;
		MatrixCollection vecnew;
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
		std::cout << k << "elements\n";
		return vecnew;
	}

	void setElement(){
		std::vector<double> tmpvec = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,1.0/2,1.0/3,1.0/4,1.0/5,1.0/6,1.0/7,1.0/8,1.0/9};
		Matrix::setElem(tmpvec);
	}

	bool generateAllToFile() {
		std::vector<bool> checkedMatrices(limit);
		for (std::vector<bool>::iterator it = checkedMatrices.begin();
			it != checkedMatrices.end(); it++) {
			*it = false;
		}

		MatrixCollection vec = getAllMatrices();

		std::cout << (vec.size() == limit) << std::endl;

		MatrixCollection vecfiltered = takeOutAllPermutations(vec, checkedMatrices);

		std::cout << vecfiltered.size() << std::endl;

		vecfiltered.saveToFile("../bin/output.mt");
		return true;
	}
}
