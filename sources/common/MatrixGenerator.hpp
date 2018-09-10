#ifndef MATRIXGENERATOR_HPP
#define MATRIXGENERATOR_HPP

#include <vector>

template<size_t N>
class Matrix;

template<size_t N>
class MatrixCollection;

namespace matrixInit {
	template<class T>
	T expn(T x, unsigned int y);

	bool getNextVector(std::vector<Ush> &input);

	template<size_t N>
	MatrixCollection<N> getAllUnfilteredMatrices();

	template<size_t N>
	void takeOutPermutationOfIt(const Matrix<N> &m, std::vector<bool> &check);

	template<size_t N>
	MatrixCollection<N> takeOutAllPermutations(MatrixCollection<N> &vecold);

	template<size_t N>
	bool generateAllToFile();
}

#include "MatrixGenerator.tpp"

#endif //MATRIXGENERATOR_HPP
