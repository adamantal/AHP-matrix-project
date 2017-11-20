#ifndef MATRIX_H
#define MATRIX_H

template<size_t N>
class LpSolution;

#include <vector>
#include <string>
#include <iostream>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

const bool verbosity = false;

typedef unsigned short Ush;

struct spair{
	int a;
	int b;

	spair(int x, int y):a(x),b(y){}
};

enum filterType {
	Inconsistency = 0,
	Consistency = 1,
	EigenVectorMethod = 2,
	AverageSpanTreeMethod = 3,
	CosineMethod = 4
};

template<size_t N>
class Matrix{
	protected:
		//we store only the indexes of the elem object:
		std::vector<Ush> data;

		Ush indexOfElement(Ush, Ush)const;
		static Ush indexOfInverse(Ush);
		Eigen::MatrixXd toEigenMatrix()const;

	protected:
		bool testPrimalEigenvectorIsParetoOptimal()const;
		bool testCosineParetoOptimal()const;
		bool testAvgSpanTreeParetoOptimal()const;

	public:
		//the same for all the matrices:
		static const std::vector<double> elem;
		//static functions:
		static void L1(std::vector<double>&);

		//constructors:
		Matrix();
		Matrix(const std::vector<Ush>&);

		//operators:
		bool operator==(const Matrix &)const;
		bool operator!=(const Matrix & rhs)const{return !(*this == rhs);};
		bool operator<(const Matrix &)const;
		bool operator>(const Matrix & rhs)const{return !(*this < rhs && *this != rhs);};
		bool operator<=(const Matrix & rhs)const{return !(*this > rhs);};
		void operator++(int);

		//getters:
		std::vector<Ush> getData()const;
		double get(Ush, Ush)const;
		unsigned long long int getIndexOfMatrix()const;
		static Matrix getMatrixOfIndex(unsigned long long int);
		bool isMinimalPermutated()const;

		//for generating 5x5 matrices:
		Matrix<N - 1> cutBottom()const;

		//IO and its necesssary conversions:
		std::string toString(bool index = false) const;
		template<size_t M> friend std::ostream& operator<<(std::ostream&, const Matrix<M> &);

		//functional part:
		Matrix<N> permutateBy(Ush p[])const;
		double largestEigenvalue()const;
		double getConsistencyRatio()const;

		std::vector<double> getPrimalEigenvector()const;
		std::vector<double> getPrimalNormEigenvector()const;
		std::vector<double> getMeanOfSpans()const;
		std::vector<double> getCosineVector()const;

		bool testParetoOptimality(filterType) const;
		bool testVectorParetoOptimal(const std::vector<double> &) const;
		LpSolution<N> LPVectorParetoOptimal(const std::vector<double> &) const;
};
#endif
