#ifndef MATRIX_H
#define MATRIX_H

class LpSolution;

#include <vector>
#include <string>
#include <iostream>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

const bool verbosity = false;

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

class Matrix{
	private:
		std::vector<int> data; //we store only the indexes of the elem object

		int indexOfElement(int, int)const;
		static int indexOfInverse(int);
		Eigen::MatrixXd toEigenMatrix()const;
	protected:
		bool testPrimalEigenvectorIsParetoOptimal()const;
		bool testCosineParetoOptimal()const;
		bool testAvgSpanTreeParetoOptimal()const;
	public:
		static std::vector<double> elem;

		static void L1(std::vector<double>&);

		Matrix();
		Matrix(const std::vector<int>&);
		Matrix(int, int, int, int, int, int);

		static bool setElem(std::vector<double>);

		bool operator==(const Matrix &)const;
		bool operator<(const Matrix &)const;
		double get(const int &, const int &)const;
		Matrix permutateBy(int p[])const;
		long int getIndexOfMatrix()const;
		std::string toString()const;
		std::string toIndexString()const;
		double largestEigenvalue()const;
		double getConsistencyRatio()const;

		std::vector<double> getPrimalEigenvector()const;
		std::vector<double> getPrimalNormEigenvector()const;

		std::vector<double> getMeanOfSpans()const;

		bool testParetoOptimality(filterType) const;

		bool testVectorParetoOptimal(const std::vector<double> &) const;

		LpSolution LPVectorParetoOptimal(const std::vector<double> &) const;

		friend std::ostream& operator<<(std::ostream&, const Matrix &);
};
#endif
