#ifndef MATRIXCOLLECTION_H
#define MATRIXCOLLECTION_H

#include <vector>
#include <string>
#include <iostream>
#include <set>
#include "Matrix.h"

const std::string PATH_4_ALL = "../res/all4x4matrices.mt";
const std::string PATH_4_CONSISTENT = "../res/consistent4x4matrices.mt";

const std::string PATH_4_EIGEN = "../res/eigen4x4matrices.mt";
const std::string PATH_4_EIGEN_CSV = "../res/eigen4x4matrices.csv";
const std::string PATH_4_SPANTREE = "../res/spantree4x4matrices.mt";
const std::string PATH_4_SPANTREE_CSV = "../res/spantree4x4matrices.csv";
const std::string PATH_4_COSINE = "../res/cosine4x4matrices.mt";
const std::string PATH_4_COSINE_CSV = "../res/cosine4x4matrices.csv";

const std::string PATH_5_ALL = "../res/all5x5matrices.mt";

template<size_t N>
class MatrixCollection {
	private:
		std::vector< Matrix<N> > data;
		size_t tmpindex;
	public:
		MatrixCollection();
		//MatrixCollection(vector<Matrix>);

		//inherited function from vector:
		size_t size();
		void add(Matrix<N> &);
		Matrix<N>& operator[](size_t);

		//utility:
		bool isIncluded(const Matrix<N>&, unsigned long long int&)const;
		void regularize ();
		void sort ();

		//core function:
		MatrixCollection applyFilter(filterType);

		//file interaction:
		bool saveToFile(std::string filename);
		static MatrixCollection readFromFile(std::string);

		void generateCsv(std::string, filterType);
		void printCSVWithAllData();
		//iterator:
		typename std::vector< Matrix<N> >::iterator begin();
		typename std::vector< Matrix<N> >::iterator end();
};
#endif
