#ifndef MATRIXCOLLECTION_H
#define MATRIXCOLLECTION_H

#include <vector>
#include <string>
#include <iostream>
#include <set>
#include "Matrix.h"

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

		//core function:
		MatrixCollection applyFilter(filterType);

		//file interaction:
		bool saveToFile(std::string filename);
		static MatrixCollection readFromFile(std::string);

		void generateCsv(std::string);
		//iterator:
		typename std::vector< Matrix<N> >::iterator begin();
		typename std::vector< Matrix<N> >::iterator end();
};
#endif
