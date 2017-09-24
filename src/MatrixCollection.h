#ifndef MATRIXCOLLECTION_H
#define MATRIXCOLLECTION_H

class Matrix;

#include <vector>
#include <string>
#include <iostream>
#include <set>

typedef std::vector<Matrix>::iterator iterator;

class MatrixCollection {
	private:
		std::vector<Matrix> data;
		size_t tmpindex;
	public:
		MatrixCollection();
		//MatrixCollection(vector<Matrix>);

		//inherited function from vector:
		size_t size();
		void add(Matrix &);
		Matrix& operator[](size_t);

		//core functions:
		MatrixCollection applyInconsistencyFilter();
		MatrixCollection applyConsistencyFilter();
		MatrixCollection applyEigenvalueMethodFilter();

		MatrixCollection applyAvgSpanTreeFilter();

		//file interaction:
		bool saveToFile(std::string filename);
		static MatrixCollection readFromFile(std::string);

		void generateCsv(std::string);
		//iterator:
		iterator begin();
		iterator end();
		//void erase(iterator);

		//bool checkAssumption();
};
#endif
