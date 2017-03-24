#ifndef MATRIXCOLLECTION_H
#define MATRIXCOLLECTION_H

#include <vector>
#include <string>
#include <iostream>
#include <set>
#include "Matrix.h"

typedef std::vector<Matrix>::iterator iterator;

class MatrixCollection {
	private:
		std::vector<Matrix> data;
		size_t tmpindex;
	public:
		MatrixCollection();
		//MatrixCollection(vector<Matrix>);
		size_t size();
		void add(Matrix &);
		//void erase(iterator);
		Matrix& operator[](size_t);
		bool saveToFile(std::string filename);
		static MatrixCollection readFromFile(std::string);
		
		iterator begin();
		iterator end();
};
#endif
