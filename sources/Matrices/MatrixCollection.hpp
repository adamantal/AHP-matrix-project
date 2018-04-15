#ifndef MATRIXCOLLECTION_HPP
#define MATRIXCOLLECTION_HPP

#include <vector>
#include <string>
#include <iostream>
#include <memory>

#include "Matrix.hpp"

//######################### Forward declarations ###############################

template<size_t N>
class MatrixCollection;

//################################# Typedefs ###################################

template <size_t N>
using MatrixCollPtr = std::shared_ptr<MatrixCollection<N>>;

//########################## class MatrixCollection ############################

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
		Matrix<N>& at(size_t);

		//utility:
		bool isIncluded(const Matrix<N>&, unsigned long long int&)const;
		void regularize ();
		void sort ();

		//core function:
		MatrixCollPtr<N> applyFilter(filterType);

		//file interaction:
		bool saveToFile(std::string filename);
		static MatrixCollPtr<N> readFromFile(std::string);

		void generateCsv(std::string, filterType);
		void printCSVWithAllData();
		//iterator:
		typename std::vector< Matrix<N> >::iterator begin();
		typename std::vector< Matrix<N> >::iterator end();
};

#include "MatrixCollection.tpp"

#endif //MATRIXCOLLECTION_HPP
