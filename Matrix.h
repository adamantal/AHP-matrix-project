#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>
#include <iostream>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues" 

class Matrix{
	private:
		std::vector<int> data; //we store only the 
			//indexes of the elem object
		
		int indexOfElement(int, int)const;
		static int indexOfInverse(int);
	public:
		static std::vector<double> elem;
		
		Matrix();
		Matrix(const std::vector<int>&);
		Matrix(int, int, int, int, int, int);
		
		static bool setElem(std::vector<double>);
		
		const bool operator==(const Matrix &)const;
		const bool operator<(const Matrix &)const;
		const double get(const int &, const int &)const;
		Matrix permutateBy(int p[])const;
		long int getIndexOfMatrix()const;
		std::string toString()const;
		std::string toIndexString()const;
		double getConsistencyRatio()const;
		
		friend std::ostream& operator<<(std::ostream&, const Matrix &); 
};
#endif
