#include "Matrix.h"

std::vector<double> Matrix::elem = {};
const double ConsistencyIndex[] = {0, 1, 2, 4.049, 6.652, 9.435, 12.245, 15.045};

int Matrix::indexOfElement(int k, int l) const{
	if (k == l) {
		return 0;
	}
	if (k == 0) {
		if (l == 1) {
			return data[0];
		}
		if (l == 2) {
			return data[1];
		}
		if (l == 3) {
			return data[2];
		}
	}
	if (k == 1) {
		if (l == 0) {
			return indexOfInverse(data[0]);
		}
		if (l == 2) {
			return data[3];
		}
		if (l == 3) {
			return data[4];
		}
	}
	if (k == 2) {
		if (l == 0) {
			return indexOfInverse(data[1]);
		}
		if (l == 1) {
			return indexOfInverse(data[3]);
		}
		if (l == 3) {
			return data[5];
		}
	}
	if (k == 3) {
		if (l == 0) {
			return indexOfInverse(data[2]);
		}
		if (l == 1) {
			return indexOfInverse(data[4]);
		}
		if (l == 2) {
			return indexOfInverse(data[5]);
		}
	}
	return -1;
}
int Matrix::indexOfInverse(int i){
	if (i == 0){ 
		return 0;
	}
	if (i > 8) {
		return i - 8;
	} else {
		return i + 8;
	}
	return 0;
}

Matrix::Matrix() {
	data = std::vector<int>();
	for (int i = 0; i < 6; i++) {
		data.push_back(0);
	}
}
Matrix::Matrix(const std::vector<int> &v) {
	if (v.size() == 6) {
		data = v;
	} else {
		Matrix();
		std::cout << "Invalid amount of elements." << std::endl;
	}
}
Matrix::Matrix(int a1, int a2, int a3, int a4, int a5, int a6) {
	data.push_back(a1);
	data.push_back(a2);
	data.push_back(a3);
	data.push_back(a4);
	data.push_back(a5);
	data.push_back(a6);
}
bool Matrix::setElem(std::vector<double> arg){
	Matrix::elem = arg;
	return true;
}
const bool Matrix::operator==(const Matrix &rhs) const{
	return data == rhs.data;
}
const bool Matrix::operator<(const Matrix &rhs) const{
	for (int i = 0; i < data.size(); i++) {
		if (data[i] < rhs.data[i]) return true;
	}
	return false;
}
const double Matrix::get(const int &i, const int &j)const{
	 return Matrix::elem[indexOfElement(i,j)];
}
Matrix Matrix::permutateBy(int p[])const {
	return Matrix(indexOfElement(p[0],p[1]),
		indexOfElement(p[0],p[2]),indexOfElement(p[0],p[3]),
		indexOfElement(p[1],p[2]),indexOfElement(p[1],p[3]),
		indexOfElement(p[2],p[3]));
}
long int Matrix::getIndexOfMatrix() const {
	return 17*(17*(17*(17*(17*data[0]+data[1])+data[2])+data[3])+data[4]) + data[5];
}
std::string Matrix::toString() const{
	std::stringstream ss;
	ss.precision(5);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			ss << Matrix::elem[indexOfElement(i,j)] << "\t";
		}
		if (i != 3) {
			ss << "\n";
		}
	}
	return ss.str();
}
std::string Matrix::toIndexString()const{
	std::stringstream ss;
	//ss.precision(5);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			ss << indexOfElement(i,j) << "\t";
		}
		if (i != 3) {
			ss << "\n";
		}
	}
	return ss.str();
}
std::ostream& operator<<(std::ostream& os, const Matrix& m){
	os << m.toString();
	return os;
}

double Matrix::getConsistencyRatio()const{
	Eigen::MatrixXd m(4,4);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			m(i,j) = elem[indexOfElement(i,j)];
		}
	}
	//std::cout << m << std::endl;
	
	Eigen::VectorXcd eigenvals = m.eigenvalues();
	double largestEigVal = eigenvals[0].real();
	for (int i = 1; i < eigenvals.rows(); i++) {
		if (largestEigVal < eigenvals[i].real()) {
			largestEigVal = eigenvals[i].real();
		}
	}
	
	return (largestEigVal-m.rows())/(ConsistencyIndex[m.rows()]-m.rows());
}
