#include <iostream>
#include <vector>
#include <set>

#include <lemon/lp.h>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

// ezek az elemek a spec PCMC-kre vonatkoznak
const double ConsistencyIndex[] = {0, 1, 2, 4.049, 6.652, 9.435, 12.245, 15.045};

template<typename T>
double getLargestEigenvalue(T m)
{
	Eigen::VectorXcd eigenvals = m.eigenvalues();
	double largestEigenvalue = eigenvals[0].real();
	for (int i = 1; i < eigenvals.rows(); i++) {
		if (largestEigenvalue < eigenvals[i].real()) {
			largestEigenvalue = eigenvals[i].real();
		}
	}
    return largestEigenvalue;
}

template<typename T>
std::vector<double> getLargestEigenvector(T m)
{
	Eigen::EigenSolver<T> eigenSolver(m);

	double largestEigenvalue = eigenSolver.eigenvalues()[0].real();

	auto largestEigenvector = eigenSolver.eigenvectors().col(0);
	for (int i = 1; i < m.rows(); i++) {
		if (largestEigenvalue < eigenSolver.eigenvalues()[i].real()) {
			largestEigenvalue = eigenSolver.eigenvalues()[i].real();
			largestEigenvector = eigenSolver.eigenvectors().col(i);
		}
	}
	std::vector<double> v;
	for (int i = 0; i < 4; i++) {
		v.push_back(largestEigenvector(i).real());
	}
    return v;
}

template<typename T>
double getConsistencyRatio(T m)
{
    return (getLargestEigenvalue(m)-m.rows())/(ConsistencyIndex[m.rows()]-m.rows());
}

template<typename T>
void vectorOutput(const std::vector<T> &v) {
	std::cout << "(";
	for (int i = 0; i < v.size() - 1; i++) {
		std::cout << v[i] << " ";
	}
	std::cout << v[v.size()-1] <<")";
}

int main()
{

	Eigen::Matrix4d M;
	M << 1.0, 7.0, 6.0, 5.0,
		1.0/7, 1.0, 1.0/2, 1.0,
		1.0/6, 2.0, 1.0, 1.0/2,
		1.0/5, 1.0, 2.0, 1.0;
	std::cout.precision(4);
	std::cout << "Lets consider the following matrix: \n";
	std::cout << M << std::endl;

	std::cout << "The consistency ratio of the matrix is: ";
	std::cout << getConsistencyRatio(M) << std::endl;
	if (getConsistencyRatio(M) < 0.1) {
		std::cout << "The matrix meets Salty's consistency criteria.\n";
		std::vector<double> e = getLargestEigenvector(M);

		//It normalizes the given eigenvector
		std::vector<double> w(4);
		for (int i = 0; i < 4; i++) {
			w[i] = std::abs(e[i]/e[0]);
		}

		std::set<int> I0a,I0b;
		std::set<int> I1a,I1b;
		std::vector<double> v(4);
		std::vector<std::vector<double>> B(4, std::vector<double>(4));
		for (int i = 0; i < 4; i++) {
			v[i] = log(w[i]);
			for (int j = 0; j < 4; j++) {
				B[i][j] = log(M(i,j));
				if (w[i]/w[j] - M(i,j) > 10^(-6)) {
					I1a.insert(i);
					I1b.insert(j);
				}
			}
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				if (abs(w[i]/w[j] - M(i,j)) < 10^(-8)) {
					I0a.insert(i);
					I0b.insert(j);
				} // valami miatt nem jó
			}
		}
		std::cout << I0a.size() << std::endl;
		std::cout << I1a.size() << std::endl;
	} else {
		std::cout << "The matrix fails Saaty's consistency criteria.\n"
		<< "Progress terminated\n";
	}

	return 0;
}
