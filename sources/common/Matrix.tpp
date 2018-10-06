#include <algorithm>

#include "lemon/lp.h"

#include "Matrix.hpp"
#include "GraphSolution.hpp"
#include "LpSolution.hpp"

template<size_t N>
const std::vector<double> Matrix<N>::elem = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,1.0/2,1.0/3,1.0/4,1.0/5,1.0/6,1.0/7,1.0/8,1.0/9};

template<size_t N>
const double Matrix<N>::ConsistencyIndex[] = {0.0, 1.0, 2.0, 4.049, 6.652, 9.435, 12.245, 15.045};

template<size_t N>
Ush Matrix<N>::indexOfElement(Ush i, Ush j) const {
	if (i == j)
		return 0;
	else if (i < j)
		return data[j - 1 + (2 * N - i - 3) * i / 2];
	else {
		return indexOfInverse(indexOfElement(j, i));
	}
}

template<size_t N>
Ush Matrix<N>::indexOfInverse(Ush i) {
	if (i == 0)
		return 0;
	if (i > 8)
		return i - 8;
	else
		return i + 8;
}

template<size_t N>
std::vector<Ush> Matrix<N>::getData()const {
	return data;
}

template<size_t N>
Matrix<N>::Matrix() {
	data = std::vector<Ush>();
	for (size_t i = 0; i < N; i++) {
		data.push_back(0);
	}
}

template<size_t N>
Matrix<N>::Matrix(const std::vector<Ush> &v) {
	if (v.size() == N * (N - 1) / 2)
		data = v;
	else
		throw "Invalid amount of elements!\n";
}

template<size_t N>
bool Matrix<N>::operator==(const Matrix &rhs) const{
	//it is OK! C++11
	return data == rhs.data;
}

template<size_t N>
bool Matrix<N>::operator<(const Matrix<N> &rhs) const{
	for (size_t i = 0; i < N * (N - 1) / 2; i++) {
		if (data[i] < rhs.data[i])
			return true;
		else if (data[i] > rhs.data[i])
			return false;
	}
	return false;
}

template<size_t N>
void Matrix<N>::operator++(int) {
	for (auto rit = data.rbegin(); rit != data.rend(); rit++) {
		if (*rit == Matrix<1>::elem.size() - 1) {
			*rit = 0;
		}
		else {
			*rit = *rit + 1;
			return;
		}
	}
}

template<size_t N>
double Matrix<N>::get(Ush i, Ush j) const {
	 return Matrix<N>::elem[indexOfElement(i,j)];
}

template<size_t N>
Matrix<N> Matrix<N>::permutateBy(Ush p[])const {
	//TODO: resolve this problem here!!!
	std::vector<Ush> v;
	for (size_t i = 0; i < N - 1; i++) {
		for (size_t j = i + 1; j < N; j++) {
			v.push_back(indexOfElement(p[i], p[j]));
		}
	}
	return Matrix<N>(v);
}

//TODO have to check and compare with INT_MAX
template<size_t N>
unsigned long long int Matrix<N>::getIndexOfMatrix() const {
	unsigned long long int x = data[0];
	for (size_t i = 1; i < N * (N - 1) / 2; i++) {
		x = Matrix<0>::elem.size() * x + data[i];
	}

	return x;
}

template<size_t N>
Matrix<N> Matrix<N>::getMatrixOfIndex(unsigned long long int x) {
	std::vector<Ush> v;
	for (size_t i = 0; i < N * (N - 1) / 2; i++) {
		Ush tmp = x % Matrix<0>::elem.size();
		x /= Matrix<0>::elem.size();
		v.push_back(tmp);
	}
	std::reverse(v.begin(), v.end());
	return Matrix(v);
}

template<size_t N>
bool Matrix<N>::isMinimalPermutated() const {
	Ush perm[N];
	for (Ush j = 0; j < N; j++) perm[j] = j;

	//quickTest procedure: IT MIGHT IMPROVE PERFORMANCE
	/*std::random_shuffle(perm, perm + N);
	Matrix<N> tmp = this -> permutateBy(perm);
	if (*this > tmp)
		return false;
	for (Ush j = 0; j < N; j++) perm[j] = j;*/

	while ( std::next_permutation(perm, perm + N) ) {
		Matrix<N> tmp = this -> permutateBy(perm);
		if (*this > tmp) {
			return false;
		}
	}
	return true;
}

template<size_t N>
Matrix<N> Matrix<N>::getItsMinimalPermutate() const {
	Matrix<N> minm (*this);

	Ush perm[N];
	for (Ush j = 0; j < N; j++) perm[j] = j;

	while ( std::next_permutation(perm, perm + N) ) {
		Matrix<N> tmp = this -> permutateBy(perm);
		if (minm > tmp) {
			minm = tmp;
		}
	}
	return minm;
}

template<size_t N>
std::string Matrix<N>::toString(bool index /*= false*/) const {
	std::stringstream ss;
	ss.precision(5);
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			if (index)
				ss << indexOfElement(i,j) << "\t";
			else
				ss << Matrix::elem[indexOfElement(i,j)] << "\t";

		}
		if (i != N - 1) {
			ss << "\n";
		}
	}
	return ss.str();
}

template<size_t N>
std::ostream& operator<<(std::ostream& os, const Matrix<N>& m) {
	os << m.toString();
	return os;
}

template<size_t N>
Eigen::MatrixXd Matrix<N>::toEigenMatrix()const{
	Eigen::MatrixXd m;
	m = Eigen::MatrixXd(N,N);
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			m(i,j) = Matrix::elem[indexOfElement(i,j)];
		}
	}
	return m;
}

template<size_t N>
double Matrix<N>::largestEigenvalue()const{
	Eigen::MatrixXd m(N, N);
	m = toEigenMatrix();

	Eigen::VectorXcd eigenvals = m.eigenvalues();
	double largestEigVal = eigenvals[0].real();
	for (int i = 1; i < eigenvals.rows(); i++) {
		if (largestEigVal < eigenvals[i].real()) {
			largestEigVal = eigenvals[i].real();
		}
	}
	return largestEigVal;
}

template<size_t N>
double Matrix<N>::getConsistencyRatio()const{
	Eigen::MatrixXd m(N, N);
	m = toEigenMatrix();

	return (largestEigenvalue()-m.rows())/(ConsistencyIndex[m.rows()]-m.rows());
}

template<size_t N>
std::vector<double> Matrix<N>::getPrimalEigenvector()const{
	Eigen::MatrixXd m(N, N);
	m = toEigenMatrix();

	Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(m);

	double largestEigenvalue = eigenSolver.eigenvalues()[0].real();
	int which = 0;
	short multiplicity = 1;

	for (int i = 1; i < eigenSolver.eigenvalues().size(); i++) {
		if (largestEigenvalue < eigenSolver.eigenvalues()[i].real()) {
			if (fabs(largestEigenvalue - eigenSolver.eigenvalues()[i].real()) < 1e-6) {
				multiplicity++;
			} else {
				largestEigenvalue = eigenSolver.eigenvalues()[i].real();
				which = i;
				multiplicity = 1;
			}
		}
	}

	std::vector<double> lambda_1;
	for (size_t i = 0; i < N; i++) {
		double tmpelement = eigenSolver.eigenvectors().col(which)[i].real();
		lambda_1.push_back(tmpelement);
	}

	return lambda_1;
}

template<size_t N>
std::vector<double> Matrix<N>::getPrimalNormEigenvector() const {
	std::vector<double> v = getPrimalEigenvector();
	std::vector<double> tmp;
	for (size_t i = 0; i < N; i++) {
		tmp.push_back(v[i] / v[0]);
	}
	/*for (int i = 0; i < 4; i ++) std::cout << tmp[i] << " ";
	std::cout << std::endl; */
	return tmp;
}

template<size_t N>
void Matrix<N>::L1(std::vector<double> &v){
	double sum = 0;
	for (size_t i = 0; i < N; i++) sum += v[i];
	for (size_t i = 0; i < N; i++) v[i] /= sum;
}

template<size_t N>
bool Matrix<N>::testPrimalEigenvectorIsParetoOptimal() const {
	std::vector<double> lambda_1 = getPrimalNormEigenvector();

	return Matrix::testVectorParetoOptimal(lambda_1);
}

template<size_t N>
bool Matrix<N>::testVectorParetoOptimal(const std::vector<double> &w) const {
	if (Matrix::LPVectorParetoOptimal(w).isOptimal()) return true;
	return false;
}

template<size_t N>
bool Matrix<N>::testParetoOptimality(filterType filter) const {
	switch (filter) {
		case (filterType::EigenVectorMethod) :
			return Matrix::testPrimalEigenvectorIsParetoOptimal();
			break;
		case (filterType::AverageSpanTreeMethod) :
			return Matrix::testAvgSpanTreeParetoOptimal();
			break;
		case (filterType::CosineMethod) :
			return Matrix::testCosineParetoOptimal();
		default :
			throw "Invalid filterType was given to the pareto optimality tester function.\n";
			break;
		}
}

template<size_t N>
LpSolution<N> Matrix<N>::LPVectorParetoOptimal(const std::vector<double> &w) const {
	//Inputs: A matrix and w vector
		if (w.size() != N) throw "invalid length of vector!\n";

		if (verbosity) std::cout << "The input matrix: \n" << *this << std::endl;
		if (verbosity) {
			std::cout << "The input w vector: \n" << std::endl;
			for (size_t i = 0; i < w.size(); i++) 	std::cout << w[i] << "\t";
			std::cout << std::endl;
		}

		std::vector<spair> I,J;
		std::vector<double> v;
		std::vector<std::vector<double> > logm;

		if (verbosity) std::cout << "Pairs in set I:" << std::endl;
		for (size_t i = 0; i < N; i++) {
			v.push_back(log(w[i]));
			std::vector<double> logmrow;
			for (size_t j = 0; j < N; j++) {
				logmrow.push_back(log(get(i,j)));
				if ((w[i]/w[j] - get(i,j)) > pow(10.0,-8.0)) {
					spair tmp(i,j);
					if (verbosity) std::cout << i << " and " << j << "\n";
					I.push_back(tmp);
				}
			}
			logm.push_back(logmrow);
		}
		if (verbosity) {
			std::cout << "The initial v vector (log w): ";
			for (size_t i = 0; i < v.size(); i++) 	std::cout << v[i] << "\t";
			std::cout << "\n";
		}

		for (size_t i = 0; i < N; i++) {
			for (size_t j = i + 1; j < N; j++) {
				if (std::abs(w[i]/w[j] - get(i,j)) < pow(10.0,-8.0)) {
					spair tmp(i,j);
					J.push_back(tmp);
				}
			}
		}
		if (verbosity) std::cout << "Sizes of I and J: " << I.size() << " and " << J.size() << "\n";

		lemon::Lp LP;
		std::vector<lemon::Lp::Col> y; //y=log(x)
		std::vector<lemon::Lp::Col> s; //s=log(t)
		if (verbosity) std::cout << "The LP:\n";
		for (size_t i = 0; i < N; i++) { //y-k definialva
			lemon::Lp::Col tmp = LP.addCol();
			y.push_back(tmp);
		}
		for (size_t i = 0; i < I.size(); i++) { //s_ij
			lemon::Lp::Col tmp = LP.addCol();
			s.push_back(tmp);
			LP.colLowerBound(tmp, 0);
		}
		for (size_t i = 0; i < I.size(); i++) {
			lemon::Lp::Expr Ex;
			Ex += y[I[i].b];
			Ex -= y[I[i].a];
			LP.addRow(Ex <= (-logm[I[i].a][I[i].b])); //logm a b
			if (verbosity) std::cout << " y" << I[i].b << " -y" << I[i].a << "<=" << -logm[I[i].a][I[i].b] <<"\n";
			lemon::Lp::Expr Ex2;
			Ex2 -= y[I[i].b];
			Ex2 += y[I[i].a];
			Ex2 += s[i];
			LP.addRow(Ex2 <= (v[I[i].a] - v[I[i].b]));
			if (verbosity) std::cout << "-y" << I[i].b << " +y" << I[i].a << " +s" << i << " <=" << v[I[i].a] <<" - " << v[I[i].b] <<"\n";
		}
		//std::cout << "\n";
		for (size_t i = 0; i < J.size(); i++) {
			lemon::Lp::Expr Ex;
			Ex += y[J[i].a];
			Ex -= y[J[i].b];
			LP.addRow(Ex == logm[J[i].a][J[i].b]);
			if (verbosity) std::cout << "y" << J[i].a << " -y" << J[i].b << "=" << logm[J[i].a][J[i].b] <<"\n";
		}
		lemon::Lp::Expr exy;
		exy += y[0];
		LP.addRow(exy == 0);

		LP.min();
		lemon::Lp::Expr e;
		for (size_t i = 0; i < I.size(); i++) {
			e -= s[i];
		}
		LP.obj(e);
		LP.solve();
		if (verbosity) std::cout << "The primal type is: " << LP.primalType() << "\n";
		if (LP.primalType() == lemon::Lp::OPTIMAL) {
			if (verbosity) {
				std::cout << "Objective function value: " << LP.primal() << std::endl;
				for (size_t i = 0; i < y.size(); i++) std::cout << "y" << i << ": " << LP.primal(y[i]) << std::endl;
				for (size_t i = 0; i < y.size(); i++) std::cout << "x" << i << ": " << exp(LP.primal(y[i])) << std::endl;

				double sum = 0;
				for (size_t i = 0; i < y.size(); i++) sum += exp(LP.primal(y[i]));
				for (size_t i = 0; i < y.size(); i++) std::cout << "x_norm_" << i << ": " << exp(LP.primal(y[i]))/sum << std::endl;

				for (size_t i = 0; i < s.size(); i++) std::cout << "s" << i << ": " << LP.primal(s[i]) << std::endl;
				std::cout << LP.primal() << std::endl;
			}
			std::vector<double> xvector;
			std::vector<double> svector;

			for (size_t i = 0; i < y.size(); i++) xvector.push_back(exp(LP.primal(y[i])));
			for (size_t i = 0; i < s.size(); i++) svector.push_back(LP.primal(s[i]));

			LpSolution<N> ret(*this, w, I, J, LP.primal(), xvector, svector);
			return ret;
		} else {
			if (verbosity) {
				std::cout << "Optimal solution not found." << std::endl;
			}
			std::cout << "There is an error during making the LP. \nERROR INFO: " << LP.primalType() << "\n";
			throw "There is an error during making the LP.";
		}
}

template<size_t N>
std::vector<double> Matrix<N>::getCosineVector() const {
	//initialising b:
	std::vector<std::vector<double> > b;

	for (size_t i = 0; i < N; i++) {
		std::vector<double> tmp;
		for (size_t j = 0; j < N; j++) {
			tmp.push_back(0.0);
		}
		b.push_back(tmp);
	}

	//calculate the b_ij-s:
	for (size_t col = 0; col < N; col++) {
		double colSum = 0.0;
		for (size_t row = 0; row < N; row++) {
			colSum += get(row, col) * get(row, col);
		}
		colSum = sqrt(colSum);
		for (size_t row = 0; row < N; row++) {
			b[row][col] = Matrix::get(row, col) / colSum;
		}
	}

	//calculate formulas:
	std::vector<double> w;

	//nominators:
	for (size_t row = 0; row < N; row++) {
		double tmpsum = 0.0;
		for (size_t col = 0; col < N; col++) {
			tmpsum += b[row][col];
		}
		w.push_back(tmpsum);
	}

	//denominator:
	double fullSum = 0.0;
	for (size_t i = 0; i < N; i++) {
		for (size_t row = 0; row < N; row++) {
			double colSum = 0.0;
			for (size_t col = 0; col < N; col++) {
				colSum += b[row][col];
			}
			colSum *= colSum;
			fullSum += colSum;
		}
		fullSum = sqrt(fullSum);
	}

	for (size_t i = 0; i < w.size(); i++){
		w[i] /= fullSum;
	}
	return w;
}

template<size_t N>
bool Matrix<N>::testCosineParetoOptimal()const{
	return Matrix::testVectorParetoOptimal(getCosineVector());
}

template<size_t N>
std::vector<double> Matrix<N>::getMeanOfSpans()const{
	if (N != 4 && N != 5 && N != 6) throw ("Currently the average spantree method in 7 or higher dimensions is not implemented!\n");
	class VectorAdd {
	private:
		std::vector<double> n;
	public:
		VectorAdd() {
			std::vector<double> k;
			for (size_t i = 0; i < N; i++) {
				k.push_back(0.0);
			}
			n = k;
		}
		VectorAdd(std::vector<double> v):n(v){}
		void add(std::vector<double> v) {
			double sum = 0;
			for (size_t i = 0; i < N; i++) v[i] = exp(v[i]);
			for (size_t i = 0; i < N; i++) sum += v[i];
			for (size_t i = 0; i < N; i++) v[i] /= sum;
			for (size_t i = 0; i < N; i++) n[i] += v[i];
		}
		void divideBy(double x) {
			for (size_t i = 0; i < N; i++) n[i] /= x;
		}
		std::vector<double> getData() {
			return n;
		}
	} b;

	if (N == 4) {
		double a12,a13,a14,a23,a24,a34;
		a12 = log(get(0,1));
		a13 = log(get(0,2));
		a14 = log(get(0,3));
		a23 = log(get(1,2));
		a24 = log(get(1,3));
		a34 = log(get(2,3));

		b.add({0,a24-a14,a34-a14,-a14});
		b.add({0,a23+a34-a14,a34-a14,-a14});
		b.add({0,a24-a14,-a23+a24-a14,-a14});
		b.add({0,a24-a13-a34,-a13,-a13-a34});
		b.add({0,a23-a13,-a13,-a13-a34});
		b.add({0,a23-a13,-a13,a23-a24-a13});
		b.add({0,a24-a14,-a13,-a14});
		b.add({0,a23-a13,-a13,-a14});
		b.add({0,-a12,-a12-a24+a34,-a12-a24});
		b.add({0,-a12,-a12-a23,-a12-a23-a34});
		b.add({0,-a12,-a12-a23,-a12-a24});
		b.add({0,-a12,a34-a14,-a14});
		b.add({0,-a12,-a12-a23,-a14});
		b.add({0,-a12,-a13,-a13-a34});
		b.add({0,-a12,-a13,-a12-a24});
		b.add({0,-a12,-a13,-a14});
		b.divideBy(16.0);

		return b.getData();
	} else if (N == 5) {
		double a12,a13,a14,a15,a23,a24,a25,a34,a35,a45;
		a12 = log(get(0,1));
		a13 = log(get(0,2));
		a14 = log(get(0,3));
		a15 = log(get(0,4));
		a23 = log(get(1,2));
		a24 = log(get(1,3));
		a25 = log(get(1,4));
		a34 = log(get(2,3));
		a35 = log(get(2,4));
		a45 = log(get(3,4));

		b.add({0,a25-a15,a35-a15,a45-a15,-a15});
		b.add({0,a25-a15,a34+a45-a15,a45-a15,-a15});
		b.add({0,a25-a15,a35-a15,-a34+a35-a15,-a15});
		b.add({0,a24+a45-a15,a35-a15,a45-a15,-a15});
		b.add({0,a24+a45-a15,a34+a45-a15,a45-a15,-a15});
		b.add({0,a24-a34+a35-a15,a35-a15,-a34+a35-a15,-a15});
		b.add({0,a25-a15,a35-a15,-a24+a25-a15,-a15});
		b.add({0,a25-a15,-a24+a25+a34-a15,-a24+a25-a15,-a15});
		b.add({0,a23+a35-a15,a35-a15,a45-a15,-a15});
		b.add({0,a23+a34+a45-a15,a34+a45-a15,a45-a15,-a15});
		b.add({0,a23+a35-a15,a35-a15,-a34+a35-a15,-a15});
		b.add({0,a25-a15,-a23+a25-a15,a45-a15,-a15});
		b.add({0,a25-a15,-a23+a25-a15,-a23+a25-a34-a15,-a15});
		b.add({0,a24+a45-a15,-a23+a24+a45-a15,a45-a15,-a15});
		b.add({0,a23+a35-a15,a35-a15,a23-a24+a35-a15,-a15});
		b.add({0,a25-a15,-a23+a25-a15,-a24+a25-a15,-a15});
		b.add({0,a25-a14-a45,a35-a14-a45,-a14,-a14-a45});
		b.add({0,a25-a14-a45,a34-a14,-a14,-a14-a45});
		b.add({0,a25+a34-a35-a14,a34-a14,-a14,a34-a35-a14});
		b.add({0,a24-a14,a35-a14-a45,-a14,-a14-a45});
		b.add({0,a24-a14,a34-a14,-a14,-a14-a45});
		b.add({0,a24-a14,a34-a14,-a14,a34-a35-a14});
		b.add({0,a24-a14,a24-a25+a35-a14,-a14,a24-a25-a14});
		b.add({0,a24-a14,a34-a14,-a14,a24-a25-a14});
		b.add({0,a23+a35-a14-a45,a35-a14-a45,-a14,-a14-a45});
		b.add({0,a23+a34-a14,a34-a14,-a14,-a14-a45});
		b.add({0,a23+a34-a14,a34-a14,-a14,a34-a35-a14});
		b.add({0,a25-a14-a45,-a23+a25-a14-a45,-a14,-a14-a45});
		b.add({0,a23+a34-a14,a34-a14,-a14,a23-a25+a34-a14});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a14-a45});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a23+a24-a35-a14});
		b.add({0,a24-a14,-a23+a24-a14,-a14,a24-a25-a14});
		b.add({0,a25-a15,a35-a15,-a14,-a15});
		b.add({0,a25-a15,a34-a14,-a14,-a15});
		b.add({0,a24-a14,a35-a15,-a14,-a15});
		b.add({0,a24-a14,a34-a14,-a14,-a15});
		b.add({0,a23+a35-a15,a35-a15,-a14,-a15});
		b.add({0,a23+a34-a14,a34-a14,-a14,-a15});
		b.add({0,a25-a15,-a23+a25-a15,-a14,-a15});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a15});
		b.add({0,a25-a13-a35,-a13,-a13-a35+a45,-a13-a35});
		b.add({0,a25-a13-a34-a45,-a13,-a13-a34,-a13-a34-a45});
		b.add({0,a25-a13-a35,-a13,-a13-a34,-a13-a35});
		b.add({0,a24-a13-a35+a45,-a13,-a13-a35+a45,-a13-a35});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a13-a34-a45});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a13-a35});
		b.add({0,a25-a13-a35,-a13,-a24+a25-a13-a35,-a13-a35});
		b.add({0,a24-a13-a34,-a13,-a13-a34,a24-a25-a13-a34});
		b.add({0,a23-a13,-a13,-a13-a35+a45,-a13-a35});
		b.add({0,a23-a13,-a13,-a13-a34,-a13-a34-a45});
		b.add({0,a23-a13,-a13,-a13-a34,-a13-a35});
		b.add({0,a23-a13,-a13,a23-a25-a13+a45,a23-a25-a13});
		b.add({0,a23-a13,-a13,-a13-a34,a23-a25-a13});
		b.add({0,a23-a13,-a13,a23-a24-a13,a23-a24-a13-a45});
		b.add({0,a23-a13,-a13,a23-a24-a13,-a13-a35});
		b.add({0,a23-a13,-a13,a23-a24-a13,a23-a25-a13});
		b.add({0,a25-a15,-a13,a45-a15,-a15});
		b.add({0,a25-a15,-a13,-a13-a34,-a15});
		b.add({0,a24+a45-a15,-a13,a45-a15,-a15});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a15});
		b.add({0,a25-a15,-a13,-a24+a25-a15,-a15});
		b.add({0,a23-a13,-a13,a45-a15,-a15});
		b.add({0,a23-a13,-a13,-a13-a34,-a15});
		b.add({0,a23-a13,-a13,a23-a24-a13,-a15});
		b.add({0,a25-a14-a45,-a13,-a14,-a14-a45});
		b.add({0,a25-a13-a35,-a13,-a14,-a13-a35});
		b.add({0,a24-a14,-a13,-a14,-a14-a45});
		b.add({0,a24-a14,-a13,-a14,-a13-a35});
		b.add({0,a24-a14,-a13,-a14,a24-a25-a14});
		b.add({0,a23-a13,-a13,-a14,-a14-a45});
		b.add({0,a23-a13,-a13,-a14,-a13-a35});
		b.add({0,a23-a13,-a13,-a14,a23-a25-a13});
		b.add({0,a25-a15,-a13,-a14,-a15});
		b.add({0,a24-a14,-a13,-a14,-a15});
		b.add({0,a23-a13,-a13,-a14,-a15});
		b.add({0,-a12,-a12-a25+a35,-a12-a25+a45,-a12-a25});
		b.add({0,-a12,-a12-a25+a34+a45,-a12-a25+a45,-a12-a25});
		b.add({0,-a12,-a12-a25+a35,-a12-a25-a34+a35,-a12-a25});
		b.add({0,-a12,-a12-a24+a35-a45,-a12-a24,-a12-a24-a45});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a24-a45});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a24+a34-a35});
		b.add({0,-a12,-a12-a25+a35,-a12-a24,-a12-a25});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a25});
		b.add({0,-a12,-a12-a23,-a12-a23-a35+a45,-a12-a23-a35});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a23-a34-a45});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a23-a35});
		b.add({0,-a12,-a12-a23,-a12-a25+a45,-a12-a25});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a25});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a24-a45});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a23-a35});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a25});
		b.add({0,-a12,a35-a15,a45-a15,-a15});
		b.add({0,-a12,a34+a45-a15,a45-a15,-a15});
		b.add({0,-a12,a35-a15,-a34+a35-a15,-a15});
		b.add({0,-a12,a35-a15,-a12-a24,-a15});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a15});
		b.add({0,-a12,-a12-a23,a45-a15,-a15});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a15});
		b.add({0,-a12,-a12-a23,-a12-a24,-a15});
		b.add({0,-a12,a35-a14-a45,-a14,-a14-a45});
		b.add({0,-a12,a34-a14,-a14,-a14-a45});
		b.add({0,-a12,a34-a14,-a14,a34-a35-a14});
		b.add({0,-a12,-a12-a25+a35,-a14,-a12-a25});
		b.add({0,-a12,a34-a14,-a14,-a12-a25});
		b.add({0,-a12,-a12-a23,-a14,-a14-a45});
		b.add({0,-a12,-a12-a23,-a14,-a12-a23-a35});
		b.add({0,-a12,-a12-a23,-a14,-a12-a25});
		b.add({0,-a12,a35-a15,-a14,-a15});
		b.add({0,-a12,a34-a14,-a14,-a15});
		b.add({0,-a12,-a12-a23,-a14,-a15});
		b.add({0,-a12,-a13,-a13-a35+a45,-a13-a35});
		b.add({0,-a12,-a13,-a13-a34,-a13-a34-a45});
		b.add({0,-a12,-a13,-a13-a34,-a13-a35});
		b.add({0,-a12,-a13,-a12-a25+a45,-a12-a25});
		b.add({0,-a12,-a13,-a13-a34,-a12-a25});
		b.add({0,-a12,-a13,-a12-a24,-a12-a24-a45});
		b.add({0,-a12,-a13,-a12-a24,-a13-a35});
		b.add({0,-a12,-a13,-a12-a24,-a12-a25});
		b.add({0,-a12,-a13,a45-a15,-a15});
		b.add({0,-a12,-a13,-a13-a34,-a15});
		b.add({0,-a12,-a13,-a12-a24,-a15});
		b.add({0,-a12,-a13,-a14,-a14-a45});
		b.add({0,-a12,-a13,-a14,-a13-a35});
		b.add({0,-a12,-a13,-a14,-a12-a25});
		b.add({0,-a12,-a13,-a14,-a15});

		b.divideBy(125.0);

		return b.getData();
	} else if (N == 6) {
		double a12,a13,a14,a15,a16,a23,a24,a25,a26,a34,a35,a36,a45,a46,a56;

		a12 = log(get(0,1));
		a13 = log(get(0,2));
		a14 = log(get(0,3));
		a15 = log(get(0,4));
		a16 = log(get(0,5));
		a23 = log(get(1,2));
		a24 = log(get(1,3));
		a25 = log(get(1,4));
		a26 = log(get(1,5));
		a34 = log(get(2,3));
		a35 = log(get(2,4));
		a36 = log(get(2,5));
		a45 = log(get(3,4));
		a46 = log(get(3,5));
		a56 = log(get(4,5));

		b.add({0,a26-a16,a36-a16,a46-a16,a56-a16,-a16});
		b.add({0,a26-a16,a36-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a26-a16,a36-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a26-a16,a35+a56-a16,a46-a16,a56-a16,-a16});
		b.add({0,a26-a16,a35+a56-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a26-a16,a35-a45+a46-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a26-a16,a36-a16,a46-a16,-a35+a36-a16,-a16});
		b.add({0,a26-a16,a36-a16,-a35+a36+a45-a16,-a35+a36-a16,-a16});
		b.add({0,a26-a16,a34+a46-a16,a46-a16,a56-a16,-a16});
		b.add({0,a26-a16,a34+a45+a56-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a26-a16,a34+a46-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a26-a16,a36-a16,-a34+a36-a16,a56-a16,-a16});
		b.add({0,a26-a16,a36-a16,-a34+a36-a16,-a34+a36-a45-a16,-a16});
		b.add({0,a26-a16,a35+a56-a16,-a34+a35+a56-a16,a56-a16,-a16});
		b.add({0,a26-a16,a34+a46-a16,a46-a16,a34-a35+a46-a16,-a16});
		b.add({0,a26-a16,a36-a16,-a34+a36-a16,-a35+a36-a16,-a16});
		b.add({0,a25+a56-a16,a36-a16,a46-a16,a56-a16,-a16});
		b.add({0,a25+a56-a16,a36-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a25-a45+a46-a16,a36-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a25+a56-a16,a35+a56-a16,a46-a16,a56-a16,-a16});
		b.add({0,a25+a56-a16,a35+a56-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a25-a45+a46-a16,a35-a45+a46-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a25-a35+a36-a16,a36-a16,a46-a16,-a35+a36-a16,-a16});
		b.add({0,a25-a35+a36-a16,a36-a16,-a35+a36+a45-a16,-a35+a36-a16,-a16});
		b.add({0,a25+a56-a16,a34+a46-a16,a46-a16,a56-a16,-a16});
		b.add({0,a25+a56-a16,a34+a45+a56-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a25-a45+a46-a16,a34+a46-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a25+a56-a16,a36-a16,-a34+a36-a16,a56-a16,-a16});
		b.add({0,a25-a34+a36-a45-a16,a36-a16,-a34+a36-a16,-a34+a36-a45-a16,-a16});
		b.add({0,a25+a56-a16,a35+a56-a16,-a34+a35+a56-a16,a56-a16,-a16});
		b.add({0,a25+a34-a35+a46-a16,a34+a46-a16,a46-a16,a34-a35+a46-a16,-a16});
		b.add({0,a25-a35+a36-a16,a36-a16,-a34+a36-a16,-a35+a36-a16,-a16});
		b.add({0,a26-a16,a36-a16,a46-a16,-a25+a26-a16,-a16});
		b.add({0,a26-a16,a36-a16,-a25+a26+a45-a16,-a25+a26-a16,-a16});
		b.add({0,a26-a16,-a25+a26+a35-a16,a46-a16,-a25+a26-a16,-a16});
		b.add({0,a26-a16,-a25+a26+a35-a16,-a25+a26+a45-a16,-a25+a26-a16,-a16});
		b.add({0,a26-a16,a34+a46-a16,a46-a16,-a25+a26-a16,-a16});
		b.add({0,a26-a16,-a25+a26+a34+a45-a16,-a25+a26+a45-a16,-a25+a26-a16,-a16});
		b.add({0,a26-a16,a36-a16,-a34+a36-a16,-a25+a26-a16,-a16});
		b.add({0,a26-a16,-a25+a26+a35-a16,-a25+a26-a34+a35-a16,-a25+a26-a16,-a16});
		b.add({0,a24+a46-a16,a36-a16,a46-a16,a56-a16,-a16});
		b.add({0,a24+a45+a56-a16,a36-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a24+a46-a16,a36-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a24+a46-a16,a35+a56-a16,a46-a16,a56-a16,-a16});
		b.add({0,a24+a45+a56-a16,a35+a56-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a24+a46-a16,a35-a45+a46-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a24+a46-a16,a36-a16,a46-a16,-a35+a36-a16,-a16});
		b.add({0,a24-a35+a36+a45-a16,a36-a16,-a35+a36+a45-a16,-a35+a36-a16,-a16});
		b.add({0,a24+a46-a16,a34+a46-a16,a46-a16,a56-a16,-a16});
		b.add({0,a24+a45+a56-a16,a34+a45+a56-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a24+a46-a16,a34+a46-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a24-a34+a36-a16,a36-a16,-a34+a36-a16,a56-a16,-a16});
		b.add({0,a24-a34+a36-a16,a36-a16,-a34+a36-a16,-a34+a36-a45-a16,-a16});
		b.add({0,a24-a34+a35+a56-a16,a35+a56-a16,-a34+a35+a56-a16,a56-a16,-a16});
		b.add({0,a24+a46-a16,a34+a46-a16,a46-a16,a34-a35+a46-a16,-a16});
		b.add({0,a24-a34+a36-a16,a36-a16,-a34+a36-a16,-a35+a36-a16,-a16});
		b.add({0,a26-a16,a36-a16,-a24+a26-a16,a56-a16,-a16});
		b.add({0,a26-a16,a36-a16,-a24+a26-a16,-a24+a26-a45-a16,-a16});
		b.add({0,a26-a16,a35+a56-a16,-a24+a26-a16,a56-a16,-a16});
		b.add({0,a26-a16,-a24+a26+a35-a45-a16,-a24+a26-a16,-a24+a26-a45-a16,-a16});
		b.add({0,a26-a16,a36-a16,-a24+a26-a16,-a35+a36-a16,-a16});
		b.add({0,a26-a16,-a24+a26+a34-a16,-a24+a26-a16,a56-a16,-a16});
		b.add({0,a26-a16,-a24+a26+a34-a16,-a24+a26-a16,-a24+a26-a45-a16,-a16});
		b.add({0,a26-a16,-a24+a26+a34-a16,-a24+a26-a16,-a24+a26+a34-a35-a16,-a16});
		b.add({0,a25+a56-a16,a36-a16,-a24+a25+a56-a16,a56-a16,-a16});
		b.add({0,a24+a46-a16,a36-a16,a46-a16,a24-a25+a46-a16,-a16});
		b.add({0,a25+a56-a16,a35+a56-a16,-a24+a25+a56-a16,a56-a16,-a16});
		b.add({0,a24+a46-a16,a24-a25+a35+a46-a16,a46-a16,a24-a25+a46-a16,-a16});
		b.add({0,a25-a35+a36-a16,a36-a16,-a24+a25-a35+a36-a16,-a35+a36-a16,-a16});
		b.add({0,a25+a56-a16,-a24+a25+a34+a56-a16,-a24+a25+a56-a16,a56-a16,-a16});
		b.add({0,a24+a46-a16,a34+a46-a16,a46-a16,a24-a25+a46-a16,-a16});
		b.add({0,a24-a34+a36-a16,a36-a16,-a34+a36-a16,a24-a25-a34+a36-a16,-a16});
		b.add({0,a26-a16,a36-a16,-a24+a26-a16,-a25+a26-a16,-a16});
		b.add({0,a26-a16,-a25+a26+a35-a16,-a24+a26-a16,-a25+a26-a16,-a16});
		b.add({0,a26-a16,-a24+a26+a34-a16,-a24+a26-a16,-a25+a26-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,a46-a16,a56-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a23+a35+a56-a16,a35+a56-a16,a46-a16,a56-a16,-a16});
		b.add({0,a23+a35+a56-a16,a35+a56-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a23+a35-a45+a46-a16,a35-a45+a46-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,a46-a16,-a35+a36-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,-a35+a36+a45-a16,-a35+a36-a16,-a16});
		b.add({0,a23+a34+a46-a16,a34+a46-a16,a46-a16,a56-a16,-a16});
		b.add({0,a23+a34+a45+a56-a16,a34+a45+a56-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a23+a34+a46-a16,a34+a46-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,-a34+a36-a16,a56-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,-a34+a36-a16,-a34+a36-a45-a16,-a16});
		b.add({0,a23+a35+a56-a16,a35+a56-a16,-a34+a35+a56-a16,a56-a16,-a16});
		b.add({0,a23+a34+a46-a16,a34+a46-a16,a46-a16,a34-a35+a46-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,-a34+a36-a16,-a35+a36-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,a46-a16,a56-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,a46-a16,-a23+a26-a35-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a23+a26-a35+a45-a16,-a23+a26-a35-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a23+a26-a34-a16,a56-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a23+a26-a34-a16,-a23+a26-a34-a45-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a23+a26-a34-a16,-a23+a26-a35-a16,-a16});
		b.add({0,a25+a56-a16,-a23+a25+a56-a16,a46-a16,a56-a16,-a16});
		b.add({0,a25+a56-a16,-a23+a25+a56-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a25-a45+a46-a16,-a23+a25-a45+a46-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,a46-a16,a23-a25+a36-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,a23-a25+a36+a45-a16,a23-a25+a36-a16,-a16});
		b.add({0,a25+a56-a16,-a23+a25+a56-a16,-a23+a25-a34+a56-a16,a56-a16,-a16});
		b.add({0,a23+a34+a46-a16,a34+a46-a16,a46-a16,a23-a25+a34+a46-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,-a34+a36-a16,a23-a25+a36-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,a46-a16,-a25+a26-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a25+a26+a45-a16,-a25+a26-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a23+a26-a34-a16,-a25+a26-a16,-a16});
		b.add({0,a24+a46-a16,-a23+a24+a46-a16,a46-a16,a56-a16,-a16});
		b.add({0,a24+a45+a56-a16,-a23+a24+a45+a56-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,a24+a46-a16,-a23+a24+a46-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,a23-a24+a36-a16,a56-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,a23-a24+a36-a16,a23-a24+a36-a45-a16,-a16});
		b.add({0,a23+a35+a56-a16,a35+a56-a16,a23-a24+a35+a56-a16,a56-a16,-a16});
		b.add({0,a24+a46-a16,-a23+a24+a46-a16,a46-a16,-a23+a24-a35+a46-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,a23-a24+a36-a16,-a35+a36-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a24+a26-a16,a56-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a24+a26-a16,-a24+a26-a45-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a24+a26-a16,-a23+a26-a35-a16,-a16});
		b.add({0,a25+a56-a16,-a23+a25+a56-a16,-a24+a25+a56-a16,a56-a16,-a16});
		b.add({0,a24+a46-a16,-a23+a24+a46-a16,a46-a16,a24-a25+a46-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,a23-a24+a36-a16,a23-a25+a36-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a24+a26-a16,-a25+a26-a16,-a16});
		b.add({0,a26-a15-a56,a36-a15-a56,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a26-a15-a56,a36-a15-a56,a45-a15,-a15,-a15-a56});
		b.add({0,a26+a45-a46-a15,a36+a45-a46-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a26-a15-a56,a35-a15,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a26-a15-a56,a35-a15,a45-a15,-a15,-a15-a56});
		b.add({0,a26+a45-a46-a15,a35-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a26+a35-a36-a15,a35-a15,a35-a36+a46-a15,-a15,a35-a36-a15});
		b.add({0,a26+a35-a36-a15,a35-a15,a45-a15,-a15,a35-a36-a15});
		b.add({0,a26-a15-a56,a34+a46-a15-a56,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a26-a15-a56,a34+a45-a15,a45-a15,-a15,-a15-a56});
		b.add({0,a26+a45-a46-a15,a34+a45-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a26-a15-a56,a36-a15-a56,-a34+a36-a15-a56,-a15,-a15-a56});
		b.add({0,a26+a34-a36+a45-a15,a34+a45-a15,a45-a15,-a15,a34-a36+a45-a15});
		b.add({0,a26-a15-a56,a35-a15,-a34+a35-a15,-a15,-a15-a56});
		b.add({0,a26-a34+a35-a46-a15,a35-a15,-a34+a35-a15,-a15,-a34+a35-a46-a15});
		b.add({0,a26+a35-a36-a15,a35-a15,-a34+a35-a15,-a15,a35-a36-a15});
		b.add({0,a25-a15,a36-a15-a56,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a25-a15,a36-a15-a56,a45-a15,-a15,-a15-a56});
		b.add({0,a25-a15,a36+a45-a46-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a25-a15,a35-a15,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a25-a15,a35-a15,a45-a15,-a15,-a15-a56});
		b.add({0,a25-a15,a35-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a25-a15,a35-a15,a35-a36+a46-a15,-a15,a35-a36-a15});
		b.add({0,a25-a15,a35-a15,a45-a15,-a15,a35-a36-a15});
		b.add({0,a25-a15,a34+a46-a15-a56,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a25-a15,a34+a45-a15,a45-a15,-a15,-a15-a56});
		b.add({0,a25-a15,a34+a45-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a25-a15,a36-a15-a56,-a34+a36-a15-a56,-a15,-a15-a56});
		b.add({0,a25-a15,a34+a45-a15,a45-a15,-a15,a34-a36+a45-a15});
		b.add({0,a25-a15,a35-a15,-a34+a35-a15,-a15,-a15-a56});
		b.add({0,a25-a15,a35-a15,-a34+a35-a15,-a15,-a34+a35-a46-a15});
		b.add({0,a25-a15,a35-a15,-a34+a35-a15,-a15,a35-a36-a15});
		b.add({0,a25-a15,a25-a26+a36-a15,a25-a26+a46-a15,-a15,a25-a26-a15});
		b.add({0,a25-a15,a25-a26+a36-a15,a45-a15,-a15,a25-a26-a15});
		b.add({0,a25-a15,a35-a15,a25-a26+a46-a15,-a15,a25-a26-a15});
		b.add({0,a25-a15,a35-a15,a45-a15,-a15,a25-a26-a15});
		b.add({0,a25-a15,a25-a26+a34+a46-a15,a25-a26+a46-a15,-a15,a25-a26-a15});
		b.add({0,a25-a15,a34+a45-a15,a45-a15,-a15,a25-a26-a15});
		b.add({0,a25-a15,a25-a26+a36-a15,a25-a26-a34+a36-a15,-a15,a25-a26-a15});
		b.add({0,a25-a15,a35-a15,-a34+a35-a15,-a15,a25-a26-a15});
		b.add({0,a24+a46-a15-a56,a36-a15-a56,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a24+a45-a15,a36-a15-a56,a45-a15,-a15,-a15-a56});
		b.add({0,a24+a45-a15,a36+a45-a46-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a24+a46-a15-a56,a35-a15,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a24+a45-a15,a35-a15,a45-a15,-a15,-a15-a56});
		b.add({0,a24+a45-a15,a35-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a24+a35-a36+a46-a15,a35-a15,a35-a36+a46-a15,-a15,a35-a36-a15});
		b.add({0,a24+a45-a15,a35-a15,a45-a15,-a15,a35-a36-a15});
		b.add({0,a24+a46-a15-a56,a34+a46-a15-a56,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a24+a45-a15,a34+a45-a15,a45-a15,-a15,-a15-a56});
		b.add({0,a24+a45-a15,a34+a45-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a24-a34+a36-a15-a56,a36-a15-a56,-a34+a36-a15-a56,-a15,-a15-a56});
		b.add({0,a24+a45-a15,a34+a45-a15,a45-a15,-a15,a34-a36+a45-a15});
		b.add({0,a24-a34+a35-a15,a35-a15,-a34+a35-a15,-a15,-a15-a56});
		b.add({0,a24-a34+a35-a15,a35-a15,-a34+a35-a15,-a15,-a34+a35-a46-a15});
		b.add({0,a24-a34+a35-a15,a35-a15,-a34+a35-a15,-a15,a35-a36-a15});
		b.add({0,a26-a15-a56,a36-a15-a56,-a24+a26-a15-a56,-a15,-a15-a56});
		b.add({0,a24+a45-a15,a24-a26+a36+a45-a15,a45-a15,-a15,a24-a26+a45-a15});
		b.add({0,a26-a15-a56,a35-a15,-a24+a26-a15-a56,-a15,-a15-a56});
		b.add({0,a24+a45-a15,a35-a15,a45-a15,-a15,a24-a26+a45-a15});
		b.add({0,a26+a35-a36-a15,a35-a15,-a24+a26+a35-a36-a15,-a15,a35-a36-a15});
		b.add({0,a26-a15-a56,-a24+a26+a34-a15-a56,-a24+a26-a15-a56,-a15,-a15-a56});
		b.add({0,a24+a45-a15,a34+a45-a15,a45-a15,-a15,a24-a26+a45-a15});
		b.add({0,a24-a34+a35-a15,a35-a15,-a34+a35-a15,-a15,a24-a26-a34+a35-a15});
		b.add({0,a25-a15,a36-a15-a56,-a24+a25-a15,-a15,-a15-a56});
		b.add({0,a25-a15,-a24+a25+a36-a46-a15,-a24+a25-a15,-a15,-a24+a25-a46-a15});
		b.add({0,a25-a15,a35-a15,-a24+a25-a15,-a15,-a15-a56});
		b.add({0,a25-a15,a35-a15,-a24+a25-a15,-a15,-a24+a25-a46-a15});
		b.add({0,a25-a15,a35-a15,-a24+a25-a15,-a15,a35-a36-a15});
		b.add({0,a25-a15,-a24+a25+a34-a15,-a24+a25-a15,-a15,-a15-a56});
		b.add({0,a25-a15,-a24+a25+a34-a15,-a24+a25-a15,-a15,-a24+a25-a46-a15});
		b.add({0,a25-a15,-a24+a25+a34-a15,-a24+a25-a15,-a15,-a24+a25+a34-a36-a15});
		b.add({0,a25-a15,a25-a26+a36-a15,-a24+a25-a15,-a15,a25-a26-a15});
		b.add({0,a25-a15,a35-a15,-a24+a25-a15,-a15,a25-a26-a15});
		b.add({0,a25-a15,-a24+a25+a34-a15,-a24+a25-a15,-a15,a25-a26-a15});
		b.add({0,a23+a36-a15-a56,a36-a15-a56,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a23+a36-a15-a56,a36-a15-a56,a45-a15,-a15,-a15-a56});
		b.add({0,a23+a36+a45-a46-a15,a36+a45-a46-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a23+a35-a15,a35-a15,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a23+a35-a15,a35-a15,a45-a15,-a15,-a15-a56});
		b.add({0,a23+a35-a15,a35-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a23+a35-a15,a35-a15,a35-a36+a46-a15,-a15,a35-a36-a15});
		b.add({0,a23+a35-a15,a35-a15,a45-a15,-a15,a35-a36-a15});
		b.add({0,a23+a34+a46-a15-a56,a34+a46-a15-a56,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a23+a34+a45-a15,a34+a45-a15,a45-a15,-a15,-a15-a56});
		b.add({0,a23+a34+a45-a15,a34+a45-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a23+a36-a15-a56,a36-a15-a56,-a34+a36-a15-a56,-a15,-a15-a56});
		b.add({0,a23+a34+a45-a15,a34+a45-a15,a45-a15,-a15,a34-a36+a45-a15});
		b.add({0,a23+a35-a15,a35-a15,-a34+a35-a15,-a15,-a15-a56});
		b.add({0,a23+a35-a15,a35-a15,-a34+a35-a15,-a15,-a34+a35-a46-a15});
		b.add({0,a23+a35-a15,a35-a15,-a34+a35-a15,-a15,a35-a36-a15});
		b.add({0,a26-a15-a56,-a23+a26-a15-a56,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a26-a15-a56,-a23+a26-a15-a56,a45-a15,-a15,-a15-a56});
		b.add({0,a26+a45-a46-a15,-a23+a26+a45-a46-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a23+a35-a15,a35-a15,a23-a26+a35+a46-a15,-a15,a23-a26+a35-a15});
		b.add({0,a23+a35-a15,a35-a15,a45-a15,-a15,a23-a26+a35-a15});
		b.add({0,a26-a15-a56,-a23+a26-a15-a56,-a23+a26-a34-a15-a56,-a15,-a15-a56});
		b.add({0,a23+a34+a45-a15,a34+a45-a15,a45-a15,-a15,a23-a26+a34+a45-a15});
		b.add({0,a23+a35-a15,a35-a15,-a34+a35-a15,-a15,a23-a26+a35-a15});
		b.add({0,a25-a15,-a23+a25-a15,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a25-a15,-a23+a25-a15,a45-a15,-a15,-a15-a56});
		b.add({0,a25-a15,-a23+a25-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a25-a15,-a23+a25-a15,-a23+a25-a36+a46-a15,-a15,-a23+a25-a36-a15});
		b.add({0,a25-a15,-a23+a25-a15,a45-a15,-a15,-a23+a25-a36-a15});
		b.add({0,a25-a15,-a23+a25-a15,-a23+a25-a34-a15,-a15,-a15-a56});
		b.add({0,a25-a15,-a23+a25-a15,-a23+a25-a34-a15,-a15,-a23+a25-a34-a46-a15});
		b.add({0,a25-a15,-a23+a25-a15,-a23+a25-a34-a15,-a15,-a23+a25-a36-a15});
		b.add({0,a25-a15,-a23+a25-a15,a25-a26+a46-a15,-a15,a25-a26-a15});
		b.add({0,a25-a15,-a23+a25-a15,a45-a15,-a15,a25-a26-a15});
		b.add({0,a25-a15,-a23+a25-a15,-a23+a25-a34-a15,-a15,a25-a26-a15});
		b.add({0,a24+a46-a15-a56,-a23+a24+a46-a15-a56,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a24+a45-a15,-a23+a24+a45-a15,a45-a15,-a15,-a15-a56});
		b.add({0,a24+a45-a15,-a23+a24+a45-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,a23+a36-a15-a56,a36-a15-a56,a23-a24+a36-a15-a56,-a15,-a15-a56});
		b.add({0,a24+a45-a15,-a23+a24+a45-a15,a45-a15,-a15,-a23+a24-a36+a45-a15});
		b.add({0,a23+a35-a15,a35-a15,a23-a24+a35-a15,-a15,-a15-a56});
		b.add({0,a23+a35-a15,a35-a15,a23-a24+a35-a15,-a15,a23-a24+a35-a46-a15});
		b.add({0,a23+a35-a15,a35-a15,a23-a24+a35-a15,-a15,a35-a36-a15});
		b.add({0,a26-a15-a56,-a23+a26-a15-a56,-a24+a26-a15-a56,-a15,-a15-a56});
		b.add({0,a24+a45-a15,-a23+a24+a45-a15,a45-a15,-a15,a24-a26+a45-a15});
		b.add({0,a23+a35-a15,a35-a15,a23-a24+a35-a15,-a15,a23-a26+a35-a15});
		b.add({0,a25-a15,-a23+a25-a15,-a24+a25-a15,-a15,-a15-a56});
		b.add({0,a25-a15,-a23+a25-a15,-a24+a25-a15,-a15,-a24+a25-a46-a15});
		b.add({0,a25-a15,-a23+a25-a15,-a24+a25-a15,-a15,-a23+a25-a36-a15});
		b.add({0,a25-a15,-a23+a25-a15,-a24+a25-a15,-a15,a25-a26-a15});
		b.add({0,a26-a16,a36-a16,a46-a16,-a15,-a16});
		b.add({0,a26-a16,a36-a16,a45-a15,-a15,-a16});
		b.add({0,a26-a16,a35-a15,a46-a16,-a15,-a16});
		b.add({0,a26-a16,a35-a15,a45-a15,-a15,-a16});
		b.add({0,a26-a16,a34+a46-a16,a46-a16,-a15,-a16});
		b.add({0,a26-a16,a34+a45-a15,a45-a15,-a15,-a16});
		b.add({0,a26-a16,a36-a16,-a34+a36-a16,-a15,-a16});
		b.add({0,a26-a16,a35-a15,-a34+a35-a15,-a15,-a16});
		b.add({0,a25-a15,a36-a16,a46-a16,-a15,-a16});
		b.add({0,a25-a15,a36-a16,a45-a15,-a15,-a16});
		b.add({0,a25-a15,a35-a15,a46-a16,-a15,-a16});
		b.add({0,a25-a15,a35-a15,a45-a15,-a15,-a16});
		b.add({0,a25-a15,a34+a46-a16,a46-a16,-a15,-a16});
		b.add({0,a25-a15,a34+a45-a15,a45-a15,-a15,-a16});
		b.add({0,a25-a15,a36-a16,-a34+a36-a16,-a15,-a16});
		b.add({0,a25-a15,a35-a15,-a34+a35-a15,-a15,-a16});
		b.add({0,a24+a46-a16,a36-a16,a46-a16,-a15,-a16});
		b.add({0,a24+a45-a15,a36-a16,a45-a15,-a15,-a16});
		b.add({0,a24+a46-a16,a35-a15,a46-a16,-a15,-a16});
		b.add({0,a24+a45-a15,a35-a15,a45-a15,-a15,-a16});
		b.add({0,a24+a46-a16,a34+a46-a16,a46-a16,-a15,-a16});
		b.add({0,a24+a45-a15,a34+a45-a15,a45-a15,-a15,-a16});
		b.add({0,a24-a34+a36-a16,a36-a16,-a34+a36-a16,-a15,-a16});
		b.add({0,a24-a34+a35-a15,a35-a15,-a34+a35-a15,-a15,-a16});
		b.add({0,a26-a16,a36-a16,-a24+a26-a16,-a15,-a16});
		b.add({0,a26-a16,a35-a15,-a24+a26-a16,-a15,-a16});
		b.add({0,a26-a16,-a24+a26+a34-a16,-a24+a26-a16,-a15,-a16});
		b.add({0,a25-a15,a36-a16,-a24+a25-a15,-a15,-a16});
		b.add({0,a25-a15,a35-a15,-a24+a25-a15,-a15,-a16});
		b.add({0,a25-a15,-a24+a25+a34-a15,-a24+a25-a15,-a15,-a16});
		b.add({0,a23+a36-a16,a36-a16,a46-a16,-a15,-a16});
		b.add({0,a23+a36-a16,a36-a16,a45-a15,-a15,-a16});
		b.add({0,a23+a35-a15,a35-a15,a46-a16,-a15,-a16});
		b.add({0,a23+a35-a15,a35-a15,a45-a15,-a15,-a16});
		b.add({0,a23+a34+a46-a16,a34+a46-a16,a46-a16,-a15,-a16});
		b.add({0,a23+a34+a45-a15,a34+a45-a15,a45-a15,-a15,-a16});
		b.add({0,a23+a36-a16,a36-a16,-a34+a36-a16,-a15,-a16});
		b.add({0,a23+a35-a15,a35-a15,-a34+a35-a15,-a15,-a16});
		b.add({0,a26-a16,-a23+a26-a16,a46-a16,-a15,-a16});
		b.add({0,a26-a16,-a23+a26-a16,a45-a15,-a15,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a23+a26-a34-a16,-a15,-a16});
		b.add({0,a25-a15,-a23+a25-a15,a46-a16,-a15,-a16});
		b.add({0,a25-a15,-a23+a25-a15,a45-a15,-a15,-a16});
		b.add({0,a25-a15,-a23+a25-a15,-a23+a25-a34-a15,-a15,-a16});
		b.add({0,a24+a46-a16,-a23+a24+a46-a16,a46-a16,-a15,-a16});
		b.add({0,a24+a45-a15,-a23+a24+a45-a15,a45-a15,-a15,-a16});
		b.add({0,a23+a36-a16,a36-a16,a23-a24+a36-a16,-a15,-a16});
		b.add({0,a23+a35-a15,a35-a15,a23-a24+a35-a15,-a15,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a24+a26-a16,-a15,-a16});
		b.add({0,a25-a15,-a23+a25-a15,-a24+a25-a15,-a15,-a16});
		b.add({0,a26-a14-a46,a36-a14-a46,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a26-a14-a45-a56,a36-a14-a45-a56,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a26-a14-a46,a36-a14-a46,-a14,-a14-a45,-a14-a46});
		b.add({0,a26-a14-a46,a35-a14-a46+a56,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a26-a14-a45-a56,a35-a14-a45,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a26-a14-a46,a35-a14-a45,-a14,-a14-a45,-a14-a46});
		b.add({0,a26-a14-a46,a36-a14-a46,-a14,-a35+a36-a14-a46,-a14-a46});
		b.add({0,a26+a35-a36-a14-a45,a35-a14-a45,-a14,-a14-a45,a35-a36-a14-a45});
		b.add({0,a26-a14-a46,a34-a14,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a26-a14-a45-a56,a34-a14,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a26-a14-a46,a34-a14,-a14,-a14-a45,-a14-a46});
		b.add({0,a26+a34-a36-a14,a34-a14,-a14,a34-a36-a14+a56,a34-a36-a14});
		b.add({0,a26+a34-a36-a14,a34-a14,-a14,-a14-a45,a34-a36-a14});
		b.add({0,a26+a34-a35-a14-a56,a34-a14,-a14,a34-a35-a14,a34-a35-a14-a56});
		b.add({0,a26-a14-a46,a34-a14,-a14,a34-a35-a14,-a14-a46});
		b.add({0,a26+a34-a36-a14,a34-a14,-a14,a34-a35-a14,a34-a36-a14});
		b.add({0,a25-a14-a46+a56,a36-a14-a46,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a25-a14-a45,a36-a14-a45-a56,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a25-a14-a45,a36-a14-a46,-a14,-a14-a45,-a14-a46});
		b.add({0,a25-a14-a46+a56,a35-a14-a46+a56,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a25-a14-a45,a35-a14-a45,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a25-a14-a45,a35-a14-a45,-a14,-a14-a45,-a14-a46});
		b.add({0,a25-a35+a36-a14-a46,a36-a14-a46,-a14,-a35+a36-a14-a46,-a14-a46});
		b.add({0,a25-a14-a45,a35-a14-a45,-a14,-a14-a45,a35-a36-a14-a45});
		b.add({0,a25-a14-a46+a56,a34-a14,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a25-a14-a45,a34-a14,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a25-a14-a45,a34-a14,-a14,-a14-a45,-a14-a46});
		b.add({0,a25+a34-a36-a14+a56,a34-a14,-a14,a34-a36-a14+a56,a34-a36-a14});
		b.add({0,a25-a14-a45,a34-a14,-a14,-a14-a45,a34-a36-a14});
		b.add({0,a25+a34-a35-a14,a34-a14,-a14,a34-a35-a14,a34-a35-a14-a56});
		b.add({0,a25+a34-a35-a14,a34-a14,-a14,a34-a35-a14,-a14-a46});
		b.add({0,a25+a34-a35-a14,a34-a14,-a14,a34-a35-a14,a34-a36-a14});
		b.add({0,a26-a14-a46,a36-a14-a46,-a14,-a25+a26-a14-a46,-a14-a46});
		b.add({0,a25-a14-a45,a25-a26+a36-a14-a45,-a14,-a14-a45,a25-a26-a14-a45});
		b.add({0,a26-a14-a46,-a25+a26+a35-a14-a46,-a14,-a25+a26-a14-a46,-a14-a46});
		b.add({0,a25-a14-a45,a35-a14-a45,-a14,-a14-a45,a25-a26-a14-a45});
		b.add({0,a26-a14-a46,a34-a14,-a14,-a25+a26-a14-a46,-a14-a46});
		b.add({0,a25-a14-a45,a34-a14,-a14,-a14-a45,a25-a26-a14-a45});
		b.add({0,a26+a34-a36-a14,a34-a14,-a14,-a25+a26+a34-a36-a14,a34-a36-a14});
		b.add({0,a25+a34-a35-a14,a34-a14,-a14,a34-a35-a14,a25-a26+a34-a35-a14});
		b.add({0,a24-a14,a36-a14-a46,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a24-a14,a36-a14-a45-a56,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a24-a14,a36-a14-a46,-a14,-a14-a45,-a14-a46});
		b.add({0,a24-a14,a35-a14-a46+a56,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a24-a14,a35-a14-a45,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a24-a14,a35-a14-a45,-a14,-a14-a45,-a14-a46});
		b.add({0,a24-a14,a36-a14-a46,-a14,-a35+a36-a14-a46,-a14-a46});
		b.add({0,a24-a14,a35-a14-a45,-a14,-a14-a45,a35-a36-a14-a45});
		b.add({0,a24-a14,a34-a14,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a24-a14,a34-a14,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a24-a14,a34-a14,-a14,-a14-a45,-a14-a46});
		b.add({0,a24-a14,a34-a14,-a14,a34-a36-a14+a56,a34-a36-a14});
		b.add({0,a24-a14,a34-a14,-a14,-a14-a45,a34-a36-a14});
		b.add({0,a24-a14,a34-a14,-a14,a34-a35-a14,a34-a35-a14-a56});
		b.add({0,a24-a14,a34-a14,-a14,a34-a35-a14,-a14-a46});
		b.add({0,a24-a14,a34-a14,-a14,a34-a35-a14,a34-a36-a14});
		b.add({0,a24-a14,a24-a26+a36-a14,-a14,a24-a26-a14+a56,a24-a26-a14});
		b.add({0,a24-a14,a24-a26+a36-a14,-a14,-a14-a45,a24-a26-a14});
		b.add({0,a24-a14,a24-a26+a35-a14+a56,-a14,a24-a26-a14+a56,a24-a26-a14});
		b.add({0,a24-a14,a35-a14-a45,-a14,-a14-a45,a24-a26-a14});
		b.add({0,a24-a14,a24-a26+a36-a14,-a14,a24-a26-a35+a36-a14,a24-a26-a14});
		b.add({0,a24-a14,a34-a14,-a14,a24-a26-a14+a56,a24-a26-a14});
		b.add({0,a24-a14,a34-a14,-a14,-a14-a45,a24-a26-a14});
		b.add({0,a24-a14,a34-a14,-a14,a34-a35-a14,a24-a26-a14});
		b.add({0,a24-a14,a24-a25+a36-a14-a56,-a14,a24-a25-a14,a24-a25-a14-a56});
		b.add({0,a24-a14,a36-a14-a46,-a14,a24-a25-a14,-a14-a46});
		b.add({0,a24-a14,a24-a25+a35-a14,-a14,a24-a25-a14,a24-a25-a14-a56});
		b.add({0,a24-a14,a24-a25+a35-a14,-a14,a24-a25-a14,-a14-a46});
		b.add({0,a24-a14,a24-a25+a35-a14,-a14,a24-a25-a14,a24-a25+a35-a36-a14});
		b.add({0,a24-a14,a34-a14,-a14,a24-a25-a14,a24-a25-a14-a56});
		b.add({0,a24-a14,a34-a14,-a14,a24-a25-a14,-a14-a46});
		b.add({0,a24-a14,a34-a14,-a14,a24-a25-a14,a34-a36-a14});
		b.add({0,a24-a14,a24-a26+a36-a14,-a14,a24-a25-a14,a24-a26-a14});
		b.add({0,a24-a14,a24-a25+a35-a14,-a14,a24-a25-a14,a24-a26-a14});
		b.add({0,a24-a14,a34-a14,-a14,a24-a25-a14,a24-a26-a14});
		b.add({0,a23+a36-a14-a46,a36-a14-a46,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a23+a36-a14-a45-a56,a36-a14-a45-a56,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a23+a36-a14-a46,a36-a14-a46,-a14,-a14-a45,-a14-a46});
		b.add({0,a23+a35-a14-a46+a56,a35-a14-a46+a56,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a23+a35-a14-a45,a35-a14-a45,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a23+a35-a14-a45,a35-a14-a45,-a14,-a14-a45,-a14-a46});
		b.add({0,a23+a36-a14-a46,a36-a14-a46,-a14,-a35+a36-a14-a46,-a14-a46});
		b.add({0,a23+a35-a14-a45,a35-a14-a45,-a14,-a14-a45,a35-a36-a14-a45});
		b.add({0,a23+a34-a14,a34-a14,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a23+a34-a14,a34-a14,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a23+a34-a14,a34-a14,-a14,-a14-a45,-a14-a46});
		b.add({0,a23+a34-a14,a34-a14,-a14,a34-a36-a14+a56,a34-a36-a14});
		b.add({0,a23+a34-a14,a34-a14,-a14,-a14-a45,a34-a36-a14});
		b.add({0,a23+a34-a14,a34-a14,-a14,a34-a35-a14,a34-a35-a14-a56});
		b.add({0,a23+a34-a14,a34-a14,-a14,a34-a35-a14,-a14-a46});
		b.add({0,a23+a34-a14,a34-a14,-a14,a34-a35-a14,a34-a36-a14});
		b.add({0,a26-a14-a46,-a23+a26-a14-a46,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a26-a14-a45-a56,-a23+a26-a14-a45-a56,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a26-a14-a46,-a23+a26-a14-a46,-a14,-a14-a45,-a14-a46});
		b.add({0,a26-a14-a46,-a23+a26-a14-a46,-a14,-a23+a26-a35-a14-a46,-a14-a46});
		b.add({0,a23+a35-a14-a45,a35-a14-a45,-a14,-a14-a45,a23-a26+a35-a14-a45});
		b.add({0,a23+a34-a14,a34-a14,-a14,a23-a26+a34-a14+a56,a23-a26+a34-a14});
		b.add({0,a23+a34-a14,a34-a14,-a14,-a14-a45,a23-a26+a34-a14});
		b.add({0,a23+a34-a14,a34-a14,-a14,a34-a35-a14,a23-a26+a34-a14});
		b.add({0,a25-a14-a46+a56,-a23+a25-a14-a46+a56,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a25-a14-a45,-a23+a25-a14-a45,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a25-a14-a45,-a23+a25-a14-a45,-a14,-a14-a45,-a14-a46});
		b.add({0,a23+a36-a14-a46,a36-a14-a46,-a14,a23-a25+a36-a14-a46,-a14-a46});
		b.add({0,a25-a14-a45,-a23+a25-a14-a45,-a14,-a14-a45,-a23+a25-a36-a14-a45});
		b.add({0,a23+a34-a14,a34-a14,-a14,a23-a25+a34-a14,a23-a25+a34-a14-a56});
		b.add({0,a23+a34-a14,a34-a14,-a14,a23-a25+a34-a14,-a14-a46});
		b.add({0,a23+a34-a14,a34-a14,-a14,a23-a25+a34-a14,a34-a36-a14});
		b.add({0,a26-a14-a46,-a23+a26-a14-a46,-a14,-a25+a26-a14-a46,-a14-a46});
		b.add({0,a25-a14-a45,-a23+a25-a14-a45,-a14,-a14-a45,a25-a26-a14-a45});
		b.add({0,a23+a34-a14,a34-a14,-a14,a23-a25+a34-a14,a23-a26+a34-a14});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a14-a45,-a14-a46});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a23+a24-a36-a14+a56,-a23+a24-a36-a14});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a14-a45,-a23+a24-a36-a14});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a23+a24-a35-a14,-a23+a24-a35-a14-a56});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a23+a24-a35-a14,-a14-a46});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a23+a24-a35-a14,-a23+a24-a36-a14});
		b.add({0,a24-a14,-a23+a24-a14,-a14,a24-a26-a14+a56,a24-a26-a14});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a14-a45,a24-a26-a14});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a23+a24-a35-a14,a24-a26-a14});
		b.add({0,a24-a14,-a23+a24-a14,-a14,a24-a25-a14,a24-a25-a14-a56});
		b.add({0,a24-a14,-a23+a24-a14,-a14,a24-a25-a14,-a14-a46});
		b.add({0,a24-a14,-a23+a24-a14,-a14,a24-a25-a14,-a23+a24-a36-a14});
		b.add({0,a24-a14,-a23+a24-a14,-a14,a24-a25-a14,a24-a26-a14});
		b.add({0,a26-a16,a36-a16,-a14,a56-a16,-a16});
		b.add({0,a26-a16,a36-a16,-a14,-a14-a45,-a16});
		b.add({0,a26-a16,a35+a56-a16,-a14,a56-a16,-a16});
		b.add({0,a26-a16,a35-a14-a45,-a14,-a14-a45,-a16});
		b.add({0,a26-a16,a36-a16,-a14,-a35+a36-a16,-a16});
		b.add({0,a26-a16,a34-a14,-a14,a56-a16,-a16});
		b.add({0,a26-a16,a34-a14,-a14,-a14-a45,-a16});
		b.add({0,a26-a16,a34-a14,-a14,a34-a35-a14,-a16});
		b.add({0,a25+a56-a16,a36-a16,-a14,a56-a16,-a16});
		b.add({0,a25-a14-a45,a36-a16,-a14,-a14-a45,-a16});
		b.add({0,a25+a56-a16,a35+a56-a16,-a14,a56-a16,-a16});
		b.add({0,a25-a14-a45,a35-a14-a45,-a14,-a14-a45,-a16});
		b.add({0,a25-a35+a36-a16,a36-a16,-a14,-a35+a36-a16,-a16});
		b.add({0,a25+a56-a16,a34-a14,-a14,a56-a16,-a16});
		b.add({0,a25-a14-a45,a34-a14,-a14,-a14-a45,-a16});
		b.add({0,a25+a34-a35-a14,a34-a14,-a14,a34-a35-a14,-a16});
		b.add({0,a26-a16,a36-a16,-a14,-a25+a26-a16,-a16});
		b.add({0,a26-a16,-a25+a26+a35-a16,-a14,-a25+a26-a16,-a16});
		b.add({0,a26-a16,a34-a14,-a14,-a25+a26-a16,-a16});
		b.add({0,a24-a14,a36-a16,-a14,a56-a16,-a16});
		b.add({0,a24-a14,a36-a16,-a14,-a14-a45,-a16});
		b.add({0,a24-a14,a35+a56-a16,-a14,a56-a16,-a16});
		b.add({0,a24-a14,a35-a14-a45,-a14,-a14-a45,-a16});
		b.add({0,a24-a14,a36-a16,-a14,-a35+a36-a16,-a16});
		b.add({0,a24-a14,a34-a14,-a14,a56-a16,-a16});
		b.add({0,a24-a14,a34-a14,-a14,-a14-a45,-a16});
		b.add({0,a24-a14,a34-a14,-a14,a34-a35-a14,-a16});
		b.add({0,a24-a14,a36-a16,-a14,a24-a25-a14,-a16});
		b.add({0,a24-a14,a24-a25+a35-a14,-a14,a24-a25-a14,-a16});
		b.add({0,a24-a14,a34-a14,-a14,a24-a25-a14,-a16});
		b.add({0,a23+a36-a16,a36-a16,-a14,a56-a16,-a16});
		b.add({0,a23+a36-a16,a36-a16,-a14,-a14-a45,-a16});
		b.add({0,a23+a35+a56-a16,a35+a56-a16,-a14,a56-a16,-a16});
		b.add({0,a23+a35-a14-a45,a35-a14-a45,-a14,-a14-a45,-a16});
		b.add({0,a23+a36-a16,a36-a16,-a14,-a35+a36-a16,-a16});
		b.add({0,a23+a34-a14,a34-a14,-a14,a56-a16,-a16});
		b.add({0,a23+a34-a14,a34-a14,-a14,-a14-a45,-a16});
		b.add({0,a23+a34-a14,a34-a14,-a14,a34-a35-a14,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a14,a56-a16,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a14,-a14-a45,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a14,-a23+a26-a35-a16,-a16});
		b.add({0,a25+a56-a16,-a23+a25+a56-a16,-a14,a56-a16,-a16});
		b.add({0,a25-a14-a45,-a23+a25-a14-a45,-a14,-a14-a45,-a16});
		b.add({0,a23+a36-a16,a36-a16,-a14,a23-a25+a36-a16,-a16});
		b.add({0,a23+a34-a14,a34-a14,-a14,a23-a25+a34-a14,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a14,-a25+a26-a16,-a16});
		b.add({0,a24-a14,-a23+a24-a14,-a14,a56-a16,-a16});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a14-a45,-a16});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a23+a24-a35-a14,-a16});
		b.add({0,a24-a14,-a23+a24-a14,-a14,a24-a25-a14,-a16});
		b.add({0,a26-a15-a56,a36-a15-a56,-a14,-a15,-a15-a56});
		b.add({0,a26-a14-a46,a36-a14-a46,-a14,-a15,-a14-a46});
		b.add({0,a26-a15-a56,a35-a15,-a14,-a15,-a15-a56});
		b.add({0,a26-a14-a46,a35-a15,-a14,-a15,-a14-a46});
		b.add({0,a26+a35-a36-a15,a35-a15,-a14,-a15,a35-a36-a15});
		b.add({0,a26-a15-a56,a34-a14,-a14,-a15,-a15-a56});
		b.add({0,a26-a14-a46,a34-a14,-a14,-a15,-a14-a46});
		b.add({0,a26+a34-a36-a14,a34-a14,-a14,-a15,a34-a36-a14});
		b.add({0,a25-a15,a36-a15-a56,-a14,-a15,-a15-a56});
		b.add({0,a25-a15,a36-a14-a46,-a14,-a15,-a14-a46});
		b.add({0,a25-a15,a35-a15,-a14,-a15,-a15-a56});
		b.add({0,a25-a15,a35-a15,-a14,-a15,-a14-a46});
		b.add({0,a25-a15,a35-a15,-a14,-a15,a35-a36-a15});
		b.add({0,a25-a15,a34-a14,-a14,-a15,-a15-a56});
		b.add({0,a25-a15,a34-a14,-a14,-a15,-a14-a46});
		b.add({0,a25-a15,a34-a14,-a14,-a15,a34-a36-a14});
		b.add({0,a25-a15,a25-a26+a36-a15,-a14,-a15,a25-a26-a15});
		b.add({0,a25-a15,a35-a15,-a14,-a15,a25-a26-a15});
		b.add({0,a25-a15,a34-a14,-a14,-a15,a25-a26-a15});
		b.add({0,a24-a14,a36-a15-a56,-a14,-a15,-a15-a56});
		b.add({0,a24-a14,a36-a14-a46,-a14,-a15,-a14-a46});
		b.add({0,a24-a14,a35-a15,-a14,-a15,-a15-a56});
		b.add({0,a24-a14,a35-a15,-a14,-a15,-a14-a46});
		b.add({0,a24-a14,a35-a15,-a14,-a15,a35-a36-a15});
		b.add({0,a24-a14,a34-a14,-a14,-a15,-a15-a56});
		b.add({0,a24-a14,a34-a14,-a14,-a15,-a14-a46});
		b.add({0,a24-a14,a34-a14,-a14,-a15,a34-a36-a14});
		b.add({0,a24-a14,a24-a26+a36-a14,-a14,-a15,a24-a26-a14});
		b.add({0,a24-a14,a35-a15,-a14,-a15,a24-a26-a14});
		b.add({0,a24-a14,a34-a14,-a14,-a15,a24-a26-a14});
		b.add({0,a23+a36-a15-a56,a36-a15-a56,-a14,-a15,-a15-a56});
		b.add({0,a23+a36-a14-a46,a36-a14-a46,-a14,-a15,-a14-a46});
		b.add({0,a23+a35-a15,a35-a15,-a14,-a15,-a15-a56});
		b.add({0,a23+a35-a15,a35-a15,-a14,-a15,-a14-a46});
		b.add({0,a23+a35-a15,a35-a15,-a14,-a15,a35-a36-a15});
		b.add({0,a23+a34-a14,a34-a14,-a14,-a15,-a15-a56});
		b.add({0,a23+a34-a14,a34-a14,-a14,-a15,-a14-a46});
		b.add({0,a23+a34-a14,a34-a14,-a14,-a15,a34-a36-a14});
		b.add({0,a26-a15-a56,-a23+a26-a15-a56,-a14,-a15,-a15-a56});
		b.add({0,a26-a14-a46,-a23+a26-a14-a46,-a14,-a15,-a14-a46});
		b.add({0,a23+a35-a15,a35-a15,-a14,-a15,a23-a26+a35-a15});
		b.add({0,a23+a34-a14,a34-a14,-a14,-a15,a23-a26+a34-a14});
		b.add({0,a25-a15,-a23+a25-a15,-a14,-a15,-a15-a56});
		b.add({0,a25-a15,-a23+a25-a15,-a14,-a15,-a14-a46});
		b.add({0,a25-a15,-a23+a25-a15,-a14,-a15,-a23+a25-a36-a15});
		b.add({0,a25-a15,-a23+a25-a15,-a14,-a15,a25-a26-a15});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a15,-a15-a56});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a15,-a14-a46});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a15,-a23+a24-a36-a14});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a15,a24-a26-a14});
		b.add({0,a26-a16,a36-a16,-a14,-a15,-a16});
		b.add({0,a26-a16,a35-a15,-a14,-a15,-a16});
		b.add({0,a26-a16,a34-a14,-a14,-a15,-a16});
		b.add({0,a25-a15,a36-a16,-a14,-a15,-a16});
		b.add({0,a25-a15,a35-a15,-a14,-a15,-a16});
		b.add({0,a25-a15,a34-a14,-a14,-a15,-a16});
		b.add({0,a24-a14,a36-a16,-a14,-a15,-a16});
		b.add({0,a24-a14,a35-a15,-a14,-a15,-a16});
		b.add({0,a24-a14,a34-a14,-a14,-a15,-a16});
		b.add({0,a23+a36-a16,a36-a16,-a14,-a15,-a16});
		b.add({0,a23+a35-a15,a35-a15,-a14,-a15,-a16});
		b.add({0,a23+a34-a14,a34-a14,-a14,-a15,-a16});
		b.add({0,a26-a16,-a23+a26-a16,-a14,-a15,-a16});
		b.add({0,a25-a15,-a23+a25-a15,-a14,-a15,-a16});
		b.add({0,a24-a14,-a23+a24-a14,-a14,-a15,-a16});
		b.add({0,a26-a13-a36,-a13,-a13-a36+a46,-a13-a36+a56,-a13-a36});
		b.add({0,a26-a13-a36,-a13,-a13-a36+a45+a56,-a13-a36+a56,-a13-a36});
		b.add({0,a26-a13-a36,-a13,-a13-a36+a46,-a13-a36-a45+a46,-a13-a36});
		b.add({0,a26-a13-a35-a56,-a13,-a13-a35+a46-a56,-a13-a35,-a13-a35-a56});
		b.add({0,a26-a13-a35-a56,-a13,-a13-a35+a45,-a13-a35,-a13-a35-a56});
		b.add({0,a26-a13-a35+a45-a46,-a13,-a13-a35+a45,-a13-a35,-a13-a35+a45-a46});
		b.add({0,a26-a13-a36,-a13,-a13-a36+a46,-a13-a35,-a13-a36});
		b.add({0,a26-a13-a36,-a13,-a13-a35+a45,-a13-a35,-a13-a36});
		b.add({0,a26-a13-a34-a46,-a13,-a13-a34,-a13-a34-a46+a56,-a13-a34-a46});
		b.add({0,a26-a13-a34-a45-a56,-a13,-a13-a34,-a13-a34-a45,-a13-a34-a45-a56});
		b.add({0,a26-a13-a34-a46,-a13,-a13-a34,-a13-a34-a45,-a13-a34-a46});
		b.add({0,a26-a13-a36,-a13,-a13-a34,-a13-a36+a56,-a13-a36});
		b.add({0,a26-a13-a36,-a13,-a13-a34,-a13-a34-a45,-a13-a36});
		b.add({0,a26-a13-a35-a56,-a13,-a13-a34,-a13-a35,-a13-a35-a56});
		b.add({0,a26-a13-a34-a46,-a13,-a13-a34,-a13-a35,-a13-a34-a46});
		b.add({0,a26-a13-a36,-a13,-a13-a34,-a13-a35,-a13-a36});
		b.add({0,a25-a13-a36+a56,-a13,-a13-a36+a46,-a13-a36+a56,-a13-a36});
		b.add({0,a25-a13-a36+a56,-a13,-a13-a36+a45+a56,-a13-a36+a56,-a13-a36});
		b.add({0,a25-a13-a36-a45+a46,-a13,-a13-a36+a46,-a13-a36-a45+a46,-a13-a36});
		b.add({0,a25-a13-a35,-a13,-a13-a35+a46-a56,-a13-a35,-a13-a35-a56});
		b.add({0,a25-a13-a35,-a13,-a13-a35+a45,-a13-a35,-a13-a35-a56});
		b.add({0,a25-a13-a35,-a13,-a13-a35+a45,-a13-a35,-a13-a35+a45-a46});
		b.add({0,a25-a13-a35,-a13,-a13-a36+a46,-a13-a35,-a13-a36});
		b.add({0,a25-a13-a35,-a13,-a13-a35+a45,-a13-a35,-a13-a36});
		b.add({0,a25-a13-a34-a46+a56,-a13,-a13-a34,-a13-a34-a46+a56,-a13-a34-a46});
		b.add({0,a25-a13-a34-a45,-a13,-a13-a34,-a13-a34-a45,-a13-a34-a45-a56});
		b.add({0,a25-a13-a34-a45,-a13,-a13-a34,-a13-a34-a45,-a13-a34-a46});
		b.add({0,a25-a13-a36+a56,-a13,-a13-a34,-a13-a36+a56,-a13-a36});
		b.add({0,a25-a13-a34-a45,-a13,-a13-a34,-a13-a34-a45,-a13-a36});
		b.add({0,a25-a13-a35,-a13,-a13-a34,-a13-a35,-a13-a35-a56});
		b.add({0,a25-a13-a35,-a13,-a13-a34,-a13-a35,-a13-a34-a46});
		b.add({0,a25-a13-a35,-a13,-a13-a34,-a13-a35,-a13-a36});
		b.add({0,a26-a13-a36,-a13,-a13-a36+a46,-a25+a26-a13-a36,-a13-a36});
		b.add({0,a26-a13-a36,-a13,-a25+a26-a13-a36+a45,-a25+a26-a13-a36,-a13-a36});
		b.add({0,a25-a13-a35,-a13,a25-a26-a13-a35+a46,-a13-a35,a25-a26-a13-a35});
		b.add({0,a25-a13-a35,-a13,-a13-a35+a45,-a13-a35,a25-a26-a13-a35});
		b.add({0,a26-a13-a34-a46,-a13,-a13-a34,-a25+a26-a13-a34-a46,-a13-a34-a46});
		b.add({0,a25-a13-a34-a45,-a13,-a13-a34,-a13-a34-a45,a25-a26-a13-a34-a45});
		b.add({0,a26-a13-a36,-a13,-a13-a34,-a25+a26-a13-a36,-a13-a36});
		b.add({0,a25-a13-a35,-a13,-a13-a34,-a13-a35,a25-a26-a13-a35});
		b.add({0,a24-a13-a36+a46,-a13,-a13-a36+a46,-a13-a36+a56,-a13-a36});
		b.add({0,a24-a13-a36+a45+a56,-a13,-a13-a36+a45+a56,-a13-a36+a56,-a13-a36});
		b.add({0,a24-a13-a36+a46,-a13,-a13-a36+a46,-a13-a36-a45+a46,-a13-a36});
		b.add({0,a24-a13-a35+a46-a56,-a13,-a13-a35+a46-a56,-a13-a35,-a13-a35-a56});
		b.add({0,a24-a13-a35+a45,-a13,-a13-a35+a45,-a13-a35,-a13-a35-a56});
		b.add({0,a24-a13-a35+a45,-a13,-a13-a35+a45,-a13-a35,-a13-a35+a45-a46});
		b.add({0,a24-a13-a36+a46,-a13,-a13-a36+a46,-a13-a35,-a13-a36});
		b.add({0,a24-a13-a35+a45,-a13,-a13-a35+a45,-a13-a35,-a13-a36});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a13-a34-a46+a56,-a13-a34-a46});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a13-a34-a45,-a13-a34-a45-a56});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a13-a34-a45,-a13-a34-a46});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a13-a36+a56,-a13-a36});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a13-a34-a45,-a13-a36});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a13-a35,-a13-a35-a56});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a13-a35,-a13-a34-a46});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a13-a35,-a13-a36});
		b.add({0,a26-a13-a36,-a13,-a24+a26-a13-a36,-a13-a36+a56,-a13-a36});
		b.add({0,a26-a13-a36,-a13,-a24+a26-a13-a36,-a24+a26-a13-a36-a45,-a13-a36});
		b.add({0,a26-a13-a35-a56,-a13,-a24+a26-a13-a35-a56,-a13-a35,-a13-a35-a56});
		b.add({0,a24-a13-a35+a45,-a13,-a13-a35+a45,-a13-a35,a24-a26-a13-a35+a45});
		b.add({0,a26-a13-a36,-a13,-a24+a26-a13-a36,-a13-a35,-a13-a36});
		b.add({0,a24-a13-a34,-a13,-a13-a34,a24-a26-a13-a34+a56,a24-a26-a13-a34});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a13-a34-a45,a24-a26-a13-a34});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a13-a35,a24-a26-a13-a34});
		b.add({0,a25-a13-a36+a56,-a13,-a24+a25-a13-a36+a56,-a13-a36+a56,-a13-a36});
		b.add({0,a24-a13-a36+a46,-a13,-a13-a36+a46,a24-a25-a13-a36+a46,-a13-a36});
		b.add({0,a25-a13-a35,-a13,-a24+a25-a13-a35,-a13-a35,-a13-a35-a56});
		b.add({0,a25-a13-a35,-a13,-a24+a25-a13-a35,-a13-a35,-a24+a25-a13-a35-a46});
		b.add({0,a25-a13-a35,-a13,-a24+a25-a13-a35,-a13-a35,-a13-a36});
		b.add({0,a24-a13-a34,-a13,-a13-a34,a24-a25-a13-a34,a24-a25-a13-a34-a56});
		b.add({0,a24-a13-a34,-a13,-a13-a34,a24-a25-a13-a34,-a13-a34-a46});
		b.add({0,a24-a13-a34,-a13,-a13-a34,a24-a25-a13-a34,-a13-a36});
		b.add({0,a26-a13-a36,-a13,-a24+a26-a13-a36,-a25+a26-a13-a36,-a13-a36});
		b.add({0,a25-a13-a35,-a13,-a24+a25-a13-a35,-a13-a35,a25-a26-a13-a35});
		b.add({0,a24-a13-a34,-a13,-a13-a34,a24-a25-a13-a34,a24-a26-a13-a34});
		b.add({0,a23-a13,-a13,-a13-a36+a46,-a13-a36+a56,-a13-a36});
		b.add({0,a23-a13,-a13,-a13-a36+a45+a56,-a13-a36+a56,-a13-a36});
		b.add({0,a23-a13,-a13,-a13-a36+a46,-a13-a36-a45+a46,-a13-a36});
		b.add({0,a23-a13,-a13,-a13-a35+a46-a56,-a13-a35,-a13-a35-a56});
		b.add({0,a23-a13,-a13,-a13-a35+a45,-a13-a35,-a13-a35-a56});
		b.add({0,a23-a13,-a13,-a13-a35+a45,-a13-a35,-a13-a35+a45-a46});
		b.add({0,a23-a13,-a13,-a13-a36+a46,-a13-a35,-a13-a36});
		b.add({0,a23-a13,-a13,-a13-a35+a45,-a13-a35,-a13-a36});
		b.add({0,a23-a13,-a13,-a13-a34,-a13-a34-a46+a56,-a13-a34-a46});
		b.add({0,a23-a13,-a13,-a13-a34,-a13-a34-a45,-a13-a34-a45-a56});
		b.add({0,a23-a13,-a13,-a13-a34,-a13-a34-a45,-a13-a34-a46});
		b.add({0,a23-a13,-a13,-a13-a34,-a13-a36+a56,-a13-a36});
		b.add({0,a23-a13,-a13,-a13-a34,-a13-a34-a45,-a13-a36});
		b.add({0,a23-a13,-a13,-a13-a34,-a13-a35,-a13-a35-a56});
		b.add({0,a23-a13,-a13,-a13-a34,-a13-a35,-a13-a34-a46});
		b.add({0,a23-a13,-a13,-a13-a34,-a13-a35,-a13-a36});
		b.add({0,a23-a13,-a13,a23-a26-a13+a46,a23-a26-a13+a56,a23-a26-a13});
		b.add({0,a23-a13,-a13,a23-a26-a13+a45+a56,a23-a26-a13+a56,a23-a26-a13});
		b.add({0,a23-a13,-a13,a23-a26-a13+a46,a23-a26-a13-a45+a46,a23-a26-a13});
		b.add({0,a23-a13,-a13,a23-a26-a13+a46,-a13-a35,a23-a26-a13});
		b.add({0,a23-a13,-a13,-a13-a35+a45,-a13-a35,a23-a26-a13});
		b.add({0,a23-a13,-a13,-a13-a34,a23-a26-a13+a56,a23-a26-a13});
		b.add({0,a23-a13,-a13,-a13-a34,-a13-a34-a45,a23-a26-a13});
		b.add({0,a23-a13,-a13,-a13-a34,-a13-a35,a23-a26-a13});
		b.add({0,a23-a13,-a13,a23-a25-a13+a46-a56,a23-a25-a13,a23-a25-a13-a56});
		b.add({0,a23-a13,-a13,a23-a25-a13+a45,a23-a25-a13,a23-a25-a13-a56});
		b.add({0,a23-a13,-a13,a23-a25-a13+a45,a23-a25-a13,a23-a25-a13+a45-a46});
		b.add({0,a23-a13,-a13,-a13-a36+a46,a23-a25-a13,-a13-a36});
		b.add({0,a23-a13,-a13,a23-a25-a13+a45,a23-a25-a13,-a13-a36});
		b.add({0,a23-a13,-a13,-a13-a34,a23-a25-a13,a23-a25-a13-a56});
		b.add({0,a23-a13,-a13,-a13-a34,a23-a25-a13,-a13-a34-a46});
		b.add({0,a23-a13,-a13,-a13-a34,a23-a25-a13,-a13-a36});
		b.add({0,a23-a13,-a13,a23-a26-a13+a46,a23-a25-a13,a23-a26-a13});
		b.add({0,a23-a13,-a13,a23-a25-a13+a45,a23-a25-a13,a23-a26-a13});
		b.add({0,a23-a13,-a13,-a13-a34,a23-a25-a13,a23-a26-a13});
		b.add({0,a23-a13,-a13,a23-a24-a13,a23-a24-a13-a46+a56,a23-a24-a13-a46});
		b.add({0,a23-a13,-a13,a23-a24-a13,a23-a24-a13-a45,a23-a24-a13-a45-a56});
		b.add({0,a23-a13,-a13,a23-a24-a13,a23-a24-a13-a45,a23-a24-a13-a46});
		b.add({0,a23-a13,-a13,a23-a24-a13,-a13-a36+a56,-a13-a36});
		b.add({0,a23-a13,-a13,a23-a24-a13,a23-a24-a13-a45,-a13-a36});
		b.add({0,a23-a13,-a13,a23-a24-a13,-a13-a35,-a13-a35-a56});
		b.add({0,a23-a13,-a13,a23-a24-a13,-a13-a35,a23-a24-a13-a46});
		b.add({0,a23-a13,-a13,a23-a24-a13,-a13-a35,-a13-a36});
		b.add({0,a23-a13,-a13,a23-a24-a13,a23-a26-a13+a56,a23-a26-a13});
		b.add({0,a23-a13,-a13,a23-a24-a13,a23-a24-a13-a45,a23-a26-a13});
		b.add({0,a23-a13,-a13,a23-a24-a13,-a13-a35,a23-a26-a13});
		b.add({0,a23-a13,-a13,a23-a24-a13,a23-a25-a13,a23-a25-a13-a56});
		b.add({0,a23-a13,-a13,a23-a24-a13,a23-a25-a13,a23-a24-a13-a46});
		b.add({0,a23-a13,-a13,a23-a24-a13,a23-a25-a13,-a13-a36});
		b.add({0,a23-a13,-a13,a23-a24-a13,a23-a25-a13,a23-a26-a13});
		b.add({0,a26-a16,-a13,a46-a16,a56-a16,-a16});
		b.add({0,a26-a16,-a13,a45+a56-a16,a56-a16,-a16});
		b.add({0,a26-a16,-a13,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a26-a16,-a13,a46-a16,-a13-a35,-a16});
		b.add({0,a26-a16,-a13,-a13-a35+a45,-a13-a35,-a16});
		b.add({0,a26-a16,-a13,-a13-a34,a56-a16,-a16});
		b.add({0,a26-a16,-a13,-a13-a34,-a13-a34-a45,-a16});
		b.add({0,a26-a16,-a13,-a13-a34,-a13-a35,-a16});
		b.add({0,a25+a56-a16,-a13,a46-a16,a56-a16,-a16});
		b.add({0,a25+a56-a16,-a13,a45+a56-a16,a56-a16,-a16});
		b.add({0,a25-a45+a46-a16,-a13,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a25-a13-a35,-a13,a46-a16,-a13-a35,-a16});
		b.add({0,a25-a13-a35,-a13,-a13-a35+a45,-a13-a35,-a16});
		b.add({0,a25+a56-a16,-a13,-a13-a34,a56-a16,-a16});
		b.add({0,a25-a13-a34-a45,-a13,-a13-a34,-a13-a34-a45,-a16});
		b.add({0,a25-a13-a35,-a13,-a13-a34,-a13-a35,-a16});
		b.add({0,a26-a16,-a13,a46-a16,-a25+a26-a16,-a16});
		b.add({0,a26-a16,-a13,-a25+a26+a45-a16,-a25+a26-a16,-a16});
		b.add({0,a26-a16,-a13,-a13-a34,-a25+a26-a16,-a16});
		b.add({0,a24+a46-a16,-a13,a46-a16,a56-a16,-a16});
		b.add({0,a24+a45+a56-a16,-a13,a45+a56-a16,a56-a16,-a16});
		b.add({0,a24+a46-a16,-a13,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a24+a46-a16,-a13,a46-a16,-a13-a35,-a16});
		b.add({0,a24-a13-a35+a45,-a13,-a13-a35+a45,-a13-a35,-a16});
		b.add({0,a24-a13-a34,-a13,-a13-a34,a56-a16,-a16});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a13-a34-a45,-a16});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a13-a35,-a16});
		b.add({0,a26-a16,-a13,-a24+a26-a16,a56-a16,-a16});
		b.add({0,a26-a16,-a13,-a24+a26-a16,-a24+a26-a45-a16,-a16});
		b.add({0,a26-a16,-a13,-a24+a26-a16,-a13-a35,-a16});
		b.add({0,a25+a56-a16,-a13,-a24+a25+a56-a16,a56-a16,-a16});
		b.add({0,a24+a46-a16,-a13,a46-a16,a24-a25+a46-a16,-a16});
		b.add({0,a25-a13-a35,-a13,-a24+a25-a13-a35,-a13-a35,-a16});
		b.add({0,a24-a13-a34,-a13,-a13-a34,a24-a25-a13-a34,-a16});
		b.add({0,a26-a16,-a13,-a24+a26-a16,-a25+a26-a16,-a16});
		b.add({0,a23-a13,-a13,a46-a16,a56-a16,-a16});
		b.add({0,a23-a13,-a13,a45+a56-a16,a56-a16,-a16});
		b.add({0,a23-a13,-a13,a46-a16,-a45+a46-a16,-a16});
		b.add({0,a23-a13,-a13,a46-a16,-a13-a35,-a16});
		b.add({0,a23-a13,-a13,-a13-a35+a45,-a13-a35,-a16});
		b.add({0,a23-a13,-a13,-a13-a34,a56-a16,-a16});
		b.add({0,a23-a13,-a13,-a13-a34,-a13-a34-a45,-a16});
		b.add({0,a23-a13,-a13,-a13-a34,-a13-a35,-a16});
		b.add({0,a23-a13,-a13,a46-a16,a23-a25-a13,-a16});
		b.add({0,a23-a13,-a13,a23-a25-a13+a45,a23-a25-a13,-a16});
		b.add({0,a23-a13,-a13,-a13-a34,a23-a25-a13,-a16});
		b.add({0,a23-a13,-a13,a23-a24-a13,a56-a16,-a16});
		b.add({0,a23-a13,-a13,a23-a24-a13,a23-a24-a13-a45,-a16});
		b.add({0,a23-a13,-a13,a23-a24-a13,-a13-a35,-a16});
		b.add({0,a23-a13,-a13,a23-a24-a13,a23-a25-a13,-a16});
		b.add({0,a26-a15-a56,-a13,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a26-a15-a56,-a13,a45-a15,-a15,-a15-a56});
		b.add({0,a26+a45-a46-a15,-a13,a45-a15,-a15,a45-a46-a15});
		b.add({0,a26-a13-a36,-a13,-a13-a36+a46,-a15,-a13-a36});
		b.add({0,a26-a13-a36,-a13,a45-a15,-a15,-a13-a36});
		b.add({0,a26-a15-a56,-a13,-a13-a34,-a15,-a15-a56});
		b.add({0,a26-a13-a34-a46,-a13,-a13-a34,-a15,-a13-a34-a46});
		b.add({0,a26-a13-a36,-a13,-a13-a34,-a15,-a13-a36});
		b.add({0,a25-a15,-a13,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a25-a15,-a13,a45-a15,-a15,-a15-a56});
		b.add({0,a25-a15,-a13,a45-a15,-a15,a45-a46-a15});
		b.add({0,a25-a15,-a13,-a13-a36+a46,-a15,-a13-a36});
		b.add({0,a25-a15,-a13,a45-a15,-a15,-a13-a36});
		b.add({0,a25-a15,-a13,-a13-a34,-a15,-a15-a56});
		b.add({0,a25-a15,-a13,-a13-a34,-a15,-a13-a34-a46});
		b.add({0,a25-a15,-a13,-a13-a34,-a15,-a13-a36});
		b.add({0,a25-a15,-a13,a25-a26+a46-a15,-a15,a25-a26-a15});
		b.add({0,a25-a15,-a13,a45-a15,-a15,a25-a26-a15});
		b.add({0,a25-a15,-a13,-a13-a34,-a15,a25-a26-a15});
		b.add({0,a24+a46-a15-a56,-a13,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a24+a45-a15,-a13,a45-a15,-a15,-a15-a56});
		b.add({0,a24+a45-a15,-a13,a45-a15,-a15,a45-a46-a15});
		b.add({0,a24-a13-a36+a46,-a13,-a13-a36+a46,-a15,-a13-a36});
		b.add({0,a24+a45-a15,-a13,a45-a15,-a15,-a13-a36});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a15,-a15-a56});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a15,-a13-a34-a46});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a15,-a13-a36});
		b.add({0,a26-a15-a56,-a13,-a24+a26-a15-a56,-a15,-a15-a56});
		b.add({0,a24+a45-a15,-a13,a45-a15,-a15,a24-a26+a45-a15});
		b.add({0,a26-a13-a36,-a13,-a24+a26-a13-a36,-a15,-a13-a36});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a15,a24-a26-a13-a34});
		b.add({0,a25-a15,-a13,-a24+a25-a15,-a15,-a15-a56});
		b.add({0,a25-a15,-a13,-a24+a25-a15,-a15,-a24+a25-a46-a15});
		b.add({0,a25-a15,-a13,-a24+a25-a15,-a15,-a13-a36});
		b.add({0,a25-a15,-a13,-a24+a25-a15,-a15,a25-a26-a15});
		b.add({0,a23-a13,-a13,a46-a15-a56,-a15,-a15-a56});
		b.add({0,a23-a13,-a13,a45-a15,-a15,-a15-a56});
		b.add({0,a23-a13,-a13,a45-a15,-a15,a45-a46-a15});
		b.add({0,a23-a13,-a13,-a13-a36+a46,-a15,-a13-a36});
		b.add({0,a23-a13,-a13,a45-a15,-a15,-a13-a36});
		b.add({0,a23-a13,-a13,-a13-a34,-a15,-a15-a56});
		b.add({0,a23-a13,-a13,-a13-a34,-a15,-a13-a34-a46});
		b.add({0,a23-a13,-a13,-a13-a34,-a15,-a13-a36});
		b.add({0,a23-a13,-a13,a23-a26-a13+a46,-a15,a23-a26-a13});
		b.add({0,a23-a13,-a13,a45-a15,-a15,a23-a26-a13});
		b.add({0,a23-a13,-a13,-a13-a34,-a15,a23-a26-a13});
		b.add({0,a23-a13,-a13,a23-a24-a13,-a15,-a15-a56});
		b.add({0,a23-a13,-a13,a23-a24-a13,-a15,a23-a24-a13-a46});
		b.add({0,a23-a13,-a13,a23-a24-a13,-a15,-a13-a36});
		b.add({0,a23-a13,-a13,a23-a24-a13,-a15,a23-a26-a13});
		b.add({0,a26-a16,-a13,a46-a16,-a15,-a16});
		b.add({0,a26-a16,-a13,a45-a15,-a15,-a16});
		b.add({0,a26-a16,-a13,-a13-a34,-a15,-a16});
		b.add({0,a25-a15,-a13,a46-a16,-a15,-a16});
		b.add({0,a25-a15,-a13,a45-a15,-a15,-a16});
		b.add({0,a25-a15,-a13,-a13-a34,-a15,-a16});
		b.add({0,a24+a46-a16,-a13,a46-a16,-a15,-a16});
		b.add({0,a24+a45-a15,-a13,a45-a15,-a15,-a16});
		b.add({0,a24-a13-a34,-a13,-a13-a34,-a15,-a16});
		b.add({0,a26-a16,-a13,-a24+a26-a16,-a15,-a16});
		b.add({0,a25-a15,-a13,-a24+a25-a15,-a15,-a16});
		b.add({0,a23-a13,-a13,a46-a16,-a15,-a16});
		b.add({0,a23-a13,-a13,a45-a15,-a15,-a16});
		b.add({0,a23-a13,-a13,-a13-a34,-a15,-a16});
		b.add({0,a23-a13,-a13,a23-a24-a13,-a15,-a16});
		b.add({0,a26-a14-a46,-a13,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a26-a14-a45-a56,-a13,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a26-a14-a46,-a13,-a14,-a14-a45,-a14-a46});
		b.add({0,a26-a13-a36,-a13,-a14,-a13-a36+a56,-a13-a36});
		b.add({0,a26-a13-a36,-a13,-a14,-a14-a45,-a13-a36});
		b.add({0,a26-a13-a35-a56,-a13,-a14,-a13-a35,-a13-a35-a56});
		b.add({0,a26-a14-a46,-a13,-a14,-a13-a35,-a14-a46});
		b.add({0,a26-a13-a36,-a13,-a14,-a13-a35,-a13-a36});
		b.add({0,a25-a14-a46+a56,-a13,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a25-a14-a45,-a13,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a25-a14-a45,-a13,-a14,-a14-a45,-a14-a46});
		b.add({0,a25-a13-a36+a56,-a13,-a14,-a13-a36+a56,-a13-a36});
		b.add({0,a25-a14-a45,-a13,-a14,-a14-a45,-a13-a36});
		b.add({0,a25-a13-a35,-a13,-a14,-a13-a35,-a13-a35-a56});
		b.add({0,a25-a13-a35,-a13,-a14,-a13-a35,-a14-a46});
		b.add({0,a25-a13-a35,-a13,-a14,-a13-a35,-a13-a36});
		b.add({0,a26-a14-a46,-a13,-a14,-a25+a26-a14-a46,-a14-a46});
		b.add({0,a25-a14-a45,-a13,-a14,-a14-a45,a25-a26-a14-a45});
		b.add({0,a26-a13-a36,-a13,-a14,-a25+a26-a13-a36,-a13-a36});
		b.add({0,a25-a13-a35,-a13,-a14,-a13-a35,a25-a26-a13-a35});
		b.add({0,a24-a14,-a13,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a24-a14,-a13,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a24-a14,-a13,-a14,-a14-a45,-a14-a46});
		b.add({0,a24-a14,-a13,-a14,-a13-a36+a56,-a13-a36});
		b.add({0,a24-a14,-a13,-a14,-a14-a45,-a13-a36});
		b.add({0,a24-a14,-a13,-a14,-a13-a35,-a13-a35-a56});
		b.add({0,a24-a14,-a13,-a14,-a13-a35,-a14-a46});
		b.add({0,a24-a14,-a13,-a14,-a13-a35,-a13-a36});
		b.add({0,a24-a14,-a13,-a14,a24-a26-a14+a56,a24-a26-a14});
		b.add({0,a24-a14,-a13,-a14,-a14-a45,a24-a26-a14});
		b.add({0,a24-a14,-a13,-a14,-a13-a35,a24-a26-a14});
		b.add({0,a24-a14,-a13,-a14,a24-a25-a14,a24-a25-a14-a56});
		b.add({0,a24-a14,-a13,-a14,a24-a25-a14,-a14-a46});
		b.add({0,a24-a14,-a13,-a14,a24-a25-a14,-a13-a36});
		b.add({0,a24-a14,-a13,-a14,a24-a25-a14,a24-a26-a14});
		b.add({0,a23-a13,-a13,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,a23-a13,-a13,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,a23-a13,-a13,-a14,-a14-a45,-a14-a46});
		b.add({0,a23-a13,-a13,-a14,-a13-a36+a56,-a13-a36});
		b.add({0,a23-a13,-a13,-a14,-a14-a45,-a13-a36});
		b.add({0,a23-a13,-a13,-a14,-a13-a35,-a13-a35-a56});
		b.add({0,a23-a13,-a13,-a14,-a13-a35,-a14-a46});
		b.add({0,a23-a13,-a13,-a14,-a13-a35,-a13-a36});
		b.add({0,a23-a13,-a13,-a14,a23-a26-a13+a56,a23-a26-a13});
		b.add({0,a23-a13,-a13,-a14,-a14-a45,a23-a26-a13});
		b.add({0,a23-a13,-a13,-a14,-a13-a35,a23-a26-a13});
		b.add({0,a23-a13,-a13,-a14,a23-a25-a13,a23-a25-a13-a56});
		b.add({0,a23-a13,-a13,-a14,a23-a25-a13,-a14-a46});
		b.add({0,a23-a13,-a13,-a14,a23-a25-a13,-a13-a36});
		b.add({0,a23-a13,-a13,-a14,a23-a25-a13,a23-a26-a13});
		b.add({0,a26-a16,-a13,-a14,a56-a16,-a16});
		b.add({0,a26-a16,-a13,-a14,-a14-a45,-a16});
		b.add({0,a26-a16,-a13,-a14,-a13-a35,-a16});
		b.add({0,a25+a56-a16,-a13,-a14,a56-a16,-a16});
		b.add({0,a25-a14-a45,-a13,-a14,-a14-a45,-a16});
		b.add({0,a25-a13-a35,-a13,-a14,-a13-a35,-a16});
		b.add({0,a26-a16,-a13,-a14,-a25+a26-a16,-a16});
		b.add({0,a24-a14,-a13,-a14,a56-a16,-a16});
		b.add({0,a24-a14,-a13,-a14,-a14-a45,-a16});
		b.add({0,a24-a14,-a13,-a14,-a13-a35,-a16});
		b.add({0,a24-a14,-a13,-a14,a24-a25-a14,-a16});
		b.add({0,a23-a13,-a13,-a14,a56-a16,-a16});
		b.add({0,a23-a13,-a13,-a14,-a14-a45,-a16});
		b.add({0,a23-a13,-a13,-a14,-a13-a35,-a16});
		b.add({0,a23-a13,-a13,-a14,a23-a25-a13,-a16});
		b.add({0,a26-a15-a56,-a13,-a14,-a15,-a15-a56});
		b.add({0,a26-a14-a46,-a13,-a14,-a15,-a14-a46});
		b.add({0,a26-a13-a36,-a13,-a14,-a15,-a13-a36});
		b.add({0,a25-a15,-a13,-a14,-a15,-a15-a56});
		b.add({0,a25-a15,-a13,-a14,-a15,-a14-a46});
		b.add({0,a25-a15,-a13,-a14,-a15,-a13-a36});
		b.add({0,a25-a15,-a13,-a14,-a15,a25-a26-a15});
		b.add({0,a24-a14,-a13,-a14,-a15,-a15-a56});
		b.add({0,a24-a14,-a13,-a14,-a15,-a14-a46});
		b.add({0,a24-a14,-a13,-a14,-a15,-a13-a36});
		b.add({0,a24-a14,-a13,-a14,-a15,a24-a26-a14});
		b.add({0,a23-a13,-a13,-a14,-a15,-a15-a56});
		b.add({0,a23-a13,-a13,-a14,-a15,-a14-a46});
		b.add({0,a23-a13,-a13,-a14,-a15,-a13-a36});
		b.add({0,a23-a13,-a13,-a14,-a15,a23-a26-a13});
		b.add({0,a26-a16,-a13,-a14,-a15,-a16});
		b.add({0,a25-a15,-a13,-a14,-a15,-a16});
		b.add({0,a24-a14,-a13,-a14,-a15,-a16});
		b.add({0,a23-a13,-a13,-a14,-a15,-a16});
		b.add({0,-a12,-a12-a26+a36,-a12-a26+a46,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,-a12-a26+a45+a56,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,-a12-a26+a46,-a12-a26-a45+a46,-a12-a26});
		b.add({0,-a12,-a12-a26+a35+a56,-a12-a26+a46,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a26+a35+a56,-a12-a26+a45+a56,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a26+a35-a45+a46,-a12-a26+a46,-a12-a26-a45+a46,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,-a12-a26+a46,-a12-a26-a35+a36,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,-a12-a26-a35+a36+a45,-a12-a26-a35+a36,-a12-a26});
		b.add({0,-a12,-a12-a26+a34+a46,-a12-a26+a46,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a26+a34+a45+a56,-a12-a26+a45+a56,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a26+a34+a46,-a12-a26+a46,-a12-a26-a45+a46,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,-a12-a26-a34+a36,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,-a12-a26-a34+a36,-a12-a26-a34+a36-a45,-a12-a26});
		b.add({0,-a12,-a12-a26+a35+a56,-a12-a26-a34+a35+a56,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a26+a34+a46,-a12-a26+a46,-a12-a26+a34-a35+a46,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,-a12-a26-a34+a36,-a12-a26-a35+a36,-a12-a26});
		b.add({0,-a12,-a12-a25+a36-a56,-a12-a25+a46-a56,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a25+a36-a56,-a12-a25+a45,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a25+a36+a45-a46,-a12-a25+a45,-a12-a25,-a12-a25+a45-a46});
		b.add({0,-a12,-a12-a25+a35,-a12-a25+a46-a56,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a25+a35,-a12-a25+a45,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a25+a35,-a12-a25+a45,-a12-a25,-a12-a25+a45-a46});
		b.add({0,-a12,-a12-a25+a35,-a12-a25+a35-a36+a46,-a12-a25,-a12-a25+a35-a36});
		b.add({0,-a12,-a12-a25+a35,-a12-a25+a45,-a12-a25,-a12-a25+a35-a36});
		b.add({0,-a12,-a12-a25+a34+a46-a56,-a12-a25+a46-a56,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a25+a34+a45,-a12-a25+a45,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a25+a34+a45,-a12-a25+a45,-a12-a25,-a12-a25+a45-a46});
		b.add({0,-a12,-a12-a25+a36-a56,-a12-a25-a34+a36-a56,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a25+a34+a45,-a12-a25+a45,-a12-a25,-a12-a25+a34-a36+a45});
		b.add({0,-a12,-a12-a25+a35,-a12-a25-a34+a35,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a25+a35,-a12-a25-a34+a35,-a12-a25,-a12-a25-a34+a35-a46});
		b.add({0,-a12,-a12-a25+a35,-a12-a25-a34+a35,-a12-a25,-a12-a25+a35-a36});
		b.add({0,-a12,-a12-a26+a36,-a12-a26+a46,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,-a12-a25+a45,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a25+a35,-a12-a26+a46,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a25+a35,-a12-a25+a45,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a26+a34+a46,-a12-a26+a46,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a25+a34+a45,-a12-a25+a45,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,-a12-a26-a34+a36,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a25+a35,-a12-a25-a34+a35,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a24+a36-a46,-a12-a24,-a12-a24-a46+a56,-a12-a24-a46});
		b.add({0,-a12,-a12-a24+a36-a45-a56,-a12-a24,-a12-a24-a45,-a12-a24-a45-a56});
		b.add({0,-a12,-a12-a24+a36-a46,-a12-a24,-a12-a24-a45,-a12-a24-a46});
		b.add({0,-a12,-a12-a24+a35-a46+a56,-a12-a24,-a12-a24-a46+a56,-a12-a24-a46});
		b.add({0,-a12,-a12-a24+a35-a45,-a12-a24,-a12-a24-a45,-a12-a24-a45-a56});
		b.add({0,-a12,-a12-a24+a35-a45,-a12-a24,-a12-a24-a45,-a12-a24-a46});
		b.add({0,-a12,-a12-a24+a36-a46,-a12-a24,-a12-a24-a35+a36-a46,-a12-a24-a46});
		b.add({0,-a12,-a12-a24+a35-a45,-a12-a24,-a12-a24-a45,-a12-a24+a35-a36-a45});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a24-a46+a56,-a12-a24-a46});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a24-a45,-a12-a24-a45-a56});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a24-a45,-a12-a24-a46});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a24+a34-a36+a56,-a12-a24+a34-a36});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a24-a45,-a12-a24+a34-a36});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a24+a34-a35,-a12-a24+a34-a35-a56});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a24+a34-a35,-a12-a24-a46});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a24+a34-a35,-a12-a24+a34-a36});
		b.add({0,-a12,-a12-a26+a36,-a12-a24,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,-a12-a24,-a12-a24-a45,-a12-a26});
		b.add({0,-a12,-a12-a26+a35+a56,-a12-a24,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a24+a35-a45,-a12-a24,-a12-a24-a45,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,-a12-a24,-a12-a26-a35+a36,-a12-a26});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a24-a45,-a12-a26});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a24+a34-a35,-a12-a26});
		b.add({0,-a12,-a12-a25+a36-a56,-a12-a24,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a24+a36-a46,-a12-a24,-a12-a25,-a12-a24-a46});
		b.add({0,-a12,-a12-a25+a35,-a12-a24,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a25+a35,-a12-a24,-a12-a25,-a12-a24-a46});
		b.add({0,-a12,-a12-a25+a35,-a12-a24,-a12-a25,-a12-a25+a35-a36});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a25,-a12-a24-a46});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a25,-a12-a24+a34-a36});
		b.add({0,-a12,-a12-a26+a36,-a12-a24,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a25+a35,-a12-a24,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a23-a36+a46,-a12-a23-a36+a56,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a23-a36+a45+a56,-a12-a23-a36+a56,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a23-a36+a46,-a12-a23-a36-a45+a46,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a23-a35+a46-a56,-a12-a23-a35,-a12-a23-a35-a56});
		b.add({0,-a12,-a12-a23,-a12-a23-a35+a45,-a12-a23-a35,-a12-a23-a35-a56});
		b.add({0,-a12,-a12-a23,-a12-a23-a35+a45,-a12-a23-a35,-a12-a23-a35+a45-a46});
		b.add({0,-a12,-a12-a23,-a12-a23-a36+a46,-a12-a23-a35,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a23-a35+a45,-a12-a23-a35,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a23-a34-a46+a56,-a12-a23-a34-a46});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a23-a34-a45,-a12-a23-a34-a45-a56});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a23-a34-a45,-a12-a23-a34-a46});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a23-a36+a56,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a23-a34-a45,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a23-a35,-a12-a23-a35-a56});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a23-a35,-a12-a23-a34-a46});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a23-a35,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a26+a46,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a26+a45+a56,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a26+a46,-a12-a26-a45+a46,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a26+a46,-a12-a23-a35,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a23-a35+a45,-a12-a23-a35,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a23-a34-a45,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a23-a35,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a25+a46-a56,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a23,-a12-a25+a45,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a23,-a12-a25+a45,-a12-a25,-a12-a25+a45-a46});
		b.add({0,-a12,-a12-a23,-a12-a23-a36+a46,-a12-a25,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a25+a45,-a12-a25,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a25,-a12-a23-a34-a46});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a25,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a26+a46,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a25+a45,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a24-a46+a56,-a12-a24-a46});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a24-a45,-a12-a24-a45-a56});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a24-a45,-a12-a24-a46});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a23-a36+a56,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a24-a45,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a23-a35,-a12-a23-a35-a56});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a23-a35,-a12-a24-a46});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a23-a35,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a24-a45,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a23-a35,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a25,-a12-a24-a46});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a25,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a25,-a12-a26});
		b.add({0,-a12,a36-a16,a46-a16,a56-a16,-a16});
		b.add({0,-a12,a36-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,-a12,a36-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,-a12,a35+a56-a16,a46-a16,a56-a16,-a16});
		b.add({0,-a12,a35+a56-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,-a12,a35-a45+a46-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,-a12,a36-a16,a46-a16,-a35+a36-a16,-a16});
		b.add({0,-a12,a36-a16,-a35+a36+a45-a16,-a35+a36-a16,-a16});
		b.add({0,-a12,a34+a46-a16,a46-a16,a56-a16,-a16});
		b.add({0,-a12,a34+a45+a56-a16,a45+a56-a16,a56-a16,-a16});
		b.add({0,-a12,a34+a46-a16,a46-a16,-a45+a46-a16,-a16});
		b.add({0,-a12,a36-a16,-a34+a36-a16,a56-a16,-a16});
		b.add({0,-a12,a36-a16,-a34+a36-a16,-a34+a36-a45-a16,-a16});
		b.add({0,-a12,a35+a56-a16,-a34+a35+a56-a16,a56-a16,-a16});
		b.add({0,-a12,a34+a46-a16,a46-a16,a34-a35+a46-a16,-a16});
		b.add({0,-a12,a36-a16,-a34+a36-a16,-a35+a36-a16,-a16});
		b.add({0,-a12,a36-a16,a46-a16,-a12-a25,-a16});
		b.add({0,-a12,a36-a16,-a12-a25+a45,-a12-a25,-a16});
		b.add({0,-a12,-a12-a25+a35,a46-a16,-a12-a25,-a16});
		b.add({0,-a12,-a12-a25+a35,-a12-a25+a45,-a12-a25,-a16});
		b.add({0,-a12,a34+a46-a16,a46-a16,-a12-a25,-a16});
		b.add({0,-a12,-a12-a25+a34+a45,-a12-a25+a45,-a12-a25,-a16});
		b.add({0,-a12,a36-a16,-a34+a36-a16,-a12-a25,-a16});
		b.add({0,-a12,-a12-a25+a35,-a12-a25-a34+a35,-a12-a25,-a16});
		b.add({0,-a12,a36-a16,-a12-a24,a56-a16,-a16});
		b.add({0,-a12,a36-a16,-a12-a24,-a12-a24-a45,-a16});
		b.add({0,-a12,a35+a56-a16,-a12-a24,a56-a16,-a16});
		b.add({0,-a12,-a12-a24+a35-a45,-a12-a24,-a12-a24-a45,-a16});
		b.add({0,-a12,a36-a16,-a12-a24,-a35+a36-a16,-a16});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,a56-a16,-a16});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a24-a45,-a16});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a24+a34-a35,-a16});
		b.add({0,-a12,a36-a16,-a12-a24,-a12-a25,-a16});
		b.add({0,-a12,-a12-a25+a35,-a12-a24,-a12-a25,-a16});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a12-a25,-a16});
		b.add({0,-a12,-a12-a23,a46-a16,a56-a16,-a16});
		b.add({0,-a12,-a12-a23,a45+a56-a16,a56-a16,-a16});
		b.add({0,-a12,-a12-a23,a46-a16,-a45+a46-a16,-a16});
		b.add({0,-a12,-a12-a23,a46-a16,-a12-a23-a35,-a16});
		b.add({0,-a12,-a12-a23,-a12-a23-a35+a45,-a12-a23-a35,-a16});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,a56-a16,-a16});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a23-a34-a45,-a16});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a23-a35,-a16});
		b.add({0,-a12,-a12-a23,a46-a16,-a12-a25,-a16});
		b.add({0,-a12,-a12-a23,-a12-a25+a45,-a12-a25,-a16});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a12-a25,-a16});
		b.add({0,-a12,-a12-a23,-a12-a24,a56-a16,-a16});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a24-a45,-a16});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a23-a35,-a16});
		b.add({0,-a12,-a12-a23,-a12-a24,-a12-a25,-a16});
		b.add({0,-a12,a36-a15-a56,a46-a15-a56,-a15,-a15-a56});
		b.add({0,-a12,a36-a15-a56,a45-a15,-a15,-a15-a56});
		b.add({0,-a12,a36+a45-a46-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,-a12,a35-a15,a46-a15-a56,-a15,-a15-a56});
		b.add({0,-a12,a35-a15,a45-a15,-a15,-a15-a56});
		b.add({0,-a12,a35-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,-a12,a35-a15,a35-a36+a46-a15,-a15,a35-a36-a15});
		b.add({0,-a12,a35-a15,a45-a15,-a15,a35-a36-a15});
		b.add({0,-a12,a34+a46-a15-a56,a46-a15-a56,-a15,-a15-a56});
		b.add({0,-a12,a34+a45-a15,a45-a15,-a15,-a15-a56});
		b.add({0,-a12,a34+a45-a15,a45-a15,-a15,a45-a46-a15});
		b.add({0,-a12,a36-a15-a56,-a34+a36-a15-a56,-a15,-a15-a56});
		b.add({0,-a12,a34+a45-a15,a45-a15,-a15,a34-a36+a45-a15});
		b.add({0,-a12,a35-a15,-a34+a35-a15,-a15,-a15-a56});
		b.add({0,-a12,a35-a15,-a34+a35-a15,-a15,-a34+a35-a46-a15});
		b.add({0,-a12,a35-a15,-a34+a35-a15,-a15,a35-a36-a15});
		b.add({0,-a12,-a12-a26+a36,-a12-a26+a46,-a15,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,a45-a15,-a15,-a12-a26});
		b.add({0,-a12,a35-a15,-a12-a26+a46,-a15,-a12-a26});
		b.add({0,-a12,a35-a15,a45-a15,-a15,-a12-a26});
		b.add({0,-a12,-a12-a26+a34+a46,-a12-a26+a46,-a15,-a12-a26});
		b.add({0,-a12,a34+a45-a15,a45-a15,-a15,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,-a12-a26-a34+a36,-a15,-a12-a26});
		b.add({0,-a12,a35-a15,-a34+a35-a15,-a15,-a12-a26});
		b.add({0,-a12,a36-a15-a56,-a12-a24,-a15,-a15-a56});
		b.add({0,-a12,-a12-a24+a36-a46,-a12-a24,-a15,-a12-a24-a46});
		b.add({0,-a12,a35-a15,-a12-a24,-a15,-a15-a56});
		b.add({0,-a12,a35-a15,-a12-a24,-a15,-a12-a24-a46});
		b.add({0,-a12,a35-a15,-a12-a24,-a15,a35-a36-a15});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a15,-a15-a56});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a15,-a12-a24-a46});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a15,-a12-a24+a34-a36});
		b.add({0,-a12,-a12-a26+a36,-a12-a24,-a15,-a12-a26});
		b.add({0,-a12,a35-a15,-a12-a24,-a15,-a12-a26});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a15,-a12-a26});
		b.add({0,-a12,-a12-a23,a46-a15-a56,-a15,-a15-a56});
		b.add({0,-a12,-a12-a23,a45-a15,-a15,-a15-a56});
		b.add({0,-a12,-a12-a23,a45-a15,-a15,a45-a46-a15});
		b.add({0,-a12,-a12-a23,-a12-a23-a36+a46,-a15,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,a45-a15,-a15,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a15,-a15-a56});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a15,-a12-a23-a34-a46});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a15,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a26+a46,-a15,-a12-a26});
		b.add({0,-a12,-a12-a23,a45-a15,-a15,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a15,-a12-a26});
		b.add({0,-a12,-a12-a23,-a12-a24,-a15,-a15-a56});
		b.add({0,-a12,-a12-a23,-a12-a24,-a15,-a12-a24-a46});
		b.add({0,-a12,-a12-a23,-a12-a24,-a15,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a12-a24,-a15,-a12-a26});
		b.add({0,-a12,a36-a16,a46-a16,-a15,-a16});
		b.add({0,-a12,a36-a16,a45-a15,-a15,-a16});
		b.add({0,-a12,a35-a15,a46-a16,-a15,-a16});
		b.add({0,-a12,a35-a15,a45-a15,-a15,-a16});
		b.add({0,-a12,a34+a46-a16,a46-a16,-a15,-a16});
		b.add({0,-a12,a34+a45-a15,a45-a15,-a15,-a16});
		b.add({0,-a12,a36-a16,-a34+a36-a16,-a15,-a16});
		b.add({0,-a12,a35-a15,-a34+a35-a15,-a15,-a16});
		b.add({0,-a12,a36-a16,-a12-a24,-a15,-a16});
		b.add({0,-a12,a35-a15,-a12-a24,-a15,-a16});
		b.add({0,-a12,-a12-a24+a34,-a12-a24,-a15,-a16});
		b.add({0,-a12,-a12-a23,a46-a16,-a15,-a16});
		b.add({0,-a12,-a12-a23,a45-a15,-a15,-a16});
		b.add({0,-a12,-a12-a23,-a12-a23-a34,-a15,-a16});
		b.add({0,-a12,-a12-a23,-a12-a24,-a15,-a16});
		b.add({0,-a12,a36-a14-a46,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,-a12,a36-a14-a45-a56,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,-a12,a36-a14-a46,-a14,-a14-a45,-a14-a46});
		b.add({0,-a12,a35-a14-a46+a56,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,-a12,a35-a14-a45,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,-a12,a35-a14-a45,-a14,-a14-a45,-a14-a46});
		b.add({0,-a12,a36-a14-a46,-a14,-a35+a36-a14-a46,-a14-a46});
		b.add({0,-a12,a35-a14-a45,-a14,-a14-a45,a35-a36-a14-a45});
		b.add({0,-a12,a34-a14,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,-a12,a34-a14,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,-a12,a34-a14,-a14,-a14-a45,-a14-a46});
		b.add({0,-a12,a34-a14,-a14,a34-a36-a14+a56,a34-a36-a14});
		b.add({0,-a12,a34-a14,-a14,-a14-a45,a34-a36-a14});
		b.add({0,-a12,a34-a14,-a14,a34-a35-a14,a34-a35-a14-a56});
		b.add({0,-a12,a34-a14,-a14,a34-a35-a14,-a14-a46});
		b.add({0,-a12,a34-a14,-a14,a34-a35-a14,a34-a36-a14});
		b.add({0,-a12,-a12-a26+a36,-a14,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,-a14,-a14-a45,-a12-a26});
		b.add({0,-a12,-a12-a26+a35+a56,-a14,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,a35-a14-a45,-a14,-a14-a45,-a12-a26});
		b.add({0,-a12,-a12-a26+a36,-a14,-a12-a26-a35+a36,-a12-a26});
		b.add({0,-a12,a34-a14,-a14,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,a34-a14,-a14,-a14-a45,-a12-a26});
		b.add({0,-a12,a34-a14,-a14,a34-a35-a14,-a12-a26});
		b.add({0,-a12,-a12-a25+a36-a56,-a14,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,a36-a14-a46,-a14,-a12-a25,-a14-a46});
		b.add({0,-a12,-a12-a25+a35,-a14,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a25+a35,-a14,-a12-a25,-a14-a46});
		b.add({0,-a12,-a12-a25+a35,-a14,-a12-a25,-a12-a25+a35-a36});
		b.add({0,-a12,a34-a14,-a14,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,a34-a14,-a14,-a12-a25,-a14-a46});
		b.add({0,-a12,a34-a14,-a14,-a12-a25,a34-a36-a14});
		b.add({0,-a12,-a12-a26+a36,-a14,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a25+a35,-a14,-a12-a25,-a12-a26});
		b.add({0,-a12,a34-a14,-a14,-a12-a25,-a12-a26});
		b.add({0,-a12,-a12-a23,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,-a12,-a12-a23,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,-a12,-a12-a23,-a14,-a14-a45,-a14-a46});
		b.add({0,-a12,-a12-a23,-a14,-a12-a23-a36+a56,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a14,-a14-a45,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a14,-a12-a23-a35,-a12-a23-a35-a56});
		b.add({0,-a12,-a12-a23,-a14,-a12-a23-a35,-a14-a46});
		b.add({0,-a12,-a12-a23,-a14,-a12-a23-a35,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a14,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a12-a23,-a14,-a14-a45,-a12-a26});
		b.add({0,-a12,-a12-a23,-a14,-a12-a23-a35,-a12-a26});
		b.add({0,-a12,-a12-a23,-a14,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a12-a23,-a14,-a12-a25,-a14-a46});
		b.add({0,-a12,-a12-a23,-a14,-a12-a25,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a14,-a12-a25,-a12-a26});
		b.add({0,-a12,a36-a16,-a14,a56-a16,-a16});
		b.add({0,-a12,a36-a16,-a14,-a14-a45,-a16});
		b.add({0,-a12,a35+a56-a16,-a14,a56-a16,-a16});
		b.add({0,-a12,a35-a14-a45,-a14,-a14-a45,-a16});
		b.add({0,-a12,a36-a16,-a14,-a35+a36-a16,-a16});
		b.add({0,-a12,a34-a14,-a14,a56-a16,-a16});
		b.add({0,-a12,a34-a14,-a14,-a14-a45,-a16});
		b.add({0,-a12,a34-a14,-a14,a34-a35-a14,-a16});
		b.add({0,-a12,a36-a16,-a14,-a12-a25,-a16});
		b.add({0,-a12,-a12-a25+a35,-a14,-a12-a25,-a16});
		b.add({0,-a12,a34-a14,-a14,-a12-a25,-a16});
		b.add({0,-a12,-a12-a23,-a14,a56-a16,-a16});
		b.add({0,-a12,-a12-a23,-a14,-a14-a45,-a16});
		b.add({0,-a12,-a12-a23,-a14,-a12-a23-a35,-a16});
		b.add({0,-a12,-a12-a23,-a14,-a12-a25,-a16});
		b.add({0,-a12,a36-a15-a56,-a14,-a15,-a15-a56});
		b.add({0,-a12,a36-a14-a46,-a14,-a15,-a14-a46});
		b.add({0,-a12,a35-a15,-a14,-a15,-a15-a56});
		b.add({0,-a12,a35-a15,-a14,-a15,-a14-a46});
		b.add({0,-a12,a35-a15,-a14,-a15,a35-a36-a15});
		b.add({0,-a12,a34-a14,-a14,-a15,-a15-a56});
		b.add({0,-a12,a34-a14,-a14,-a15,-a14-a46});
		b.add({0,-a12,a34-a14,-a14,-a15,a34-a36-a14});
		b.add({0,-a12,-a12-a26+a36,-a14,-a15,-a12-a26});
		b.add({0,-a12,a35-a15,-a14,-a15,-a12-a26});
		b.add({0,-a12,a34-a14,-a14,-a15,-a12-a26});
		b.add({0,-a12,-a12-a23,-a14,-a15,-a15-a56});
		b.add({0,-a12,-a12-a23,-a14,-a15,-a14-a46});
		b.add({0,-a12,-a12-a23,-a14,-a15,-a12-a23-a36});
		b.add({0,-a12,-a12-a23,-a14,-a15,-a12-a26});
		b.add({0,-a12,a36-a16,-a14,-a15,-a16});
		b.add({0,-a12,a35-a15,-a14,-a15,-a16});
		b.add({0,-a12,a34-a14,-a14,-a15,-a16});
		b.add({0,-a12,-a12-a23,-a14,-a15,-a16});
		b.add({0,-a12,-a13,-a13-a36+a46,-a13-a36+a56,-a13-a36});
		b.add({0,-a12,-a13,-a13-a36+a45+a56,-a13-a36+a56,-a13-a36});
		b.add({0,-a12,-a13,-a13-a36+a46,-a13-a36-a45+a46,-a13-a36});
		b.add({0,-a12,-a13,-a13-a35+a46-a56,-a13-a35,-a13-a35-a56});
		b.add({0,-a12,-a13,-a13-a35+a45,-a13-a35,-a13-a35-a56});
		b.add({0,-a12,-a13,-a13-a35+a45,-a13-a35,-a13-a35+a45-a46});
		b.add({0,-a12,-a13,-a13-a36+a46,-a13-a35,-a13-a36});
		b.add({0,-a12,-a13,-a13-a35+a45,-a13-a35,-a13-a36});
		b.add({0,-a12,-a13,-a13-a34,-a13-a34-a46+a56,-a13-a34-a46});
		b.add({0,-a12,-a13,-a13-a34,-a13-a34-a45,-a13-a34-a45-a56});
		b.add({0,-a12,-a13,-a13-a34,-a13-a34-a45,-a13-a34-a46});
		b.add({0,-a12,-a13,-a13-a34,-a13-a36+a56,-a13-a36});
		b.add({0,-a12,-a13,-a13-a34,-a13-a34-a45,-a13-a36});
		b.add({0,-a12,-a13,-a13-a34,-a13-a35,-a13-a35-a56});
		b.add({0,-a12,-a13,-a13-a34,-a13-a35,-a13-a34-a46});
		b.add({0,-a12,-a13,-a13-a34,-a13-a35,-a13-a36});
		b.add({0,-a12,-a13,-a12-a26+a46,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a13,-a12-a26+a45+a56,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a13,-a12-a26+a46,-a12-a26-a45+a46,-a12-a26});
		b.add({0,-a12,-a13,-a12-a26+a46,-a13-a35,-a12-a26});
		b.add({0,-a12,-a13,-a13-a35+a45,-a13-a35,-a12-a26});
		b.add({0,-a12,-a13,-a13-a34,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a13,-a13-a34,-a13-a34-a45,-a12-a26});
		b.add({0,-a12,-a13,-a13-a34,-a13-a35,-a12-a26});
		b.add({0,-a12,-a13,-a12-a25+a46-a56,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a13,-a12-a25+a45,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a13,-a12-a25+a45,-a12-a25,-a12-a25+a45-a46});
		b.add({0,-a12,-a13,-a13-a36+a46,-a12-a25,-a13-a36});
		b.add({0,-a12,-a13,-a12-a25+a45,-a12-a25,-a13-a36});
		b.add({0,-a12,-a13,-a13-a34,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a13,-a13-a34,-a12-a25,-a13-a34-a46});
		b.add({0,-a12,-a13,-a13-a34,-a12-a25,-a13-a36});
		b.add({0,-a12,-a13,-a12-a26+a46,-a12-a25,-a12-a26});
		b.add({0,-a12,-a13,-a12-a25+a45,-a12-a25,-a12-a26});
		b.add({0,-a12,-a13,-a13-a34,-a12-a25,-a12-a26});
		b.add({0,-a12,-a13,-a12-a24,-a12-a24-a46+a56,-a12-a24-a46});
		b.add({0,-a12,-a13,-a12-a24,-a12-a24-a45,-a12-a24-a45-a56});
		b.add({0,-a12,-a13,-a12-a24,-a12-a24-a45,-a12-a24-a46});
		b.add({0,-a12,-a13,-a12-a24,-a13-a36+a56,-a13-a36});
		b.add({0,-a12,-a13,-a12-a24,-a12-a24-a45,-a13-a36});
		b.add({0,-a12,-a13,-a12-a24,-a13-a35,-a13-a35-a56});
		b.add({0,-a12,-a13,-a12-a24,-a13-a35,-a12-a24-a46});
		b.add({0,-a12,-a13,-a12-a24,-a13-a35,-a13-a36});
		b.add({0,-a12,-a13,-a12-a24,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a13,-a12-a24,-a12-a24-a45,-a12-a26});
		b.add({0,-a12,-a13,-a12-a24,-a13-a35,-a12-a26});
		b.add({0,-a12,-a13,-a12-a24,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a13,-a12-a24,-a12-a25,-a12-a24-a46});
		b.add({0,-a12,-a13,-a12-a24,-a12-a25,-a13-a36});
		b.add({0,-a12,-a13,-a12-a24,-a12-a25,-a12-a26});
		b.add({0,-a12,-a13,a46-a16,a56-a16,-a16});
		b.add({0,-a12,-a13,a45+a56-a16,a56-a16,-a16});
		b.add({0,-a12,-a13,a46-a16,-a45+a46-a16,-a16});
		b.add({0,-a12,-a13,a46-a16,-a13-a35,-a16});
		b.add({0,-a12,-a13,-a13-a35+a45,-a13-a35,-a16});
		b.add({0,-a12,-a13,-a13-a34,a56-a16,-a16});
		b.add({0,-a12,-a13,-a13-a34,-a13-a34-a45,-a16});
		b.add({0,-a12,-a13,-a13-a34,-a13-a35,-a16});
		b.add({0,-a12,-a13,a46-a16,-a12-a25,-a16});
		b.add({0,-a12,-a13,-a12-a25+a45,-a12-a25,-a16});
		b.add({0,-a12,-a13,-a13-a34,-a12-a25,-a16});
		b.add({0,-a12,-a13,-a12-a24,a56-a16,-a16});
		b.add({0,-a12,-a13,-a12-a24,-a12-a24-a45,-a16});
		b.add({0,-a12,-a13,-a12-a24,-a13-a35,-a16});
		b.add({0,-a12,-a13,-a12-a24,-a12-a25,-a16});
		b.add({0,-a12,-a13,a46-a15-a56,-a15,-a15-a56});
		b.add({0,-a12,-a13,a45-a15,-a15,-a15-a56});
		b.add({0,-a12,-a13,a45-a15,-a15,a45-a46-a15});
		b.add({0,-a12,-a13,-a13-a36+a46,-a15,-a13-a36});
		b.add({0,-a12,-a13,a45-a15,-a15,-a13-a36});
		b.add({0,-a12,-a13,-a13-a34,-a15,-a15-a56});
		b.add({0,-a12,-a13,-a13-a34,-a15,-a13-a34-a46});
		b.add({0,-a12,-a13,-a13-a34,-a15,-a13-a36});
		b.add({0,-a12,-a13,-a12-a26+a46,-a15,-a12-a26});
		b.add({0,-a12,-a13,a45-a15,-a15,-a12-a26});
		b.add({0,-a12,-a13,-a13-a34,-a15,-a12-a26});
		b.add({0,-a12,-a13,-a12-a24,-a15,-a15-a56});
		b.add({0,-a12,-a13,-a12-a24,-a15,-a12-a24-a46});
		b.add({0,-a12,-a13,-a12-a24,-a15,-a13-a36});
		b.add({0,-a12,-a13,-a12-a24,-a15,-a12-a26});
		b.add({0,-a12,-a13,a46-a16,-a15,-a16});
		b.add({0,-a12,-a13,a45-a15,-a15,-a16});
		b.add({0,-a12,-a13,-a13-a34,-a15,-a16});
		b.add({0,-a12,-a13,-a12-a24,-a15,-a16});
		b.add({0,-a12,-a13,-a14,-a14-a46+a56,-a14-a46});
		b.add({0,-a12,-a13,-a14,-a14-a45,-a14-a45-a56});
		b.add({0,-a12,-a13,-a14,-a14-a45,-a14-a46});
		b.add({0,-a12,-a13,-a14,-a13-a36+a56,-a13-a36});
		b.add({0,-a12,-a13,-a14,-a14-a45,-a13-a36});
		b.add({0,-a12,-a13,-a14,-a13-a35,-a13-a35-a56});
		b.add({0,-a12,-a13,-a14,-a13-a35,-a14-a46});
		b.add({0,-a12,-a13,-a14,-a13-a35,-a13-a36});
		b.add({0,-a12,-a13,-a14,-a12-a26+a56,-a12-a26});
		b.add({0,-a12,-a13,-a14,-a14-a45,-a12-a26});
		b.add({0,-a12,-a13,-a14,-a13-a35,-a12-a26});
		b.add({0,-a12,-a13,-a14,-a12-a25,-a12-a25-a56});
		b.add({0,-a12,-a13,-a14,-a12-a25,-a14-a46});
		b.add({0,-a12,-a13,-a14,-a12-a25,-a13-a36});
		b.add({0,-a12,-a13,-a14,-a12-a25,-a12-a26});
		b.add({0,-a12,-a13,-a14,a56-a16,-a16});
		b.add({0,-a12,-a13,-a14,-a14-a45,-a16});
		b.add({0,-a12,-a13,-a14,-a13-a35,-a16});
		b.add({0,-a12,-a13,-a14,-a12-a25,-a16});
		b.add({0,-a12,-a13,-a14,-a15,-a15-a56});
		b.add({0,-a12,-a13,-a14,-a15,-a14-a46});
		b.add({0,-a12,-a13,-a14,-a15,-a13-a36});
		b.add({0,-a12,-a13,-a14,-a15,-a12-a26});
		b.add({0,-a12,-a13,-a14,-a15,-a16});

		b.divideBy(1296.0);

		return b.getData();
	} else {
		throw "Incorrect dimension while calculating means of spans!";
	}
}

template<size_t N>
bool Matrix<N>::testAvgSpanTreeParetoOptimal() const {
	return Matrix::testVectorParetoOptimal(getMeanOfSpans());
}

template<size_t N>
Matrix<N - 1> Matrix<N>::cutBottom() const {
	std::vector<Ush> v(data.begin() + N - 1, data.end());
	return Matrix<N - 1>(v);
}

/*template<size_t N>
Ush Matrix<N>::countParetoVectorsByAlgorithm() const {
	MatrixProcessor<N> mp(*this);
	return mp.run();
}*/

template<size_t N>
void Matrix<N>::regularize() {
	//first lets order by the number of >=1 element from the upper triangle
	std::vector<Matrix<N>> temporaryMatrices;
	Ush maxNumOfIntegerElements = 0;

	Ush perm[N];
	for (Ush j = 0; j < N; j++) perm[j] = j;

	while (std::next_permutation(perm, perm + N)) {
		Matrix<N> tmp = permutateBy(perm);
		if (tmp.countIndexOfUpperTriangle () > maxNumOfIntegerElements) {
			temporaryMatrices.clear ();
			temporaryMatrices.push_back(tmp);
			maxNumOfIntegerElements = countIndexOfUpperTriangle ();
		}
	}

	this->data = std::min_element( std::begin(temporaryMatrices), std::end(temporaryMatrices) )->data;

	return;
}

template<size_t N>
Ush Matrix<N>::countIndexOfUpperTriangle () const {
	Ush k = 0;
	for (Ush i = 0; i < N * (N - 1) / 2; i++) {
		if (data[i] < 9) {
			k++;
		}
	}
	return k;
}
