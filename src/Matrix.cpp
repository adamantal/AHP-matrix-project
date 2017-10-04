#include "Matrix.h"
#include "LpSolution.h"

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
	if (i == 0) {
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
bool Matrix::operator==(const Matrix &rhs) const{
	return data == rhs.data;
}
bool Matrix::operator<(const Matrix &rhs) const{
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i] < rhs.data[i]) return true;
	}
	return false;
}
double Matrix::get(const int &i, const int &j)const{
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

Eigen::MatrixXd Matrix::toEigenMatrix()const{
	Eigen::MatrixXd m(4,4);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			m(i,j) = elem[indexOfElement(i,j)];
		}
	}
	return m;
}

double Matrix::largestEigenvalue()const{
	Eigen::MatrixXd m(4,4);
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

double Matrix::getConsistencyRatio()const{
	Eigen::MatrixXd m(4,4);
	m = toEigenMatrix();

	return (largestEigenvalue()-m.rows())/(ConsistencyIndex[m.rows()]-m.rows());
}


std::vector<double> Matrix::getPrimalEigenvector()const{
	Eigen::MatrixXd m(4,4);
	m = toEigenMatrix();

	Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(m);

	double largestEigenvalue = eigenSolver.eigenvalues()[0].real();
	int which = 0;
	short multiplicity = 1;

	for (int i = 1; i < eigenSolver.eigenvalues().size(); i++) {
		if (largestEigenvalue < eigenSolver.eigenvalues()[i].real()) {
			if (abs(largestEigenvalue - eigenSolver.eigenvalues()[i].real()) < 1e-6) {
				multiplicity++;
			} else {
				largestEigenvalue = eigenSolver.eigenvalues()[i].real();
				which = i;
				multiplicity = 1;
			}
		}
	}

	std::vector<double> lambda_1;
	for (int i = 0; i < 4; i++) {
		double tmpelement = eigenSolver.eigenvectors().col(which)[i].real();
		lambda_1.push_back(tmpelement);
	}

	/*for (int i = 0; i < 4; i ++) std::cout << lambda_1[i] << " ";
	std::cout << std::endl;*/

	return lambda_1;
}

std::vector<double> Matrix::getPrimalNormEigenvector()const{
	std::vector<double> v = getPrimalEigenvector();
	std::vector<double> tmp;
	for (int i = 0; i < 4; i++) {
		tmp.push_back(v[i] / v[0]);
	}
	/*for (int i = 0; i < 4; i ++) std::cout << tmp[i] << " ";
	std::cout << std::endl; */
	return tmp;
}

void Matrix::L1(std::vector<double> &v){
	double sum = 0;
	for (int i = 0; i < 4; i++) sum += v[i];
	for (int i = 0; i < 4; i++) v[i] /= sum;
}

bool Matrix::testPrimalEigenvectorIsParetoOptimal()const{
	std::vector<double> lambda_1 = getPrimalNormEigenvector();

	return Matrix::testVectorParetoOptimal(lambda_1);
}

bool Matrix::testVectorParetoOptimal(const std::vector<double> &w) const {
	//std::cout << "Return of the LP is " << Matrix::testVectorParetoOptimal(lambda_1) << std::endl;
	if (Matrix::LPVectorParetoOptimal(w).isOptimal()) return true;
	return false;
}

bool Matrix::testParetoOptimality(filterType filter) const{
	switch (filter){
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

LpSolution Matrix::LPVectorParetoOptimal(const std::vector<double> &w) const {
	//Inputs: A matrix and w vector

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
		for (int i = 0; i < 4; i++) {
			v.push_back(log(w[i]));
			std::vector<double> logmrow;
			for (int j = 0; j < 4; j++) {
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

		for (int i = 0; i < 4; i++) {
			for (int j = i + 1; j < 4; j++) {
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
		for (int i = 0; i < 4; i++) { //y-k definialva
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

			LpSolution ret(*this, w, I, J, LP.primal(), xvector, svector);
			return ret;
		} else {
			if (verbosity) {
				std::cout << "Optimal solution not found." << std::endl;
			}
			std::cout << "There is an error during making the LP. \nERROR INFO: " << LP.primalType() << "\n";
			throw "There is an error during making the LP.";
		}
}
bool Matrix::testCosineParetoOptimal()const{
	const int n = 4;

	//initialising b:
	std::vector<std::vector<double> > b;

	for (size_t i = 0; i < n; i++) {
		std::vector<double> tmp;
		for (size_t j = 0; j < n; j++) {
			tmp.push_back(0.0);
		}
		b.push_back(tmp);
	}

	//calculate the b_ij-s:
	for (size_t col = 0; col < n; col++) {
		double colSum = 0.0;
		for (size_t row = 0; row < n; row++) {
			colSum += get(row, col) * get(row, col);
		}
		colSum = sqrt(colSum);
		for (size_t row = 0; row < n; row++) {
			b[row][col] = Matrix::get(row, col) / colSum;
		}
	}

	//calculate formulas:
	std::vector<double> w;

	//nominators:
	for (size_t row = 0; row < n; row++) {
		double tmpsum = 0.0;
		for (size_t col = 0; col < n; col++) {
			tmpsum += b[row][col];
		}
		w.push_back(tmpsum);
	}

	//denominator:
	double fullSum = 0.0;
	for (size_t i = 0; i < n; i++) {
		for (size_t row = 0; row < n; row++) {
			double colSum = 0.0;
			for (size_t col = 0; col < n; col++) {
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

	return Matrix::testVectorParetoOptimal(w);
}

std::vector<double> Matrix::getMeanOfSpans()const{
	struct VectorAdd{
		std::vector<double> n;

		VectorAdd(){
			std::vector<double> k;
			for (int i = 0; i < 4; i++){
				k.push_back(0.0);
			}
			n = k;
		}
		VectorAdd(std::vector<double> v):n(v){}
		void add(std::vector<double> v) {
			double sum = 0;
			for (int i = 0; i < 4; i++) v[i] = exp(v[i]);
			for (int i = 0; i < 4; i++) sum += v[i];
			for (int i = 0; i < 4; i++) v[i] /= sum;

			for (int i = 0; i < 4; i++){
				n[i] += v[i];
			}
		}
		void divideBy(double x){
			for (int i = 0; i < 4; i++) n[i] /= x;
		}
		std::vector<double> getData(){
			return n;
		}
	};

	double a12,a13,a14,a23,a24,a34;
	a12 = log(get(0,1));
	a13 = log(get(0,2));
	a14 = log(get(0,3));
	a23 = log(get(1,2));
	a24 = log(get(1,3));
	a34 = log(get(2,3));
	std::vector<double> v;
	VectorAdd b,c;

	v = {0,a24-a14,a34-a14,-a14};
	b.add(v);
	v = {0,a23+a34-a14,a34-a14,-a14};
	b.add(v);
	v = {0,a24-a14,-a23+a24-a14,-a14};

	b.add(v);
	v = {0,a24-a13-a34,-a13,-a13-a34};

	b.add(v);
	v = {0,a23-a13,-a13,-a13-a34};

	b.add(v);
	v = {0,a23-a13,-a13,a23-a24-a13};

	b.add(v);
	v = {0,a24-a14,-a13,-a14};

	b.add(v);
	v = {0,a23-a13,-a13,-a14};

	b.add(v);
	v = {0,-a12,-a12-a24+a34,-a12-a24};

	b.add(v);
	v = {0,-a12,-a12-a23,-a12-a23-a34};

	b.add(v);
	v = {0,-a12,-a12-a23,-a12-a24};

	b.add(v);
	v = {0,-a12,a34-a14,-a14};

	b.add(v);
	v = {0,-a12,-a12-a23,-a14};

	b.add(v);
	v = {0,-a12,-a13,-a13-a34};

	b.add(v);
	v = {0,-a12,-a13,-a12-a24};

	b.add(v);
	v = {0,-a12,-a13,-a14};

	b.add(v);

	b.divideBy(16.0);
	return b.getData();
}

bool Matrix::testAvgSpanTreeParetoOptimal()const{
		std::vector<double> v = getMeanOfSpans();

		return Matrix::testVectorParetoOptimal(v);
}
