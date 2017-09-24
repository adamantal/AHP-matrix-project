#include "Matrix.h"
#include "MatrixCollection.h"

MatrixCollection::MatrixCollection():
	data(std::vector<Matrix>()),
	tmpindex(0){}

void MatrixCollection::add(Matrix &m){
	data.push_back(m);
}

Matrix& MatrixCollection::operator[](size_t ind){
	return data[ind];
}

size_t MatrixCollection::size(){
	return data.size();
}

bool MatrixCollection::saveToFile(std::string filename){
	std::ofstream F;
	F.open(filename);
	size_t c = 0;
	for (auto i = data.begin(); i != data.end(); i++) {
		F << "#" << c++ << std::endl;
		F << i->toIndexString() << std::endl;
	}
	F.close();
	return true;
}

MatrixCollection MatrixCollection::readFromFile(std::string filename) {
	std::ifstream I;
	I.open(filename);

	if (!I.is_open()) throw "FILE CAN NOT BE OPENED.";

	MatrixCollection mc;
	std::string x;

	while (I >> x) {
		std::vector<int> vv;
		std::string tmp1;
		I >> tmp1;
		for (int i = 0; i < 3; i++) {
			std::string cha;
			I >> cha;
			int tmpint = std::stoi(cha);
			vv.push_back(tmpint);
		}
		for (int i = 0; i < 2; i++) {
			std::string cha;
			I >> cha;
		}
		for (int i = 0; i < 2; i++) {
			std::string cha;
			I >> cha;
			int tmpint = std::stoi(cha);
			vv.push_back(tmpint);
		}
		for (int i = 0; i < 3; i++) {
			std::string cha;
			I >> cha;
		}
		I >> tmp1;
		int tmpint = std::stoi(tmp1);
		vv.push_back(tmpint);
		for (int i = 0; i < 4; i++) {
			std::string cha;
			I >> cha;
		}
		Matrix m = Matrix(vv);
		mc.add(m);
	}
	I.close();

	return mc;
}

typedef std::vector<Matrix>::iterator iterator;

iterator MatrixCollection::begin() { return data.begin(); }
iterator MatrixCollection::end() { return data.end(); }

MatrixCollection MatrixCollection::applyInconsistencyFilter(){
	MatrixCollection tmp;
	for (auto it = data.begin(); it != data.end(); it++) {
		double consistencyRatio = it->getConsistencyRatio();
		if (consistencyRatio > 0.1) tmp.add(*it);
	}
	return tmp;
}
MatrixCollection MatrixCollection::applyConsistencyFilter(){
	MatrixCollection tmp;
	for (auto it = data.begin(); it != data.end(); it++) {
		//if (std::distance(data.begin(), it) % 1000 == 0) std::cout << std::distance(data.begin(), it) << std::endl;
		double consistencyRatio = it->getConsistencyRatio();
		if (consistencyRatio <= 0.1) tmp.add(*it);
	}
	return tmp;
}
MatrixCollection MatrixCollection::applyEigenvalueMethodFilter(){
	MatrixCollection tmp;
	for (auto it = data.begin(); it != data.end(); it++) {
		if (!(it->testPrimalEigenvectorIsParetoOptimal())) tmp.add(*it);
	}
	return tmp;
}
MatrixCollection MatrixCollection::applyAvgSpanTreeFilter(){
	MatrixCollection tmp;
	for (auto it = data.begin(); it != data.end(); it++) {
		if (!(it->testAvgSpanTreeParetoOptimal())) tmp.add(*it);
	}
	return tmp;
}
void MatrixCollection::generateCsv(std::string filename){
	std::ofstream F;
	F.open(filename);
	//TODO: header?
	F << "Index, Matrix, , , , , , Eig.val, CR, Eig.vec, , , , , LP. vec, , , , , Vec2, , , , , Vec3, , , , , \n";

	for (auto it = data.begin(); it != data.end(); it++){
		F << std::distance(data.begin(), it) + 1 << ", "
			<< it->get(0,1) << ", " << it->get(0,2) << ", " << it->get(0,3) << ", " << it->get(1,2) << ", " << it->get(1,3) << ", " << it->get(2,3) << ", "
			<< it->largestEigenvalue() << ", "
			<< it->getConsistencyRatio() << ", , "
			<< it->getPrimalNormEigenvector()[0] << ", " << it->getPrimalNormEigenvector()[1] << ", " << it->getPrimalNormEigenvector()[2] << ", " << it->getPrimalNormEigenvector()[3] << ", , ";

				LpSolution lp = it->LPVectorParetoOptimal(it->getPrimalNormEigenvector());
				std::vector<double> vec = lp.getxnorm();
		F << vec[0] << ", " << vec[1] << ", " << vec[2] << ", " << vec[3] << ", , ";
				std::vector<double> oth = lp.getOtherTwoVector();
		F << oth[0] << ", " << oth[1] << ", " << oth[2] << ", " << oth[3] << ", , "
			<< oth[4] << ", " << oth[5] << ", " << oth[6] << ", " << oth[7]
			<< "\n";
	}
	F.close();
}

/*bool MatrixCollection::checkAssumption(){
	for (auto it = data.begin(); it != std::end(data); it++){
		std::cout << "#" << std::distance(data.begin(), it) << std::endl;
		LpSolution lp = it->LPVectorParetoOptimal(it->getPrimalNormEigenvector());
		std::cout << "Hi";
		if (std::distance(data.begin(), it) == 6) {
			lp.printOutData();
		}
		if (lp.getOtherTwoVector() != 2) return false;
	}
	return true;
}*/
