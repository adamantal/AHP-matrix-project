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

MatrixCollection MatrixCollection::readFromFile(std::string filename){
	std::ifstream I;
	I.open(filename);
	
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
