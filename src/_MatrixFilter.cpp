#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>
#include <ctime>

#include <lemon/lp.h>

#include "Matrix.cpp"
#include "MatrixCollection.cpp"
#include "generateAll4.cpp"

/*
//Generating all 4x4 matrices:
void generatingAll4x4Matrices() {
	begin = clock();
	matrixInit::generateAllToFile();
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Generating the file takes " << elapsed_secs << " secs." << std::endl;
}

void readAll4x4Matrcies(){
	begin = clock();
	MatrixCollection m = matrixInit::readAllFromFile("../res/allMatrices.mt");
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Reading the file takes " << elapsed_secs << " secs." << std::endl;
}

void filteringMatricesByConsistencyRatio(){
	//Checking the consistency ratio for all the matrices:
	begin = clock();
	MatrixCollection coll;

	for (int i = 0; i < m.size(); i++){
		if (m[i].getConsistencyRatio() < 0.1) {
			coll.add(m[i]);
		}
	}
	coll.saveToFile("consistents.mt");
	std::cout << coll.size() << "\n";
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Writing and getting all consistency in " << elapsed_secs << " secs." << std::endl;
}
*/
int main() {
	//It is always important, to set this:
	matrixInit::setElement();
	clock_t begin,end;
	double elapsed_secs;

	//Checking the consistency ratio for all the matrices:
	begin = clock();
	MatrixCollection m = MatrixCollection::readFromFile("../res/allMatrices.mt");
	std::cout << "Read!\n";
	std::vector<double> ratios(m.size());

	for (size_t i = 0; i < m.size(); i++) {
		ratios[i] = m[i].getConsistencyRatio();
	}
	std::cout << "Done!\n";

	return 0;
}
