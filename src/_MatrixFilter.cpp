#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>
#include <ctime>

#include <lemon/lp.h>

#include "Matrix.cpp"
#include "MatrixCollection.cpp"
#include "MatrixGenerator.cpp"

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
int main(int argc, char* argv[]) {
	if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [parameter] " << std::endl;
  	return 1;
  }

	//It is always important, to set this:
	clock_t begin,end;
	double elapsed_secs;

	begin = clock();

	std::string s = std::string(argv[1]);
	if (s == "ALL") {
		matrixInit::generateAllToFile<4>();
	}

	else if (s == "CONS") {
		MatrixCollection<4> m = MatrixCollection<4>::readFromFile("../res/allMatrices.mt");
		std::cout << "Matrices are read.\n";
		MatrixCollection<4> m1 = m.applyFilter(filterType::Consistency);
		std::cout << "The number of consistent (<0.1) matrices: " << m1.size() << "\n";
		m1.saveToFile("../res/consistents.mt");
		std::cout << "Matrices saved.\n";
	}

	else if (s == "FILTER") {
		MatrixCollection<4> m = MatrixCollection<4>::readFromFile("../res/consistents.mt");
		std::cout << "Matrices are read.\n";
		MatrixCollection<4> m2 = m.applyFilter(filterType::EigenVectorMethod);
		m2.saveToFile("../res/eigenPareto.mt");
		std::cout << "The number of consistent matrices by eigenvector method: " << m2.size() << "\n";
		MatrixCollection<4> m3 = m.applyFilter(filterType::AverageSpanTreeMethod);
		m3.saveToFile("../res/avgspanPareto.mt");
		std::cout << "The number of consistent matrices by avg spantree method: " << m3.size() << "\n";

		m = MatrixCollection<4>::readFromFile("../res/consistents.mt");
		MatrixCollection<4> m4 = m.applyFilter(filterType::CosineMethod);
		m4.saveToFile("../res/cosinePareto.mt");
		std::cout << "The number of consistent matrices by cosine method: " << m4.size() << "\n";
	}

	else if (s == "5TEST") {
		try {
			matrixInit::generateAllToFile<5>();
		} catch (const char * s) {
			std::cout << s << "\n";
		}
	}

	else {
		std::cerr << "Parameter not recognised.\n" << std::endl;
	}

	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	std::cout << "Procedure finished in " << elapsed_secs << " secs.\n";
	return 0;
}
