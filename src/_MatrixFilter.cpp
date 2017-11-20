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
	if (s == "ALL4") {
		matrixInit::generateAllToFile<4>();
	}
	else if (s == "COUNT4") {
		MatrixCollection<4> m = MatrixCollection<4>::readFromFile(PATH_4_ALL);
		std::cout << "The number of 4x4 filtered matrices: " << m.size() << std::endl;
	}
	else if (s == "CONS") {
		MatrixCollection<4> m = MatrixCollection<4>::readFromFile(PATH_4_ALL);
		std::cout << "Matrices are read.\n";
		MatrixCollection<4> m1 = m.applyFilter(filterType::Consistency);
		std::cout << "The number of consistent (<0.1) matrices: " << m1.size() << "\n";
		m1.saveToFile(PATH_4_CONSISTENT);
		std::cout << "Matrices saved.\n";
	}
	else if (s == "FILTER") {
		MatrixCollection<4> m = MatrixCollection<4>::readFromFile(PATH_4_CONSISTENT);
		std::cout << "Matrices are read.\n";
		MatrixCollection<4> m2 = m.applyFilter(filterType::EigenVectorMethod);
		m2.saveToFile(PATH_4_EIGEN);
		std::cout << "The number of consistent matrices by eigenvector method: " << m2.size() << "\n";
		MatrixCollection<4> m3 = m.applyFilter(filterType::AverageSpanTreeMethod);
		m3.saveToFile(PATH_4_SPANTREE);
		std::cout << "The number of consistent matrices by avg spantree method: " << m3.size() << "\n";
		MatrixCollection<4> m4 = m.applyFilter(filterType::CosineMethod);
		m4.saveToFile(PATH_4_COSINE);
		std::cout << "The number of consistent matrices by cosine method: " << m4.size() << "\n";
	}
	else if (s == "CSV") {
		std::cout << "Generating csvs...\n";
		MatrixCollection<4> meigen = MatrixCollection<4>::readFromFile(PATH_4_EIGEN);
		meigen.generateCsv(PATH_4_EIGEN_CSV, filterType::EigenVectorMethod);
		std::cout << PATH_4_EIGEN_CSV << " generated.\n";

		MatrixCollection<4> mst = MatrixCollection<4>::readFromFile(PATH_4_SPANTREE);
		mst.generateCsv(PATH_4_SPANTREE_CSV, filterType::AverageSpanTreeMethod);
		std::cout << PATH_4_SPANTREE_CSV << " generated.\n";

		MatrixCollection<4> mcos = MatrixCollection<4>::readFromFile(PATH_4_COSINE);
		mcos.generateCsv(PATH_4_COSINE_CSV, filterType::CosineMethod);
		std::cout << PATH_4_COSINE_CSV << " generated.\n";
	}
	else if (s == "ALL5") {
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
