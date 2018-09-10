#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>
#include <ctime>

#include <lemon/lp.h>

#include "Matrix.hpp"
#include "MatrixGenerator.hpp"
#include "MatrixCollection.hpp"
#include "Paths.hpp"

int main(int argc, char* argv[]) {
	//It is always important, to set this:
	clock_t begin,end;
	double elapsed_secs;

	begin = clock();

	try {
		if (argc != 2) {
	    std::cerr << "Usage: " << argv[0] << " [parameter] " << std::endl;
	    throw "Invalid parameters for the program!\n";
	  }
		std::string s = std::string(argv[1]);

		if (s == "ALL4") {
			matrixInit::generateAllToFile<4>();
		}
		else if (s == "SORT") {
			MatrixCollPtr<4> m = MatrixCollection<4>::readFromFile(PATH_4_ALL);
			m->sort();
			m->saveToFile(PATH_4_ALL);
			std::cout << "The number of 4x4 filtered matrices: " << m->size() << std::endl;
		}
		else if (s == "CONS") {
			MatrixCollPtr<4> m = MatrixCollection<4>::readFromFile(PATH_4_ALL);
			std::cout << "Matrices are read.\n";
			MatrixCollPtr<4> m1 = m->applyFilter(filterType::Consistency);
			std::cout << "The number of consistent (<0.1) matrices: " << m1->size() << "\n";
			m1->saveToFile(PATH_4_CONSISTENT);
			std::cout << "Matrices saved.\n";
		}
		else if (s == "FILTER") {
			MatrixCollPtr<4> m = MatrixCollection<4>::readFromFile(PATH_4_CONSISTENT);
			std::cout << "Matrices are read.\n";
			MatrixCollPtr<4> m2 = m->applyFilter(filterType::EigenVectorMethod);
			m2->saveToFile(PATH_4_EIGEN);
			std::cout << "The number of inconsistent matrices by eigenvector method: " << m2->size() << "\n";
			MatrixCollPtr<4> m3 = m->applyFilter(filterType::AverageSpanTreeMethod);
			m3->saveToFile(PATH_4_SPANTREE);
			std::cout << "The number of inconsistent matrices by avg spantree method: " << m3->size() << "\n";
			MatrixCollPtr<4> m4 = m->applyFilter(filterType::CosineMethod);
			m4->saveToFile(PATH_4_COSINE);
			std::cout << "The number of inconsistent matrices by cosine method: " << m4->size() << "\n";
		}
		else if (s == "CSV") {
			std::cout << "Generating csvs...\n";
			MatrixCollPtr<4> meigen = MatrixCollection<4>::readFromFile(PATH_4_EIGEN);
			meigen->regularize();
			meigen->sort();
			meigen->generateCsv(PATH_4_EIGEN_CSV, filterType::EigenVectorMethod);
			std::cout << PATH_4_EIGEN_CSV << " generated.\n";

			MatrixCollPtr<4> mst = MatrixCollection<4>::readFromFile(PATH_4_SPANTREE);
			mst->regularize();
			mst->sort();
			mst->generateCsv(PATH_4_SPANTREE_CSV, filterType::AverageSpanTreeMethod);
			std::cout << PATH_4_SPANTREE_CSV << " generated.\n";

			MatrixCollPtr<4> mcos = MatrixCollection<4>::readFromFile(PATH_4_COSINE);
			mcos->regularize();
			mcos->sort();
			mcos->generateCsv(PATH_4_COSINE_CSV, filterType::CosineMethod);
			std::cout << PATH_4_COSINE_CSV << " generated.\n";
		}
		else if (s == "ALL5") {
			std::vector<bool> v;
			std::cout << v.max_size() << std::endl;
			std::cout << 84113665016 << std::endl;
			throw "Insufficient funds.\n";

			matrixInit::generateAllToFile<5>();
		}
		else if (s == "CSVALL") {
			MatrixCollPtr<4> m = MatrixCollection<4>::readFromFile(PATH_4_CONSISTENT);
			std::cout << "Matrices are read.\n";
			std::cout << "Sorting the data...\n";
			m->regularize();
			m->sort();
			std::cout << "Printing csv...\n";
			m->printCSVWithAllData();
		}
		else {
			std::cerr << "Parameter not recognised.\n" << std::endl;
		}
	} catch (const char * e) {
		//TODO https://stackoverflow.com/questions/77005/how-to-generate-a-stacktrace-when-my-gcc-c-app-crashes
		std::cerr << e << std::endl;
	}

	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	std::cout << "Procedure finished in " << elapsed_secs << " secs.\n";
	return 0;
}
