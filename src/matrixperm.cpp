#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>
#include <ctime>

#include "Matrix.cpp"
#include "MatrixCollection.cpp"
#include "generateAll4.cpp"

int main() {
	//It is always important, to set this:
	matrixInit::setElement();

	//Generating all 4x4 matrices:
	/*clock_t begin = clock();
	generateAllToFile();
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Generating the file takes " << elapsed_secs << " secs." << std::endl;
	*/


	//Reading file:
	/*clock_t begin = clock();
	MatrixCollection m = matrixInit::readAllFromFile("output.mt");
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Reading the file takes " << elapsed_secs << " secs." << std::endl;*/


	//Checking the consistency ratio for all the matrices:
	/*begin = clock();
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
	std::cout << "Writing and getting all consistency in " << elapsed_secs << " secs." << std::endl;*/

	clock_t begin = clock();
	std::string s = "consistents.mt";
	MatrixCollection coll;
	coll.readAllFromFile(s);
	clock_t end = clock();
	clock_t elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Reading all consistencies in " << elapsed_secs << " secs." << std::endl;

	return 0;
}
