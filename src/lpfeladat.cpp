#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <math.h>

#include <lemon/lp.h>

#include "Matrix.cpp"
#include "MatrixCollection.cpp"
#include "generateAll4.cpp"
#include "lpformatrix.cpp"

//g++ --std=c++14 -o3 -lglpk -lemon file.cpp -o file

int main() {
	matrixInit::setElement();
	Matrix A(6,5,4,9,0,9);
	std::vector<double> w = {0.64737128, 0.09361541, 0.11588291 , 0.14313040};
	
	clock_t begin = clock();
	lpTask(A,w,false);
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Without log it takes " << elapsed_secs << " secs to run the lp." << std::endl;
	
	begin = clock();
	lpTask(A,w,true);
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "With log it takes " << elapsed_secs << " secs to run the lp." << std::endl;
	
	begin = clock();
	std::string fileName = "consistents.mt";
	MatrixCollection mc = matrixInit::readAllFromFile(fileName);
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Matrixes loaded from file in " << elapsed_secs << " secs." << std::endl;
	begin = clock();
	for (int i = 0; i < mc.size(); i++) {
		lpTask(mc[i], w, false);
	}
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Running the LP takes " << elapsed_secs << " secs." << std::endl;
	
	return 0;	
}
