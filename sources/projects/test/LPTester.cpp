#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <math.h>

#include <lemon/lp.h>

#include "Matrix.hpp"
#include "MatrixCollection.hpp"
#include "LpSolution.hpp"

int main() {
	Matrix<4> M1({6, 5, 4, 9, 0, 9}); //It is the main example in the article
	std::vector<double> w = {0.64737128, 0.09361541, 0.11588291 , 0.14313040}; //An example vector testing its consistency
	Matrix<4> M2({0,5,5,5,5,0}); //1,6,6,6,6,1 - note, that it has CR = 0 !
	Matrix<4> M3({0, 0, 2, 1, 1, 4}); //1 1	3	2	2	5
	Matrix<4> M4({0, 3, 8, 6, 4, 3});
	std::vector<double> w4 = {0.404518, 0.436173, 0.110295, 0.049014};

	Matrix<6> M6(
	{
		3,  2,  0,  2,  3,
		    6,  2, 12,  0,
			     12, 12, 13,
					      0, 10,
						 	 	    2
	});

	/*std::cout << M6 << std::endl;
	std::cout << M6.testParetoOptimality(filterType::CosineMethod) << std::endl;*/

	Matrix<4> buggyMatrix(
	{
			0, 0, 4,
			   4, 7,
				    2
	});

	std::vector<double> span = buggyMatrix.getMeanOfSpans();

	for (size_t i = 0; i < span.size(); i++) std::cout << span[i] << " ";
	std::cout << std::endl;



	//testing input/output of matrixcollection
	/*
	std::cout << M1 << std::endl;
	MatrixCollection<4> ex;
	ex.add(M1);
	std::cout << ex.size() << "\n";
	std::cout << "Hello!\n";
	ex.saveToFile("../res/MatrixTest.mt");
	std::cout << "It does not look so great!\n";

	MatrixCollection<4> re = MatrixCollection<4>::readFromFile("../res/MatrixTest.mt");
	std::cout << "Heyyy\n";
	std::cout << re[0] << std::endl;
	*/

//Testing rerun:
	/*std::cout << (M4.testPrimalEigenvectorIsParetoOptimal()?"true":"false") << std::endl;
	LpSolution e = M4.LPVectorParetoOptimal(w4);
	e.printOutData();
	e.getOtherTwoVector();

	std::cout << (M1.testPrimalEigenvectorIsParetoOptimal()?"true":"false") << std::endl;
	LpSolution e2 = M1.LPVectorParetoOptimal(w);
	e2.printOutData();
	e2.getOtherTwoVector();*/




	/*std::vector<double> x = e.getxnorm();
	for (int i = 0; i < 4; i++) std::cout << x[i] << " ";
	std::cout << "\n";*/

	//measuring time:
	clock_t begin,end,elapsed_secs;

	/*
	std::vector<double> t;
	t = A.getMeanOfSpans();
	std::cout << std::endl;
	for (int i = 0; i < 4; i++) std::cout << t[i] << " ";
	std::cout << std::endl;
	*/

	//std::cout << A.LPVectorParetoOptimal(w) << std::endl;
	//std::cout << B.testPrimalEigenvectorIsParetoOptimal() << std::endl;


	//testing whether a vector is Pareto optimal
	/*
	begin = clock();
	A.testVectorParetoOptimal(w);
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "" << elapsed_secs << std::endl;
	*/


	/*
	Matrix P(0, 0, 2, 1, 1, 4);
	std::vector<double> p = {1, 1.16084648844746, 1.0267302695589, 0.34671572217346};
	std::vector<double> p2 = {1, 1.13062458843216, 1.00000000000256, 0.337689198861405};

	std::cout << P.testVectorParetoOptimal(p) << std::endl;
	std::cout << P.testVectorParetoOptimal(p2) << std::endl;
	std::cout << P.testPrimalEigenvectorIsParetoOptimal() << std::endl;
	*/


	/*
	begin = clock();
	MatrixCollection C = MatrixCollection::readFromFile("../res/weaklyEffEigens.mt");
	std::cout << "The size is " << C.size() << std::endl;
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "It takes " << elapsed_secs << " secs." << std::endl;
	*/


	begin = clock();
	MatrixCollPtr<4> C = MatrixCollection<4>::readFromFile("../res/consistents.mt");
	std::cout << "The size is " << C->size() << " (consistents)" << std::endl;
	size_t x = C->size();
	C->applyFilter(filterType::Consistency);
	std::cout << "The size is " << C->size() << " after double check." << std::endl;
	if (x != C->size()) throw "OMG";

	C = C->applyFilter(filterType::CosineMethod);
	std::cout << "The size is " << C->size() << " after applying cosine method filter" << std::endl;
	//C.saveToFile("../res/consistentCosines.mt");
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "It takes " << elapsed_secs << " secs." << std::endl;
	/*
	*/

	//Writing en excel file:
	/*MatrixCollection MC = MatrixCollection::readFromFile("../res/consistents.mt");
	std::cout << MC.size() << std::endl;
	MC.generateCsv("../res/consistentsCsv.csv");
	std::cout << "File has been generated to ../res/consistentsCsv.csv" << std::endl;
	*/

	// GENERATING CSV
/*
	MatrixCollection C = MatrixCollection::readFromFile("../res/weaklyEffsBySpanTree.mt");
	std::cout << C.size() << std::endl;
	C.generateCsv("../res/avgSpanTree.csv");
*/

	//MAIN PART:

  /*
	begin = clock();
	MatrixCollection C = MatrixCollection::readFromFile("res/weaklyEffsByEig.mt");
	std::cout << C.size() << std::endl;
	//MatrixCollection weakeffs = C.applyEigenvalueMethodFilter();
	//std::cout << weakeffs.size() << std::endl;
	C.generateCsv("res/eigenvalNonEff.csv");

	std::cout << weakeffs.size() << std::endl;
	weakeffs.saveToFile("../res/weaklyEffsByEig.mt");
	weakeffs.generateCsv("../res/eigenvalEff.csv");

	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Filtering the matrices takes " << elapsed_secs << " secs." << std::endl;
*/
	return 0;
}
