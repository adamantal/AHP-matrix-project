#include "Matrix.h"
#include "MatrixCollection.h"

template<size_t N>
MatrixCollection<N>::MatrixCollection():
	data(std::vector<Matrix<N>>()),
	tmpindex(0){}

template<size_t N>
void MatrixCollection<N>::add(Matrix<N> &m){
	data.push_back(m);
}

template<size_t N>
Matrix<N>& MatrixCollection<N>::operator[](size_t ind) {
	return data[ind];
}

template<size_t N>
size_t MatrixCollection<N>::size(){
	return data.size();
}

template<size_t N>
bool MatrixCollection<N>::saveToFile(std::string filename){
	std::ofstream F;
	F.open(filename);

	if (!F.is_open()) {
		std::cout << "File can not be opened, please check file!\n";
		throw "FILE CAN NOT BE OPENED.";
	}

	size_t c = 0;
	for (auto i = data.begin(); i != data.end(); i++) {
		F << "#" << c++ << std::endl;
		F << i->toString(true) << std::endl;
	}
	F.close();
	return true;
}

template<size_t N>
MatrixCollection<N> MatrixCollection<N>::readFromFile(std::string filename) {
	if (N != 4) throw "The readFromFile method has not been implemented for 5 or greater dimensions.";
	std::ifstream I;
	I.open(filename);

	if (!I.is_open()) {
		std::cout << "File can not be opened, please check file!\n";
		throw "FILE CAN NOT BE OPENED.";
	}

	MatrixCollection<N> mc;
	std::string s;

	while (I >> s) {
		std::vector<Ush> elements;
		for (size_t i = 0; i < N; i++) {
			for (size_t j = 0; j < N; j++) {
				std::string tmp;
				I >> tmp;
				if (i < j) {
					Ush elem = std::stoi(tmp);
					elements.push_back(elem);
				}
			}
		}
		Matrix<N> tmpm(elements);
		mc.add(tmpm);
	}
	I.close();

	return mc;
}

template<size_t N>
typename std::vector< Matrix<N> >::iterator MatrixCollection<N>::begin() { return data.begin(); }

template<size_t N>
typename std::vector< Matrix<N> >::iterator MatrixCollection<N>::end() { return data.end(); }

template<size_t N>
MatrixCollection<N> MatrixCollection<N>::applyFilter(filterType filter) {
	MatrixCollection tmp;
	switch(filter) {
		case (filterType::Inconsistency) :
			{
				for (auto it = data.begin(); it != data.end(); it++) {
					double consistencyRatio = it->getConsistencyRatio();
					if (consistencyRatio > 0.1) tmp.add(*it);
				}
			}
			break;
		case (filterType::Consistency) :
			{
				for (auto it = data.begin(); it != data.end(); it++) {
					double consistencyRatio = it->getConsistencyRatio();
					if (consistencyRatio <= 0.1) tmp.add(*it);
				}
			}
			break;
		case (filterType::EigenVectorMethod) :
		case (filterType::AverageSpanTreeMethod) :
		case (filterType::CosineMethod) :
			{
				for (auto it = data.begin(); it != data.end(); it++) {
					if (!(it->testParetoOptimality(filter))) tmp.add(*it);
				}
			}
			break;
		default:
			throw "Unknown filter applied.\n";
			break;
	}
	return tmp;
}

template<size_t N>
void MatrixCollection<N>::generateCsv(std::string filename){
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

				LpSolution<N> lp = it->LPVectorParetoOptimal(it->getPrimalNormEigenvector());
				std::vector<double> vec = lp.getxnorm();
		F << vec[0] << ", " << vec[1] << ", " << vec[2] << ", " << vec[3] << ", , ";
				std::vector<double> oth = lp.getOtherTwoVector();
		F << oth[0] << ", " << oth[1] << ", " << oth[2] << ", " << oth[3] << ", , "
			<< oth[4] << ", " << oth[5] << ", " << oth[6] << ", " << oth[7]
			<< "\n";
	}
	F.close();
}
