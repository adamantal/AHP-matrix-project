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

//TODO: delete this function
template<size_t N>
Matrix<N>& MatrixCollection<N>::operator[](size_t ind) {
	return data.at(ind);
}

template<size_t N>
Matrix<N>& MatrixCollection<N>::at(size_t ind) {
	return data.at(ind);
}

template<size_t N>
size_t MatrixCollection<N>::size(){
	return data.size();
}

template<size_t N>
bool MatrixCollection<N>::saveToFile(std::string filename) {
	//regularize ();

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
MatrixCollPtr<N> MatrixCollection<N>::readFromFile(std::string filename) {
	if (N != 4) throw "The readFromFile method has not been implemented for 5 or greater dimensions.";
	std::ifstream I;
	I.open(filename);

	if (!I.is_open()) {
		std::cout << "File can not be opened, please check file!\n";
		throw "FILE CAN NOT BE OPENED.";
	}

	MatrixCollPtr<N> mc = std::make_shared<MatrixCollection<N>>();
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
		mc->add(tmpm);
	}
	I.close();

	return mc;
}

template<size_t N>
typename std::vector< Matrix<N> >::iterator MatrixCollection<N>::begin() { return data.begin(); }

template<size_t N>
typename std::vector< Matrix<N> >::iterator MatrixCollection<N>::end() { return data.end(); }

template<size_t N>
MatrixCollPtr<N> MatrixCollection<N>::applyFilter(filterType filter) {
	MatrixCollPtr<N> tmp = std::make_shared<MatrixCollection<N>>();
	switch(filter) {
		case (filterType::Inconsistency) :
			{
				for (auto it = data.begin(); it != data.end(); it++) {
					double consistencyRatio = it->getConsistencyRatio();
					if (consistencyRatio > 0.1) tmp->add(*it);
				}
			}
			break;
		case (filterType::Consistency) :
			{
				for (auto it = data.begin(); it != data.end(); it++) {
					double consistencyRatio = it->getConsistencyRatio();
					if (consistencyRatio <= 0.1) tmp->add(*it);
				}
			}
			break;
		default :
			{
				for (auto it = data.begin(); it != data.end(); it++) {
					if (!(it->testParetoOptimality(filter))) tmp->add(*it);
				}
			}
			break;
	}
	return tmp;
}

template<>
void MatrixCollection<4>::generateCsv(std::string filename, filterType filter) {
	std::ofstream F;
	F.open(filename);
	//TODO: header?
	if (filter == filterType::EigenVectorMethod) {
		F << "Index, Matrix, , , , , , Eig.val, CR, Eig.vec, , , , , LP. vec, , , , , Vec2, , , , , Vec3, , , , , \n";
	} else {
		F << "Index, Matrix, , , , , , Eig.val, CR, Eig.vec, , , , , Found vector, , , , , LP. vec, , , , , Vec2, , , , , Vec3, , , , , \n";
	}

	for (auto it = data.begin(); it != data.end(); it++){
		F << std::defaultfloat << std::distance(data.begin(), it) + 1 << ", "
			<< it->get(0,1) << ", " << it->get(0,2) << ", " << it->get(0,3) << ", " << it->get(1,2) << ", " << it->get(1,3) << ", " << it->get(2,3) << ", ";
		F << std::fixed << std::setprecision(10) <<  it->largestEigenvalue() << ", "
			<< it->getConsistencyRatio() << ", , "
			<< it->getPrimalNormEigenvector()[0] << ", " << it->getPrimalNormEigenvector()[1] << ", " << it->getPrimalNormEigenvector()[2] << ", " << it->getPrimalNormEigenvector()[3] << ", , ";

				std::vector<double> V;
				switch (filter) {
					case (filterType::EigenVectorMethod):
						V = it->getPrimalNormEigenvector();
						break;

					case (filterType::CosineMethod):
						V = it->getCosineVector();
						break;
					case (filterType::AverageSpanTreeMethod):
						V = it->getMeanOfSpans();
						break;
					default :
						throw "For this filter the CSV generator file has not yet been implemented.\n";
				}
				if (filter != filterType::EigenVectorMethod) {
					for (size_t i = 0; i < 4; i++) {
						F << V[i] << ", ";
					}
					F << ", ";
				}
				LpSolution<4> lp = it->LPVectorParetoOptimal(V);
				std::vector<double> vec = lp.getxnorm();
		F << vec[0] << ", " << vec[1] << ", " << vec[2] << ", " << vec[3] << ", , ";
				std::vector<double> oth = lp.getOtherTwoVector();
		F << oth[0] << ", " << oth[1] << ", " << oth[2] << ", " << oth[3] << ", , "
			<< oth[4] << ", " << oth[5] << ", " << oth[6] << ", " << oth[7]
			<< "\n";
	}
	F.close();
}

template<>
void MatrixCollection<4>::printCSVWithAllData() {
	std::ofstream F;
	F.open("../res/all4x4OptimalityData.csv");
	std::vector<filterType> filterTypes = {filterType::EigenVectorMethod, filterType::AverageSpanTreeMethod, filterType::CosineMethod};
	//header:
	F << "#,Matrix, , , , , ,Larg.eig.value,Efficiency,Methods:,Eigenvector,Spantree,Cosine,Vectors:,Eigenvector, , , , ,Spantree, , , , ,Cosine, , , ,LP return if inefficient:,Eigenvector, , , , , Spantree, , , , ,Cosine , , , , ,\n";

	for (auto it = data.begin(); it != data.end(); it++) {
		if ((std::distance(data.begin(), it) + 1) % 1000 == 0) {
			std::cout << "[" << data.size() << "//" << std::distance(data.begin(), it) + 1 << "]" << std::endl;
		}

		//base data:
		F << std::defaultfloat << std::distance(data.begin(), it) + 1 << ", "
			<< it->get(0,1) << ", " << it->get(0,2) << ", " << it->get(0,3) << ", " << it->get(1,2) << ", " << it->get(1,3) << ", " << it->get(2,3) << ", "
			<< std::fixed << std::setprecision(10) << it->largestEigenvalue() << ", "
			<< it->getConsistencyRatio() << ", , ";

		for (auto const& filter : filterTypes) {
			//for all filtering type test the matrix
			if (it->testParetoOptimality(filter)) {
				F << "Y, ";
			} else {
				F << "N, ";
			}
		}
		std::vector<std::vector<double>> vectors;
		vectors.push_back(it->getPrimalNormEigenvector());
		vectors.push_back(it->getMeanOfSpans());
		vectors.push_back(it->getCosineVector());
		std::for_each(vectors.begin(), vectors.end(), &Matrix<4>::L1);

		F << ", ";
		for (const auto & vec : vectors) {
			for (auto it = vec.begin(); it != vec.end(); it++) {
				F << *it << ", ";
			}
			F << ", ";
		}

		assert(vectors.size () == filterTypes.size());
		for (size_t i = 0; i < filterTypes.size(); i++) {
			if (!(it->testParetoOptimality(filterTypes[i]))) {
				LpSolution<4> lp = it->LPVectorParetoOptimal(vectors[i]);
				std::vector<double> opt = lp.getxnorm();
				for (size_t j = 0; j < opt.size(); j++) {
					F << opt[j] << ", ";
				}
				F << ", ";
			} else {
				F << ", , , , , ";
			}
		}

		F << "\n";
	}
}

template<size_t N>
bool MatrixCollection<N>::isIncluded(const Matrix<N>& M, unsigned long long int& j)const{
	//WE ASSUME THE COLLECTION IS ORDERED
	size_t begin, end;
	begin = 0;
	end = data.size() - 1;
	while (begin != end) { //logarithmic search
		size_t middle = (begin + end) / 2;
		if (begin == middle || middle == end) {
			if (M == data[begin]) {
				j = begin;
				return true;
			} else if (M == data[end]) {
				j = end;
				return true;
			} else {
				break;
			}
		}
		if (M < data[middle]) {
			end = middle;
		} else {
			begin = middle;
		}
	}
	return false;
}

template<size_t N>
void MatrixCollection<N>::regularize () {
	for (auto i = data.begin(); i != data.end(); i++) {
		i->regularize ();
	}
}

template<size_t N>
void MatrixCollection<N>::sort () {
	std::sort(data.begin(), data.end(),
		[] (const Matrix<N> & a, const Matrix<N> & b) -> bool
			{
				Ush aT = a.countIndexOfUpperTriangle ();
				Ush bT = b.countIndexOfUpperTriangle ();
				if (aT != bT) {
					return aT > bT;
				} else {
					return a < b;
				}
			}
	);
}
