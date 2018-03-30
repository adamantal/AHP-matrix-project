#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

#include <lemon/lp.h>

#include "Matrix.cpp"
#include "MatrixCollection.cpp"
#include "MatrixGenerator.cpp"

/*
This program aims to create a histogram to get average and total number of integers in the upper triangle part of the considered matrices.
*/

int main() {
  double steps = 0.07;
  double maxConsistency = 3.8;
  Ush maxIndex = (int)floor(maxConsistency / steps);

  std::vector<unsigned int> bucketsAll(maxIndex);
  std::vector<unsigned int> bucketsTriang(maxIndex);
  std::vector<unsigned int> bucketsTriang0(maxIndex);
  std::vector<unsigned int> bucketsTriang1(maxIndex);
  std::vector<unsigned int> bucketsTriang2(maxIndex);
  std::vector<unsigned int> bucketsTriang3(maxIndex);

  MatrixCollPtr<4> m = MatrixCollection<4>::readFromFile(PATH_4_ALL);
  std::cout << "Read finished! The number of matrices is " << m->size() <<"\n";
  std::cout << "Regularizing\n";
  m->regularize ();

  std::cout << "Calculating upper triangles\n";
  std::vector<double> ratiosAll;
  std::vector<Ush> numberTriang;

  //max search:
  Matrix<4> maxMatrix0;
  double maxInconsistency0 = 0.0;
  Matrix<4> minMatrix2;
  double minInconsistency2 = 3.9;

  for (size_t i = 0; i < m->size(); i++) {
    if (i % 10070 == 0) {
        std::cout << std::setprecision(2) << ((double(i)) / m->size() * 100) << "%\n";
    }
    ratiosAll.push_back(m->at(i).getConsistencyRatio());
    numberTriang.push_back(6 - m->at(i).countIndexOfUpperTriangle ());
    if (m->at(i).countIndexOfUpperTriangle () == 6 && m->at(i).getConsistencyRatio() > maxInconsistency0) {
      maxMatrix0 = m->at(i);
      maxInconsistency0 = m->at(i).getConsistencyRatio();
    }
    if (m->at(i).countIndexOfUpperTriangle () == 4 && m->at(i).getConsistencyRatio() < minInconsistency2) {
      minMatrix2 = m->at(i);
      minInconsistency2 = m->at(i).getConsistencyRatio();
    }
  }

  std::cout << maxInconsistency0 << std::endl << maxMatrix0 << std::endl;
  std::cout << minInconsistency2 << std::endl << minMatrix2 << std::endl;

  std::cout << "Creating histogram...\n";

  if (ratiosAll.size() != numberTriang.size())
    throw "Sizes mismtach";

  std::cout << ratiosAll.size() << "\n";
  for (size_t i = 0; i < ratiosAll.size(); i++) {
    if (ratiosAll.at(i) < maxConsistency) {
      unsigned int index = (unsigned int)floor(ratiosAll.at(i) / steps);
      if (index < bucketsAll.size()) {
        bucketsAll.at(index) += 1;
        bucketsTriang.at(index) += numberTriang.at(i);
        switch (numberTriang.at(i)) {
            case 0:
                bucketsTriang0.at(index) += 1;
            break;
            case 1:
                bucketsTriang1.at(index) += 1;
            break;
            case 2:
                bucketsTriang2.at(index) += 1;
            break;
            case 3:
                bucketsTriang3.at(index) += 1;
            break;
            default:
                throw "ERROR!";
            break;
        }
      }
    }
  }

  std::cout << "Bucketing\n";

  std::vector<double> outputAll(maxIndex);
  std::vector<double> outputTriang(maxIndex);
  std::vector<double> outputTriang0(maxIndex);
  std::vector<double> outputTriang1(maxIndex);
  std::vector<double> outputTriang2(maxIndex);
  std::vector<double> outputTriang3(maxIndex);

  for (size_t i = 0; i < bucketsAll.size(); i++) {
    outputAll.at(i) = bucketsAll.at(i);
    if (bucketsAll.at(i) != 0) {
      outputTriang.at(i) = ((double)bucketsTriang.at(i)) / bucketsAll.at(i);
      outputTriang0.at(i) = ((double)bucketsTriang0.at(i)) / bucketsAll.at(i);
      outputTriang1.at(i) = ((double)bucketsTriang1.at(i)) / bucketsAll.at(i);
      outputTriang2.at(i) = ((double)bucketsTriang2.at(i)) / bucketsAll.at(i);
      outputTriang3.at(i) = ((double)bucketsTriang3.at(i)) / bucketsAll.at(i);
    } else {
      outputTriang.at(i) = 0.0;
    }
  }

  std::cout << "Bucketing done, writing file\n";

  std::ofstream iss("../res/basicHistogram2.csv");
  //Header:
  iss << "consistency\t#all\tavg upper\t6 upper\t5 upper\t4 upper\t3 upper\n";
  for (size_t i = 0; i < bucketsAll.size(); i++) {
    iss << i * steps << "\t" << outputAll.at(i) << "\t" <<
        outputTriang.at(i) << "\t" << outputTriang0.at(i) << "\t" << outputTriang1.at(i) << "\t" << outputTriang2.at(i) << "\t" << outputTriang3.at(i) << "\n";
  }
  iss.close();
  std::cout << "Procedure finished.\n";

  return 0;
}
