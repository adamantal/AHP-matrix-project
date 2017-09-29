#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <ctime>

#include <lemon/lp.h>

#include "Matrix.cpp"
#include "MatrixCollection.cpp"
#include "MatrixGenerator.cpp"

int main() {
  matrixInit::setElement();

  //Checking the consistency ratio for all the matrices:
  MatrixCollection m = MatrixCollection::readFromFile("../res/allMatrices.mt");
  std::cout << "Read finished! The number of matrices is " << m.size() <<"\n";

  std::vector<double> ratiosAll;
  std::vector<bool> ratiosEigen,ratiosSpantree,ratiosCosine;

  for (size_t i = 0; i < m.size(); i++) {
    if (i % 10070 == 0) std::cout << ((double(i)) / m.size() * 100) << "%\n";
    ratiosAll.push_back(m[i].getConsistencyRatio());
    ratiosEigen.push_back(m[i].testPrimalEigenvectorIsParetoOptimal());
    ratiosSpantree.push_back(m[i].testAvgSpanTreeParetoOptimal());
    ratiosCosine.push_back(m[i].testCosineParetoOptimal());
  }

  std::cout << "Consistencies calulcated.\n";
  std::cout << "Creating histogram...\n";

  double steps = 0.01;
  double maxConsistency = 1.0;
  //fmod
  std::vector<unsigned int> bucketsAll((int)floor(maxConsistency / steps));
  std::vector<unsigned int> bucketsEigen((int)floor(maxConsistency / steps));
  std::vector<unsigned int> bucketsSpantree((int)floor(maxConsistency / steps));
  std::vector<unsigned int> bucketsCommon((int)floor(maxConsistency / steps));
  std::vector<unsigned int> bucketsCosine((int)floor(maxConsistency / steps));

  for (size_t i = 0; i < bucketsAll.size(); i++) bucketsAll[i] = 0;
  for (size_t i = 0; i < bucketsEigen.size(); i++) bucketsEigen[i] = 0;
  for (size_t i = 0; i < bucketsSpantree.size(); i++) bucketsSpantree[i] = 0;
  for (size_t i = 0; i < bucketsCommon.size(); i++) bucketsCommon[i] = 0;
  for (size_t i = 0; i < bucketsCommon.size(); i++) bucketsCosine[i] = 0;

  for (size_t i = 0; i < ratiosAll.size(); i++) {
    if (ratiosAll[i] < maxConsistency) {
      if ((unsigned int)floor(ratiosAll[i] / steps) <= bucketsAll.size() - 1) {
        bucketsAll[(int)floor(ratiosAll[i] / steps)]++;
        if (!ratiosEigen[i]) bucketsEigen[(int)floor(ratiosAll[i] / steps)]++;
        if (!ratiosSpantree[i]) bucketsSpantree[(int)floor(ratiosAll[i] / steps)]++;
        if (!ratiosEigen[i] && !ratiosSpantree[i]) bucketsCommon[(int)floor(ratiosAll[i] / steps)]++;
        if (!ratiosCosine[i]) bucketsCosine[(int)floor(ratiosAll[i] / steps)]++;
      }
    }
  }
  std::cout << "Bucketing done, writing file...\n";

  std::ofstream I("../res/basicHistogram.csv");
  //Header:
  I << "bucket,all,eigen,spantree,common\n";
  for (size_t i = 0; i < bucketsAll.size(); i++) {
    I << i << "," << bucketsAll[i] << "," << bucketsEigen[i] << "," << bucketsSpantree[i] << "," << bucketsCommon[i] << std::endl;
  }
  I.close();
  std::cout << "Procedure finished.\n";

  return 0;
}
