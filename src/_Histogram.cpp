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
  //Checking the consistency ratio for all the matrices:
  MatrixCollection<4> m = MatrixCollection<4>::readFromFile(PATH_4_ALL);
  std::cout << "Read finished! The number of matrices is " << m.size() <<"\n";

  std::vector<double> ratiosAll;
  std::vector<bool> ratiosEigen,ratiosSpantree,ratiosCosine;
  std::vector<Ush> numberOfPareto;

  for (size_t i = 0; i < m.size(); i++) {
    if (i % 1007 == 0) std::cout << std::setprecision(2) << ((double(i)) / m.size() * 100) << "%\n";
    ratiosAll.push_back(m[i].getConsistencyRatio());
    ratiosEigen.push_back(m[i].testParetoOptimality(filterType::EigenVectorMethod));
    ratiosSpantree.push_back(m[i].testParetoOptimality(filterType::AverageSpanTreeMethod));
    ratiosCosine.push_back(m[i].testParetoOptimality(filterType::CosineMethod));
    numberOfPareto.push_back(0);//m[i].countParetoVectorsByAlgorithm());
  }

  std::cout << "Consistencies calulcated.\n";
  std::cout << "Creating histogram...\n";

  double steps = 0.091;
  double maxConsistency = 3.8;
  //fmod
  std::vector<unsigned int> bucketsAll((int)floor(maxConsistency / steps));
  std::vector<unsigned int> bucketsEigen((int)floor(maxConsistency / steps));
  std::vector<unsigned int> bucketsSpantree((int)floor(maxConsistency / steps));
  std::vector<unsigned int> bucketsCommon((int)floor(maxConsistency / steps));
  std::vector<unsigned int> bucketsCosine((int)floor(maxConsistency / steps));
  std::vector<unsigned int> bucketsOfPareto((int)floor(maxConsistency / steps));

  std::vector<double> outputAll((int)floor(maxConsistency / steps));
  std::vector<double> outputEigen((int)floor(maxConsistency / steps));
  std::vector<double> outputSpantree((int)floor(maxConsistency / steps));
  std::vector<double> outputCommon((int)floor(maxConsistency / steps));
  std::vector<double> outputCosine((int)floor(maxConsistency / steps));
  std::vector<double> outputAvgOfPareto((int)floor(maxConsistency / steps));

  for (size_t i = 0; i < bucketsAll.size(); i++) bucketsAll[i] = 0;
  for (size_t i = 0; i < bucketsEigen.size(); i++) bucketsEigen[i] = 0;
  for (size_t i = 0; i < bucketsSpantree.size(); i++) bucketsSpantree[i] = 0;
  for (size_t i = 0; i < bucketsCommon.size(); i++) bucketsCommon[i] = 0;
  for (size_t i = 0; i < bucketsCosine.size(); i++) bucketsCosine[i] = 0;
  for (size_t i = 0; i < bucketsOfPareto.size(); i++) bucketsOfPareto[i] = 0;

  for (size_t i = 0; i < ratiosAll.size(); i++) {
    if (ratiosAll[i] < maxConsistency) {
      if ((unsigned int)floor(ratiosAll[i] / steps) <= bucketsAll.size() - 1) {
        bucketsAll[(int)floor(ratiosAll[i] / steps)]++;
        if (!ratiosEigen[i]) bucketsEigen[(int)floor(ratiosAll[i] / steps)]++;
        if (!ratiosSpantree[i]) bucketsSpantree[(int)floor(ratiosAll[i] / steps)]++;
        if (!ratiosEigen[i] && !ratiosSpantree[i] && !ratiosCosine[i]) bucketsCommon[(int)floor(ratiosAll[i] / steps)]++;
        if (!ratiosCosine[i]) bucketsCosine[(int)floor(ratiosAll[i] / steps)]++;
        bucketsOfPareto[(int)floor(ratiosAll[i] / steps)] += numberOfPareto[i];
      }
    }
  }

  for (size_t i = 0; i < bucketsAll.size(); i++) {
    outputAll[i] = bucketsAll[i];
    if (bucketsAll[i] != 0) {
      outputEigen[i] = ((double)bucketsEigen[i]) / bucketsAll[i];
      outputSpantree[i] = ((double)bucketsSpantree[i]) / bucketsAll[i];
      outputCommon[i] = ((double)bucketsCommon[i]) / bucketsAll[i];
      outputCosine[i] = ((double)bucketsCosine[i]) / bucketsAll[i];
      outputAvgOfPareto[i] = ((double)bucketsOfPareto[i]) / bucketsAll[i];
    }
  }

  std::cout << "Bucketing done, writing file...\n";

  std::ofstream I("../res/basicHistogram.csv");
  //Header:
  I << "interval\t#all\t%eigen\t%spantree\t%cosine\t%common\t\n";
  for (size_t i = 0; i < bucketsAll.size(); i++) {
    I << i * steps << "\t" << outputAll[i] << "\t"
      << outputEigen[i] << "\t" << outputSpantree[i] << "\t" << outputCosine[i] << "\t" << outputCommon[i] << "\t"
      << outputAvgOfPareto[i] << std::endl;
  }
  I.close();
  std::cout << "Procedure finished.\n";

  return 0;
}
