#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

#include <lemon/lp.h>

#include "Matrix.cpp"
#include "MatrixCollection.cpp"
#include "MatrixGenerator.cpp"

int main() {
  double steps = 0.11;
  double maxConsistency = 3.8;
  Ush maxIndex = (int)floor(maxConsistency / steps);

  std::vector<unsigned int> bucketsAll(maxIndex);
  std::vector<unsigned int> bucketsTriang(maxIndex);
  std::vector<unsigned int> bucketsTriang0(maxIndex);
  std::vector<unsigned int> bucketsTriang1(maxIndex);
  std::vector<unsigned int> bucketsTriang2(maxIndex);
  std::vector<unsigned int> bucketsTriang3(maxIndex);

  std::vector<double> outputAll(maxIndex);
  std::cout << "line2.1\n";
  std::vector<double> outputTriang(maxIndex);
  std::cout << "line2.2\n";
  std::vector<double> outputTriang0(maxIndex);
  std::cout << "line2.3\n";
  std::vector<double> outputTriang1(maxIndex);
  std::cout << "line2.4\n";
  std::vector<double> outputTriang2(maxIndex);
  std::cout << "line2.5\n";
  std::vector<double> outputTriang3(maxIndex);
  std::cout << "line3\n";

  std::ofstream iss("../res/basicHistogram2.csv");

  MatrixCollection<4> m = MatrixCollection<4>::readFromFile(PATH_4_ALL);
  std::cout << "Read finished! The number of matrices is " << m.size() <<"\n";
  m.regularize ();

  std::vector<double> ratiosAll;
  std::vector<Ush> numberTriang;
  for (size_t i = 0; i < m.size(); i++) {
    if (i % 10070 == 0) {
        std::cout << std::setprecision(2) << ((double(i)) / m.size() * 100) << "%\n";
    }
    ratiosAll.push_back(m[i].getConsistencyRatio());
    numberTriang.push_back(6 - m[i].countIndexOfUpperTriangle ());
  }

  std::cout << "Consistencies calulcated.\n";
  std::cout << "Creating histogram...\n";

  std::cout << "line1\n";

  if (ratiosAll.size() != numberTriang.size())
    throw "Sizes mismtach";

  std::cout << ratiosAll.size() << "\n";
  for (size_t i = 0; i < ratiosAll.size(); i++) {
    std::cout << i << "\n";
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
                std::cout << "OH MY GOD! Something really went wrong.\n";
                throw "ERROR!";
            break;
        }
      }
    }
  }

  std::cout << "line2\n";

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

  std::cout << "line4\n";

  std::cout << "Bucketing done, writing file...\n";

  std::cout << "Program ok\n";
  //Header:
  //I << "interval\t#all\t%eigen\t%spantree\t%cosine\t%common\t\n";
  for (size_t i = 0; i < bucketsAll.size(); i++) {
    std::cout << i << " " << bucketsAll.size() << "\n";
    iss << i * steps << "\t" << outputAll.at(i) << "\t" <<
        outputTriang.at(i) << "\t" << outputTriang0.at(i) << "\t" << outputTriang1.at(i) << "\t" << outputTriang2.at(i) << "\t" << outputTriang3.at(i) << "\n";
  }
  std::cout << "About to close the connection.\n";
  iss.close();
  std::cout << "Procedure finished.\n";

  return 0;
}
