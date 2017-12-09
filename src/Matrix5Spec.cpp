#include <iostream>
#include <algorithm>
#include "Matrix.h"

typedef unsigned short Ush;

class Matrix5Spec : public Matrix<4> {
  private:
    std::vector<Ush> predata;
  public:
    //Warning to the user: please make sure which constructor you're using!!
    Matrix5Spec (const std::vector<Ush> &v):Matrix(v) {
      for (Ush i = 0; i < 4; i++) {
        predata.push_back(0);
      }
    }
    Matrix5Spec (const std::vector<Ush>&pre, const std::vector<Ush> &v):Matrix(v),predata(pre) {
    }

    static const Matrix5Spec getMatrixOutOfIndex(const unsigned long long int& i, MatrixCollection<4> &mc) {
      unsigned long int preindex = i % 83521; // 1007096
      std::vector<Ush> pre;
      for (Ush i = 0; i < 4; i++) {
        pre.push_back(preindex % 17);
        preindex /= 17;
      }
      std::reverse(pre.begin(), pre.end()); // is this really necessary? YES - otherwise it throws error!

      std::vector<Ush> post = mc[i / 83521].getData();
      return (Matrix5Spec(pre, post));
    }

    static bool getIndexOutOfMatrix(const Matrix<5>& m, unsigned long long int &i, const MatrixCollection<4> &mc) {
      unsigned long long int j = 0;
      if (!mc.isIncluded(m.cutBottom(), j))
        return false;
      else {
        std::vector<Ush> data = m.getData();
        unsigned long long int x = data[0];
      	for (size_t i = 1; i < 4; i++) {
      		x = Matrix<0>::elem.size() * x + data[i];
      	}
        i = x + 1007096 * j;
        return true;
      }
    }

    Matrix<5> getMatrix()const{
      std::vector<Ush> r = predata;
      r.insert(r.end(), data.begin(), data.end());
      return Matrix<5>(r);
    }
};
