#ifndef LPSOLUTION_H
#define LPSOLUTION_H

#include<vector>
#include "Matrix.h"
#include<set>
#include<algorithm>

template<size_t N>
class LpSolution {
  private:
    //which does not requires the LP:
    Matrix<N> matrix;
    std::vector<double> initv;
    std::vector<spair> I,J;

    //the output of the LP:
    double value;
    std::vector<double> x;
    std::vector<double> s;
  public:
    LpSolution(Matrix<N> m, std::vector<double> v, std::vector<spair> i,  std::vector<spair> j, double val, std::vector<double> xv, std::vector<double> sv)
      :matrix(m), initv(v), I(i), J(j), value(val), x(xv), s(sv)  {}
    bool isOptimal() {
      return (value > -1e-8);
    }

    void printOutData()const{
      std::cout << "Matrix: \n";
      std::cout << matrix << std::endl;

      std::cout << "Vector: \n";
      for (size_t i = 0; i < N; i++) std::cout << initv[i] << " ";
      std::cout << std::endl << std::endl;

      std::cout << "Elements of I: ";
      for (size_t i = 0; i < I.size(); i++) std::cout << "(" << I[i].a << ", " << I[i].b << ") ";
      std::cout << std::endl;
      std::cout << "Elements of J: ";
      for (size_t i = 0; i < J.size(); i++) std::cout << "(" << J[i].a << ", " << J[i].b << ") ";
      std::cout << std::endl << std::endl;

      std::cout << "Value: \n";
      std::cout << value << std::endl;

      std::cout << "x: \n";
      for (size_t i = 0; i < N; i++) std::cout << x[i] << " ";
      std::cout << std::endl;

      std::cout << "s: \n";
      for (size_t i = 0; i < s.size(); i++) std::cout << s[i] << " ";
      std::cout << std::endl;
    }

    Matrix<N> getMatrix()const{ return matrix; }
    std::vector<double> getw(){ return initv; }
    std::vector<double> getx(){ return x; }
    std::vector<double> getxnorm(){
      std::vector<double> tmp = x;
      double sum = 0;
      for (size_t i = 0; i < N; i++) sum += x[i];
      if (abs(sum) > 1e-6) {
        for (size_t i = 0; i < N; i++) tmp[i] /= sum;
      }
      return tmp;
    }
    std::vector<double> gets(){ return s; }
    std::vector<double> getOtherTwoVector(){
      //determining the index:
      std::set<int> indexes;
      for (size_t i = 0; i < I.size(); i++) {
        if (s[i] < 1e-8) {
          indexes.insert(I[i].a);
          indexes.insert(I[i].b);
        }
      }
      int specialIndex = -1;
      for (size_t i = 0; i < N; i++){
        if (std::find(indexes.begin(), indexes.end(), i) == indexes.end()) {
          specialIndex = i;
          break;
        }
      }
      if (specialIndex == -1) throw "Fundamental problem in calculation - specialIndex can not be determined.";
      //leírni: 4x4-esre csak ez a lehetőség fordulhat elő
      //std::cout << "Ind: " << specialIndex << "\n";

      std::vector<std::vector<double>> returns;
      for (size_t i = 0; i < N; i++) {
        if (i != specialIndex) {
          std::vector<double> tmp = x;
          tmp[specialIndex] = matrix.get(specialIndex, i) * x[i];
          Matrix<N>::L1(tmp);
          returns.push_back(tmp);
        }
      }
      //std::cout << "Its ok until now" << std::endl;

      std::vector<double> xnorm = getxnorm();
      for (auto it = begin(returns); it != end(returns); it++) {
        bool equal = true;
        for (size_t i = 0; i < N; i++) {
          if (abs((*it)[i] - xnorm[i]) > 1e-6) {
            equal = false;
            break;
          }
        }
        if (equal) {
          returns.erase(it);
          break;
        }
      }

      if (returns.size() != 2) throw "Invalid amouont of vectors - assumption not hold.";
      //mi van ha kettő megegyezik? - ezen szűk mátrixlistákban nem fordul elő - mi van a nagyobb halmazzal?

      if (returns[0][0] < returns[1][0]) {
        for (size_t i = 0; i < N; i++) returns[0].push_back(returns[1][i]);
        return returns[0];
      } else {
        for (size_t i = 0; i < N; i++) returns[1].push_back(returns[0][i]);
        return returns[1];
      }
    }
};

#endif //LPSOLUTION_H
