#ifndef GRAPHSOLUTION_HPP
#define GRAPHSOLUTION_HPP

#include <iostream>
#include <vector>
#include <map>
#include "Matrix.hpp"

typedef unsigned short Ush;
typedef std::vector<std::pair<Ush, Ush>> EdgeList;

/*
inputs:
order of edges, we'll work on them by "contracting" one by one of those nodes
the matrix
*/

class MatrixNode {
  private:
    const Ush label;
    double X;
    Ush group;

  public:
    MatrixNode(Ush l):label(l),X(1.0),group(0){}

    bool isGroupped()const;
    Ush getGroup()const;
    void groupIt(Ush g);
    void ungroupIt();

    Ush getLabel()const{ return label;}
    //note that you can not set label after constructed

    void setX(double x){X = x; }
    double getX()const{ return X;}

    MatrixNode& operator*= (const double &rhs);
    MatrixNode& operator/= (const double &rhs);
};

template<size_t N>
class MatrixProcessor {
  private:
    std::map<Ush, MatrixNode*> nodes;
    Matrix<N> m;

    static const std::vector<std::pair<Ush, Ush>> edges;

    std::vector<double> runForEdgeList(EdgeList);
    bool plus(std::vector<Ush> &, EdgeList &) const;
  public:
    MatrixProcessor<N>(Matrix<N> m);
    ~MatrixProcessor<N>();
    Ush run();
    void printNodesValue();
};

class VectorEps {
  private:
    std::vector<std::vector<double>> data;

    static const double EPS;
    static bool isEpsClose(double x, double y) {
      return fabs(x - y) < EPS;
    }
    static bool isEpsClose(std::vector<double> v1, std::vector<double> v2) {
      if (v1.size () != v2.size()) throw "Size of vectors mismtach!\n";
      for (size_t i = 0; i < v1.size(); i++) {
        if (!isEpsClose(v1[i], v2[i])) {
          return false;
        }
      }
      return true;
    }
  public:
    VectorEps(){}
    void add(std::vector<double> x) {
      bool foundVectorEps = false;
      for (auto it = data.begin(); it != data.end(); it++) {
        if (isEpsClose(*it, x)) foundVectorEps = true;
      }
      if (!foundVectorEps) data.push_back(x);
    }
    std::vector<std::vector<double>> getVectors()const{
      return data;
    }
    size_t size() const;
};

#endif //GRAPHSOLUTION_HPP
