#include<iostream>
#include<vector>
#include<map>

#include <lemon/lp.h>

#include "Matrix.cpp"

#ifndef GRAPHSOLUTION_H
#define GRAPHSOLUTION_H

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
  public:
    MatrixProcessor<N>(Matrix<N> m);
    ~MatrixProcessor<N>();
    std::vector<double> run(EdgeList);
};

#endif //GRAPHSOLUTION_H
