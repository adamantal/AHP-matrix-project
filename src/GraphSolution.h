#include<iostream>
#include<vector>
#include "Matrix.h"

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
    Ush group; //zero means non groupped

    double X;
  public:
    MatrixNode(Ush l):label(l),group(0),x(-1.0){}
    Ush getLabel()const{ return label}
    bool isGroupped()const{return (group != 0);}
    void setGroup(Ush g){group = g;}
    Ush getGroup()const{return group;}

    void setX(double x){X = x;}
    double getX()const{}
};

template<size_t N>
class AlgorithmProcessor {
  private:
    std::map<Ush, MatrixNode> nodeMap;
    Matrix<N> m;

    //utility:
    Ush nextGroupLabel = 1;

    void groupNodes(Ush label1, Ush label2) {
      MatrixNode * v1 = & nodeMap[label1];
      MatrixNode * v2 = & nodeMap[label2];

      if ((!v1 -> isGroupped()) && (!v2 -> isGroupped())) {
        v1 -> setGroup(nextGroupLabel);
        v2 -> setGroup(nextGroupLabel);
        nextGroupLabel++;

        //HERE GOES THE CALCAULATION
      } else {
        if (v1 -> isGroupped() && v2 -> isGroupped()) {
          Ush minGroup = min(v1 -> getGroup(), v2 -> getGroup());

          v1 -> setGroup(minGroup);
          v2 -> setGroup(minGroup);

          //HERE GOES ANOTHER CALCAULATION
        } else if (v1 -> isGroupped()) {
          v2 -> setGroup(v1 -> getGroup());

          //HERE GOES ANOTHER CALCAULATION
        } else if (v2 -> isGroupped()) {
          v1 -> setGroup(v2 -> getGroup());

          //HERE GOES ANOTHER CALCAULATION
        }
      }
    }
  public:
    template<size_t N>
    AlgorithmProcessor<N>(Matrix<N> n):m(n){}

    std::vector<double> run();
    std::vector<double> run(EdgeList e) {
      for (auto it = e.begin(); it != e.end(); it++) {
        groupNodes(it -> first(), it -> second());
      }
      std::vector<double> v;
      for (auto it = nodeMap.begin(); it != nodeMap.end(); it++) {
        v.push_back(it -> second -> getX());
      }
      return v;
    }
};

#endif //GRAPHSOLUTION_H
