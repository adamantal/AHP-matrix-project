#include "GraphSolution.h"

bool MatrixNode::isGroupped()const{
  return (group != 0);
}
Ush MatrixNode::getGroup()const{
  return group;
}
void MatrixNode::groupIt(Ush g) {
  if (group != 0)
    throw "Node is already groupped!\n";
  else
    group = g;
}
void MatrixNode::ungroupIt(){
  if (group == 0)
    throw "Node is not groupped.\n";
  else
    group = 0;
}
MatrixNode& MatrixNode::operator*= (const double &rhs){
  X *= rhs;
  return *this;
}


template<size_t N>
MatrixProcessor<N>::MatrixProcessor(Matrix<N> m):m(m) {
  for (size_t i = 0; i < N; i++){
    //MatrixNode* n = new MatrixNode(i);
    nodes[i] = new MatrixNode(i);
  }
}
template<size_t N>
MatrixProcessor<N>::~MatrixProcessor() {
  for (size_t i = 0; i < N; i++){
    delete nodes[i];
  }
}


template<size_t N>
std::vector<double> MatrixProcessor<N>::run(EdgeList e) {
  if (e.size() + 1 != N) throw "Invalid edgelist given.\n";
  for (auto it = e.begin(); it != e.end(); it++) {
    Ush index = std::distance(e.begin(), it) + 1;
    if (it->first >= N || it->second >= N) throw "Invalid index in the edgelist.\n";
    if (it->first == it->second) throw "Indexes match, can't go on.\n";
    double elem = m.get(it->second, it->first);

    MatrixNode* n1 = nodes[it->first];
    MatrixNode* n2 = nodes[it->second];

    if (!n1->isGroupped() && !n2->isGroupped()) {
      //std::cout << "C1" << index << std::endl;
      n1->groupIt(index);
      n2->groupIt(index);

      *n2 *= elem;

    } else if (!(n1->isGroupped() && n2->isGroupped())) {
      //std::cout << "C2" << index << std::endl;
      //assert n1 groupped and n2 not:
      if ((!n1->isGroupped()) && n2->isGroupped()) {
        //std::swap(n1, n2);
        MatrixNode* tmp = n1;
        n1 = n2;
        n2 = tmp;
      }
      n2->groupIt(n1->getGroup());
      *n2 *= (n1->getX() * elem);

    } else if (n1->isGroupped() && n2->isGroupped()) {
      //std::cout << "C3" << index << std::endl;
      elem *= n1->getX();
      Ush oldGroup = n2->getGroup();

      for (size_t i = 0; i < nodes.size(); i++) {
        if (nodes[i]->getGroup() == oldGroup) {
          *(nodes[i]) *= elem;
          nodes[i]->groupIt(n1->getGroup());
        }
      }
    } else {
      throw "Invalid ifelse branch.\n";
    }
  }
  //checking if solution is valid:
  auto it = nodes.begin();
  Ush allGroup = it->second->getGroup();
  //std::cout << allGroup << std::endl;
  it++;

  for (; it != nodes.end(); it++) {
    if (it->second->getGroup() != allGroup) {
      std::cout << it->second->getGroup() << std::endl;
      throw "The edgeslist not spanned all the nodes.\n";
    }
  }

  std::vector<double> r;
  for (size_t i = 0; i < N; i++) {
    r.push_back(nodes[i]->getX());
  }
  return r;
}

bool plus(std::vector<Ush> &v) {
  for (auto rit = v.rbegin(); rit != v.rend(); rit++) {
    if (*rit != 7) {
      (*rit)++;
      return true;
    } else {
      *rit = 0;
    }
  }
  return false;
}

int main() {
  //Matrix<5> m({1, 3, 9, 3, 10, 11, 14, 0, 6, 3});
  Matrix<4> m({1, 3, 7, 1, 3, 1});
  std::cout << m << std::endl << std::endl;
  std::cout << m.getConsistencyRatio() << std::endl;
  unsigned long int bad = 0, good = 0;

  std::vector<Ush> v;
  for (Ush i = 0; i < 6; i++) {
    v.push_back(0);
  }
  while (plus(v)) {
    MatrixProcessor<4> mp = MatrixProcessor<4>(m);
    std::vector<double> r;
    try {
      std::pair<Ush, Ush> p1(v[0],v[1]);
      std::pair<Ush, Ush> p2(v[2],v[3]);
      std::pair<Ush, Ush> p3(v[4],v[5]);
      //std::pair<Ush, Ush> p4(v[6],v[7]);
      EdgeList e;
      e.push_back(p1);
      e.push_back(p2);
      e.push_back(p3);
      //e.push_back(p4);

      r = mp.run(e);
      Matrix<4>::L1(r);
      for (size_t i = 0; i < r.size(); i++) {
        std::cout << r[i] << " ";
      }
      std::cout << " /////// " << ((m.testVectorParetoOptimal(r))?"OK":"BAD") << std::endl;
      if (!(m.testVectorParetoOptimal(r))) {
        bad++;
      } else {
        good++;
      }
      for (size_t i = 0; i < v.size(); i++) {
          std::cout << v[i] << " ";
        }
        std::cout << std::endl;
      
    } catch (const char* e){
      //std::cout << e;
    }
  }
  std::cout << "Good: " << good << std::endl;
  std::cout << "Bad: " << bad << std::endl;


  /*try {
    r = mp.run(e);
  } catch(const char * c){
    std::cout << c << std::endl;
    throw "";
  }*/
  return 0;
}
