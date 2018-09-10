#include "GraphSolution.hpp"

const double VectorEps::EPS = 1e-8;

bool MatrixNode::isGroupped()const{
  return (group != 0);
}
Ush MatrixNode::getGroup()const{
  return group;
}
void MatrixNode::groupIt(Ush g) {
  if (isGroupped())
    throw "Node is already groupped!\n";
  else
    group = g;
}
void MatrixNode::ungroupIt() {
  if (group == 0)
    throw "Node is not groupped.\n";
  else
    group = 0;
}
MatrixNode& MatrixNode::operator*= (const double &rhs){
  X *= rhs;
  return *this;
}
MatrixNode& MatrixNode::operator/= (const double &rhs){
  X /= rhs;
  return *this;
}

template<size_t N>
MatrixProcessor<N>::MatrixProcessor(Matrix<N> m):m(m) {
  for (size_t i = 0; i < N; i++) {
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
void MatrixProcessor<N>::printNodesValue() {
  for (size_t i = 0; i < N; i++) {
    std::cout << nodes[i]->getX() << " ";
  }
  std::cout << std::endl;
}

template<size_t N>
std::vector<double> MatrixProcessor<N>::runForEdgeList(EdgeList e) {
  if (e.size() + 1 != N) throw "Invalid edgelist given.\n";
  //printNodesValue();
  for (auto it = e.begin(); it != e.end(); it++) {
    Ush index = std::distance(e.begin(), it) + 1;
    if (it->first >= N || it->second >= N) throw "Invalid index in the edgelist.\n";
    if (it->first == it->second) throw "Indexes match, can't go on.\n";

    MatrixNode* n1 = nodes[it->first];
    MatrixNode* n2 = nodes[it->second];
    //std::cout << "Groupping " << n1->getLabel() << " and " << n2->getLabel() << "\n";

    if (!n1->isGroupped() && !n2->isGroupped()) {
      //std::cout << "\tNone is groupped.\n";
      n1->groupIt(index);
      n2->groupIt(index);

      *n2 *= m.get(it->second, it->first);
      //std::cout << "\t The second node's value is " << m.get(it->second, it->first) <<"\n";

    } else if (!(n1->isGroupped() && n2->isGroupped())) {
      //assert n1 is groupped and n2 is not:
      if (n2->isGroupped() && !n1->isGroupped()) {
        //std::swap(n1, n2);
        MatrixNode* tmp = n1;
        n1 = n2;
        n2 = tmp;
      }
      //std::cout << "\tOne is groupped.\n";
      double minX = 100001;
      for (auto nodeit = nodes.begin(); nodeit != nodes.end(); nodeit++) {
        if (nodeit->second->getGroup() == n1->getGroup()) {
          //std::cout << "\t\tGroup member: " << nodeit->second->getLabel() << " with value " << m.get(n2->getLabel(), nodeit->first) << "\n";
          minX = std::min(minX, nodeit->second->getX() * m.get(n2->getLabel(), nodeit->first) ); //nodeit->second->getX());
        }
      }
      //std::cout << "\tThe min is " << minX << " assigning it to " << n2->getLabel() << "\n";
      if (minX > 100000) throw "Calculation error - no entry found in group.\n";

      *n2 *= minX;
      n2->groupIt(n1->getGroup());

    } else if (n1->isGroupped() && n2->isGroupped()) {
      if (n1->getGroup() == n2->getGroup()) throw "Circle attempted.\n";
      //we're going to merge n2's group to n1's
      double minX = 100001;
      for (auto n2it = nodes.begin(); n2it != nodes.end(); n2it++) {
        MatrixNode * elementOfn2 = n2it->second;
        if (elementOfn2->getGroup() == n2->getGroup()) {
          for (auto n1it = nodes.begin(); n1it != nodes.end(); n1it++) {
            MatrixNode* elementOfn1 = n1it->second;
            if (elementOfn1->getGroup() == n1->getGroup()) {
              minX = std::min(minX, elementOfn1->getX() * m.get(elementOfn2->getLabel(), elementOfn1->getLabel()));
            }
          }
        }
      }
      if (minX > 100000) throw "Calculation error - no entry found in group.\n";

      for (auto n2it = nodes.begin(); n2it != nodes.end(); n2it++) {
        MatrixNode * elementOfn2 = n2it->second;
        if (elementOfn2->getGroup() == n2->getGroup()) {
          *elementOfn2 *= minX;
          elementOfn2->groupIt(n1->getGroup());
        }
      }
    } else {
      throw "Invalid ifelse branch.\n";
    }
    //printNodesValue();
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

template<>
const std::vector<std::pair<Ush, Ush>> MatrixProcessor<4>::edges = {std::make_pair<Ush, Ush>(0,1),
  std::make_pair<Ush, Ush>(0,2), std::make_pair<Ush, Ush>(0,3), std::make_pair<Ush, Ush>(1,2),
  std::make_pair<Ush, Ush>(1,3), std::make_pair<Ush, Ush>(2,3)};

template<>
bool MatrixProcessor<4>::plus(std::vector<Ush> &indexes, EdgeList &e) const {
  /*for (size_t i = 0; i < indexes.size(); i++) {
    std::cout << indexes[i] << " ";
  }
  std::cout << "; ";*/

  size_t i = 3;
  do {
    i--;
    if (indexes[i] != 5) {
      indexes[i]++;
      break;
    } else {
      indexes[i] = 0;
    }
  } while (i != 0);

  if (i == 0 && indexes[0] == 0) return false;

  //if the tree is not spanning:
  if ( ((indexes[0] != 2 && indexes[0] != 4 && indexes[0] != 5) &&
       (indexes[1] != 2 && indexes[1] != 4 && indexes[1] != 5)) || //3
       (indexes[0] % 2 != 1 && indexes[1] % 2 != 1) || //2
       ((indexes[0] != 0 && indexes[0] != 3 && indexes[0] != 4) &&
       (indexes[1] != 0 && indexes[1] != 3 && indexes[1] != 4)) || //1
       indexes[0] == indexes[1] || indexes[1] == indexes[2] || indexes[0] == indexes[2]
     )
  {
    return plus(indexes, e);
  }

  size_t ri = 3;
  do {
    ri--;
    e[ri] = edges[indexes[ri]];
  } while (ri != 0);

  return true;
}

template<>
Ush MatrixProcessor<4>::run() {
  VectorEps vectors;

  std::vector<Ush> indexes(3);
  EdgeList e(3);

  while (plus(indexes, e)) {
    MatrixProcessor<4> mp = MatrixProcessor<4>(m);
    std::vector<double> r;
    try {
      r = mp.runForEdgeList(e);
      Matrix<4>::L1(r);
      if (!(m.testVectorParetoOptimal(r))) {
        throw "Error in calulcation - the received vector is not optimal.\n";
        break;
      }
      vectors.add(r);
    } catch (const char* e) {
      //std::cout << e;
    }
  }
  return vectors.getVectors().size();

  /*std::cout << "Good vectors:\n";
  for (size_t i = 0; i < results.size(); i++) {
    std::cout << "(";
    for (size_t j = 0; j < results[i].size() - 1; j++) {
      std::cout << results[i][j] << ", ";
    }
    std::cout << results[i][results[i].size() - 1] << ")\n";
  }*/
}
