#include <iostream>
#include "Dense"
#include <vector>
#include <algorithm>

using Eigen::MatrixXd;
using namespace std;

const vector<double> e = {1,2,3,4,5,6,7,8,9,1/2,1/3,1/4,1/5,1/6,1/7,1/8,1/9};

const int elemsize = 17;

MatrixXd permutateOneMatrix(const MatrixXd &M, int p[]) {
	MatrixXd mp(4,4);
	mp = M;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			mp(i,j) = M(p[i],p[j]);
		}
	}
	return mp;
}

void takeOutPermutations(vector<MatrixXd> &vec, long int &i) {
	int cols[] = {0,1,2,3};
	
  do {
	  bool test = false;
		for (long int k = 0; k < i; k++) {
			if (permutateOneMatrix(vec[i],cols) == vec[k]) {
				i--;
				test = true;
				break;
//				cout << "permutation found" <<endl;
			}
		}
		if (test) {
			break;
		}
  } while ( std::next_permutation(cols,cols+4) );
}

void getAllMatrices(vector<MatrixXd> &vec) {
	long int i = 0;
  for (int i1 = 0; i1 < elemsize; i1++) {
  	for (int i2 = 0; i2 < elemsize; i2++) {
	  	std::cout << "Loading datas: " << i1*17+i2 << "/" << elemsize*elemsize << std::endl;
 	  	for (int i3 = 0; i3 < 17; i3++) {
 		  	for (int i4 = 0; i4 < 17; i4++) {
 	 		  	for (int i5 = 0; i5 < 17; i5++) {
 	 	 		  	for (int i6 = 0; i6 < 17; i6++) {
		 		  		MatrixXd tmp(4,4);
		 		  		tmp << 1,e[i1],e[i2],e[i3],1/e[i1],1,e[i4],e[i5],1/e[i2],1/e[i4],1,e[i6],1/e[i3],1/e[i5],1/e[i6],1;
		 		  		vec[i] = tmp;
		 		  		takeOutPermutations(vec,i);
		 		  		i++;
	 		  		}
 		  		}
 		  	}
  		}
	  }
  }
}

int main() {
  std::vector<MatrixXd> vec(24137569);
  getAllMatrices(vec);
  std::cout << vec.size() << std::endl; 
/*  for (long int i = 0; i < vec.size(); i++) {
  	if (i % 1000000 == 0) {
  		std::cout << i << std::endl;
  	}
  	vec[i] = vec[i].inverse();
  }*/
  return 0;
}

