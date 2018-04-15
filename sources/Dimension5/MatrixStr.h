#ifndef MATRIXS_H
#define MATRIXS_H

typedef unsigned short Ush;

#include <stdlib.h>
#include <vector>
#include <string>

template<size_t N>
class MatrixStr {
private:
    std::vector<short> inner;

    Ush indexOfInverse(Ush i) {
    	if (i == 0)
    		return 0;
    	if (i > 8)
    		return i - 8;
    	else
    		return i + 8;
    }

public:
    MatrixStr() {
        std::vector<short> letters = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        for (Ush i = 0; i < N * (N - 1) / 2; i++)
            inner.push_back(letters.at(i));
    }
    MatrixStr(std::vector<short> l):inner(l) {}

    short indexOfElement(Ush i, Ush j) const {
    	if (i == j)
    		return 0; //should probably never happen
    	else if (i < j)
    		return inner[j - 1 + (2 * N - i - 3) * i / 2];
    	else {
    		return (-1) * indexOfElement(j, i);
    	}
    }
    MatrixStr<N> permutateBy(Ush p[]) const {
        std::vector<short> v;
    	for (size_t i = 0; i < N - 1; i++) {
    		for (size_t j = i + 1; j < N; j++) {
    			v.push_back(indexOfElement(p[i], p[j]));
    		}
    	}
    	return MatrixStr<N>(v);
    }
    std::vector<short> get() const{
        return inner;
    }
    std::string toString() const {
        std::string ret;
        for (Ush i = 0; i < N * (N - 1) / 2; i++)
            ret += std::to_string(inner.at(i)) + ",";
        return ret;
    }
    std::vector<Ush> applyTo (std::vector<Ush> matrixElems) {
        if (matrixElems.size() != N * (N - 1) / 2)
            throw "Oh noo!\n";

        std::vector<Ush> ret;
        for (Ush i = 0; i < N * (N - 1) / 2; i++) {
            if (inner[i] > 0)
                ret.push_back (matrixElems[inner[i] - 1]);
            else if (inner[i] < 0)
                ret.push_back (indexOfInverse (matrixElems[-1 * inner[i] - 1]));
            else
                throw "Can not be zero!\n";
        }
        return ret;
    }
};

#endif //MATRIXS_H
