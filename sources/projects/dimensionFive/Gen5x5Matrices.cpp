#include <iostream>
#include <string>
#include <algorithm>

#include "MatrixStr.h"
#include "Matrix.hpp"

typedef unsigned long long int Bnum;
const Bnum SUM5 = 2015993900449;

const std::vector<double> elements = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,1.0/2,1.0/3,1.0/4,1.0/5,1.0/6,1.0/7,1.0/8,1.0/9};

template<size_t N>
std::vector<Ush> indexToVector (Bnum x) {
    std::vector<Ush> v;
    for (Ush i = 0; i < N * (N - 1) / 2; i++) {
        Ush tmp = x % elements.size();
        x /= elements.size();
        v.push_back(tmp);
    }
    std::reverse(v.begin(), v.end());
    if (v.size() != 10)
        throw "The size does not match in indexToVector function!\n";
    return v;
}

Bnum vectorToIndex (std::vector<Ush> v) {
    if (v.size() != 10)
        throw "The size does not match in vectorToIndex function!\n";
    Bnum ret;
    for (Ush i = 0; i < v.size(); i++) {
        ret *= elements.size();
        ret += v.at(i);
    }
    return ret;
}

int main() {
    //Generating all the permutations:
    std::vector<MatrixStr<5>> perms;
    {
        MatrixStr<5> ms;
        Ush perm[] = {0, 1, 2, 3, 4};
        do {
            MatrixStr<5> tmp = ms.permutateBy(perm);
            perms.push_back(tmp);
        } while (std::next_permutation (perm, perm + 5));
    }

    //Start vector:
    std::vector<Ush> vec;
    for (Ush i = 0; i < 10; i++) {
        vec.push_back(0);
    }

    Bnum slices = 100;
    Bnum sliceIndex = SUM5 / slices;

    //Input from user -- which file to be generated
    short which;
    std::cout << "Which part of the file do you want to generate (/" << std::to_string(slices) << "): ";
    std::cin >> which;

    if (which != 0) return 0;
    else {
        //Bnum indexDelta = which * SUM5;
        //Bnum indexMax = SUM5 / slices;

        std::vector<bool>* filter = nullptr;

        filter = new std::vector<bool>(sliceIndex);

        Bnum percent0_01 = SUM5 / 10000;
        Bnum i = 0;
        clock_t cycleBegin = clock(), cycleEnd;
        for (auto it = filter->begin(); it != filter->end(); it++) {
            if (i % percent0_01 == 0) {
                cycleEnd = clock();
                double elapsed_secs = double(cycleEnd - cycleBegin) / CLOCKS_PER_SEC;
                cycleBegin = cycleEnd;
                std::cout << elapsed_secs << " seconds while progress increased by 0.01%\n";
                std::cout << (double)(std::distance (filter->begin(), it) / percent0_01) * 0.01 << "%\n";
            }
            if (!*it) { //has not yet been filtered
                std::vector<Ush> mVector = indexToVector<5>(i);
                for (auto it = perms.begin(); it != perms.end(); it++) {
                    std::vector<Ush> mPermed = it->applyTo(mVector);
                    try {
                        filter->at(vectorToIndex(mPermed)) = true;
                    } catch (...) {}
                }
            }
            i++;
        }
        std::cout << std::count(filter->begin(), filter->end(), true) << " out of " << sliceIndex << std::endl;
        std::cout << "success" << std::endl;
    }

    /*std::vector<bool>* filter = nullptr;
    bool notSucceeded = true;
    Bnum N = 30; //1e4 OK
    while (notSucceeded) {
        try {
            clock_t begin,end;
            double elapsed_secs;
            begin = clock();

            Bnum fin = SUM5 / N;
            filter = new std::vector<bool>(fin);
            std::cout << "Allocated!\n";
            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            std::cout << "Has taken " << elapsed_secs << " seconds.\n";

            begin = clock();
            end = clock();
            for (Bnum index = 0; index < SUM5 / N; index++) {
                if (index % (fin / 100) == 0) {
                    end = clock();
                    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    std::cout << "Index: " << index << " - percent:" << index / fin << std::endl <<
                                 "Secs from loop start: " << elapsed_secs << std::endl;
                }

                filter->at(index) = false;
            }
            std::cout << "Indexes set!\n";
            notSucceeded = false;
        } catch (std::bad_alloc) {
            std::cout << N++ << std::endl;
            notSucceeded = true;
        }
    }
    std::cout << "success:" << N << std::endl;*/

    return 0;
}
