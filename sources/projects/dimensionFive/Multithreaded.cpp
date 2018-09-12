// I/O streams
#include <iostream>
#include <fstream>

// Datastructures
#include <vector>
#include <list>

// Algorithm utilities
#include <algorithm>

// Multithreading support
#include <thread>
#include <mutex>

// Matrix permutation calculator utility class
#include "MatrixStr.h"
#include "File.hpp"

// typedefs
typedef unsigned long long int Ulli;

// forward declarations
class ThreadController;
void excludePermutations(ThreadController* controller, Ulli index, std::vector<MatrixStr<5>>* perms);

// constants
const Ulli ALL = 2015993900449;
const std::vector<double> elements = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,1.0/2,1.0/3,1.0/4,1.0/5,1.0/6,1.0/7,1.0/8,1.0/9};

// globals
std::vector<bool>* filter;

template<size_t N>
std::vector<unsigned short> indexToVector (Ulli x) {
    std::vector<unsigned short> v;
    for (unsigned short i = 0; i < N * (N - 1) / 2; i++) {
        unsigned short tmp = x % elements.size();
        x /= elements.size();
        v.push_back(tmp);
    }
    std::reverse(v.begin(), v.end());
    if (v.size() != 10)
        throw "The size does not match in indexToVector function!\n";
    return v;
}

Ulli vectorToIndex (std::vector<unsigned short> v) {
    if (v.size() != 10)
        throw "The size does not match in vectorToIndex function!\n";
    Ulli ret = 0;
    for (Ush i = 0; i < v.size(); i++) {
        ret *= elements.size();
        ret += v.at(i);
    }
    return ret;
}

inline void findProperNext(Ulli& next) {
    while (filter->at(next)) {
        std::cout << "next\n" << next << std::endl;
        next++;
    }
}

class ThreadController {
private:
    Ulli finished;
    Ulli next;
    Ulli maxIndex;

    unsigned int maxThreads;
    unsigned int threads;

    std::mutex mut;

    std::vector<MatrixStr<5>> perms;
    Ulli percent0_01;

    clock_t cycleBegin, cycleEnd;

    std::ofstream* file;

private:
    void startThread() {
        findProperNext(next);
        *file << next << ",";
        std::thread thread(excludePermutations, this, next++, &perms);
        thread.detach();
        threads++;
    }

public:
    ThreadController(Ulli maxInd = ALL / 40):
        finished(0),
        next(0),
        maxIndex(maxInd),
        maxThreads(std::thread::hardware_concurrency() * 2), // TODO to be tested
        threads(0)
    {
        std::cout << "ThreadController setting up.\n";

        filter = new std::vector<bool>(maxInd); //filled with falses
        percent0_01 = maxInd / 10000000;
        std::cout << "Vector of size " << filter->size() << " initialized!\n";

        {
            MatrixStr<5> ms;
            Ush perm[] = {0, 1, 2, 3, 4};
            do {
                MatrixStr<5> tmp = ms.permutateBy(perm);
                perms.push_back(tmp);
            } while (std::next_permutation (perm, perm + 5));
        }
        std::cout << "Permutations calculated.\n";

        std::vector<Ush> vec;
        for (Ush i = 0; i < 10; i++) {
            vec.push_back(0);
        }

        cycleBegin = clock();
        file = new std::ofstream("../results/minimalPermuated.txt");
    }
    void start() {
        {
            std::lock_guard<std::mutex> localLock(mut);
            for (unsigned int i = 0; i < maxThreads; i++) {
                startThread();
            }
        }
        while (true) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            {
                std::lock_guard<std::mutex> localLock(mut);
                *file << "\n";
                if (finished == maxIndex) {
                    break;
                }
            }
        }
    }
    void threadExited() {
        std::lock_guard<std::mutex> localLock(mut);
        if (next % percent0_01 == 0) {
            cycleEnd = clock();
            double elapsed_secs = double(cycleEnd - cycleBegin) / CLOCKS_PER_SEC;
            cycleBegin = cycleEnd;
            std::cout << elapsed_secs << " seconds while progress increased by 0.01% -- " << next / percent0_01 * 0.00001 << "%\n";
        }
        finished++;
        if (--threads >= maxThreads) {
            std::cout << "ERROR! More thread than allowed!\n";
            std::terminate();
        }
        if (next < maxIndex) {
            startThread();
        }
    }
};

void excludePermutations(ThreadController* controller, Ulli index, std::vector<MatrixStr<5>>* perms) {
    std::vector<Ush> mVector = indexToVector<5>(index);
    for (auto it = perms->begin(); it != perms->end(); it++) {
        std::vector<Ush> mPermed = it->applyTo(mVector);
        try {
            filter->at(vectorToIndex(mPermed)) = true;
        } catch (...) {}
    }
    controller->threadExited();
}

int main() {
    std::cout << "The number of threads that can run concurrently is "
        << std::thread::hardware_concurrency() << std::endl;

    ThreadController controller;
    controller.start();
    std::cout << "Procedure successfully completed!\n";
    return 0;
}
