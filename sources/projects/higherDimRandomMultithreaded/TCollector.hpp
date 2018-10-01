#ifndef COLLECTOR_HPP
#define COLLECTOR_HPP

#include <list>

#include "TRoutine.hpp"

template<size_t N>
class Collector {
private:
    std::list<RoutinePtr<N>> routines;

public:
    void registerRoutine(RoutinePtr<N> r) {
        routines.push_back(r);
    }
    bool cycle(const Matrix<N>& m) {
        bool ret = true;
        for (auto routine = routines.begin(); routine != routines.end(); routine++) {
            ret &= (*routine)->cycle(m);
        }
        return ret;
    }
    /*void updateHistories() {
        for (auto routine = routines.begin(); routine != routines.end(); routine++) {
            (*routine)->updateHistory();
        }
    }
    bool exitConidtion() {
        for (auto routine = routines.begin(); routine != routines.end(); routine++) {
            if (!(*routine)->testExitCondition()) {
                return false;
            }
        }
        return true;
    }*/
    void printResults(std::ostream* out) {
        for (auto routine = routines.begin(); routine != routines.end(); routine++) {
            (*routine)->printResult((*routine)->getCounter(), out);
        }
    }
};

#endif
