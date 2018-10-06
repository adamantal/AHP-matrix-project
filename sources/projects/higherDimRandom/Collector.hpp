#ifndef COLLECTOR_HPP
#define COLLECTOR_HPP

#include <list>

#include "routines/Routine.hpp"

template<size_t N>
class Collector {
private:
    std::list<RoutinePtr<N>> routines;

public:
    void registerRoutine(RoutinePtr<N> r) {
        routines.push_back(r);
    }
    void collectRoutineOutputs(Matrix<N> m) {
        for (auto routine = routines.begin(); routine != routines.end(); routine++) {
            (*routine)->calculate(m);
        }
    }
    void updateHistories() {
        for (auto routine = routines.begin(); routine != routines.end(); routine++) {
            (*routine)->updateHistory();
        }
    }
    void loadFromFolder(std::string s) {
        for (auto routine = routines.begin(); routine != routines.end(); routine++) {
            (*routine)->loadFromFolder(s);
        }
    }
    bool exitConidtion() {
        for (auto routine = routines.begin(); routine != routines.end(); routine++) {
            if (!(*routine)->testExitCondition()) {
                return false;
            }
        }
        return true;
    }
    void printResults(std::ostream* out) {
        for (auto routine = routines.begin(); routine != routines.end(); routine++) {
            (*routine)->printResult(out);
        }
    }
};

#endif
