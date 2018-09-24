#ifndef ROUTINE_HPP
#define ROUTINE_HPP

#include <memory>

#include "Matrix.hpp"

typedef unsigned long long int Ulli;

template<size_t N>
class Routine;
template<size_t N>
using RoutinePtr = std::shared_ptr<Routine<N>>;

template<size_t N>
class Routine {
protected:
    Routine() {
    }
public:
    virtual void calculate(const Matrix<N>&) = 0;
    virtual void updateHistory() = 0;
    virtual bool testExitCondition() = 0;
    virtual void printResult(std::ostream*) = 0;
};

#endif
