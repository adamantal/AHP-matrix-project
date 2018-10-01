#ifndef TESTROUTINE_HPP
#define TESTROUTINE_HPP

#include "TRoutine.hpp"
#include "TCounter.hpp"

typedef unsigned long long int Ulli;

template<size_t N>
class TestRoutine : public Routine<N> {
private:
    Ulli maxCounter;

public:
    TestRoutine(Ulli max):maxCounter(max) {
        CounterPtr c = std::make_shared<Counter> (max);
        Routine<N>::setCounter(c);
    }
    virtual void calculate(Ulli /*count*/, const Matrix<N>& /*m*/) override {
        //do nothing
    }
    virtual void updateHistory(Ulli /*c*/) override {
        //do nothing
    }
    virtual bool testExitCondition(Ulli c) override {
        if (c >= maxCounter)
            return true;
        else
            return false;
    }
    virtual void printResult(Ulli c, std::ostream* out) override {
        *out << c << " matrix generated!\n";
    }
};

#endif // TESTROUTINE_HPP
