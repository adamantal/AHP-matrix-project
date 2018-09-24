#ifndef TESTROUTINE_HPP
#define TESTROUTINE_HPP

#include "Routine.hpp"
#include "counters/Counter.hpp"

typedef unsigned long long int Ulli;

template<size_t N>
class TestRoutine : public Routine<N>, public Counter {
private:
    Ulli maxCounter;

public:
    TestRoutine(Ulli max):Counter(max),maxCounter(max) {
    }

    void calculate(const Matrix<N>& /*m*/) override {
        //do nothing
    }
    virtual bool testExitCondition() override {
        if (getCounter() >= maxCounter)
            return true;
        else
            return false;
    }
    void updateHistory() override {
        incrementCounter();
    }
    virtual void printResult(std::ostream* out) override {
        *out << getCounter() << " matrix generated!\n";
    }
};

#endif // TESTROUTINE_HPP
