#ifndef LOGGINGROUTINE_HPP
#define LOGGINGROUTINE_HPP

typedef unsigned long long int Ulli;

#include "Routine.hpp"

template<size_t N>
class LoggingRoutine : public Routine<N>, private Counter {
public:
    LoggingRoutine(Ulli num):Counter(num) {}

    void calculate(const Matrix<N>& /*m*/) override {
        if (isActive()) {
            std::cout << getCounter() << " of matrix has been processed so far.\n";
        }
    }
    void updateHistory() override {
        incrementCounter();
    }
    virtual bool testExitCondition() override {
        return true;
    }
    virtual void printResult(std::ostream*) override {
        // do nothing
    }
};

#endif // LOGGINGROUTINE_HPP
