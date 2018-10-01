#ifndef LOGGINGROUTINE_HPP
#define LOGGINGROUTINE_HPP

typedef unsigned long long int Ulli;

#include "TRoutine.hpp"

template<size_t N>
class LoggingRoutine : public Routine<N> {
public:
    LoggingRoutine(Ulli num) {
        Routine<N>::setCounter(std::make_shared<Counter> (num));
    }

    void calculate(const Ulli count, const Matrix<N>& /*m*/) override {
        if (Routine<N>::isActive(count)) {
            std::cout << count << " of matrix has been processed so far.\n";
        }
    }
    void updateHistory(Ulli /*count*/) override {
        //do nothing
    }
    virtual bool testExitCondition(Ulli /*count*/) override {
        return true;
    }
    virtual void printResult(Ulli /*count*/, std::ostream*) override {
        // do nothing
    }
};

#endif // LOGGINGROUTINE_HPP
