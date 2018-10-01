#ifndef EXPONENTIALCOUNTER_HPP
#define EXPONENTIALCOUNTER_HPP

#include "TCounter.hpp"

typedef unsigned long long int Ulli;

class ExponentialCounter : public Counter {
private:
    double lambda;

public:
    ExponentialCounter(Ulli un, double unitIncrease = 2.0):Counter(un),lambda(unitIncrease) {
    }

    virtual void increaseUnit() override {
        std::lock_guard<std::mutex> lock(unitMutex);
        unit *= lambda;
    }
};

#endif // EXPONENTIALCOUNTER_HPP
