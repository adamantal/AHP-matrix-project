#ifndef EXPONENTIALCOUNTER_HPP
#define EXPONENTIALCOUNTER_HPP

#include "Counter.hpp"

typedef unsigned long long int Ulli;

class ExponentialCounter : public Counter {
private:
    double lambda;

public:
    ExponentialCounter(Ulli un, double unitIncrease = 2.0):Counter(un),lambda(unitIncrease) {
    }

    void increaseUnit() {
        unit *= lambda;
    }
};

#endif // EXPONENTIALCOUNTER_HPP
