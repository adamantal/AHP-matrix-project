#ifndef ROUTINE_HPP
#define ROUTINE_HPP

#include <memory>

#include "Matrix.hpp"
#include "TCounter.hpp"

typedef unsigned long long int Ulli;

template<size_t N>
class Routine;
template<size_t N>
using RoutinePtr = std::shared_ptr<Routine<N>>;

template<size_t N>
class Routine {
private:
    CounterPtr counter = nullptr;

protected:
    Routine() {
    }
    virtual ~Routine() {
    }
    void setCounter(CounterPtr c) {
        counter = c;
    }
    bool isCounterSet() {
        return counter != nullptr;
    }

    // counter accessors
    virtual bool isActive(Ulli count) final {
        return counter->isActive(count);
    }
    virtual Ulli getUnit() final {
        return counter->getUnit();
    }
    virtual void increaseUnit() final {
        counter->increaseUnit();
    }
public:
    virtual Ulli getCounter() final {
        return counter->get();
    }

    // virtual functions:
public:
    virtual void calculate(const Ulli count, const Matrix<N>&) = 0;
    virtual void updateHistory(const Ulli count) = 0;
    virtual bool testExitCondition(const Ulli count) = 0;

    virtual bool cycle(const Matrix<N> m) final {
        Ulli c = counter->getAndIncrement();
        calculate(c, m);
        updateHistory(c);
        return testExitCondition(c);
    }
    virtual void printResult(const Ulli count, std::ostream*) = 0;
};

#endif
