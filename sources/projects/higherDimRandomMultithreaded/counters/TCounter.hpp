#ifndef COUNTER_HPP
#define COUNTER_HPP

#include <memory>
#include <mutex>

class Counter;
using CounterPtr = std::shared_ptr<Counter>;

typedef unsigned long long int Ulli;

class Counter {
private:
    std::mutex counterMutex;
    Ulli counter = 0;

protected:
    std::mutex unitMutex;
    Ulli unit;

public:
    Counter(Ulli un):unit(un) {
    }

    Ulli get() {
        std::lock_guard<std::mutex> lock(counterMutex);
        return counter;
    }
    /*void increment() {
        std::lock_guard<std::mutex> lock(counterMutex);
        counter++;
    }*/
    Ulli getAndIncrement() {
        std::lock_guard<std::mutex> lock(counterMutex);
        Ulli c = counter;
        counter++;
        return c;
    }

    Ulli getUnit() {
        std::lock_guard<std::mutex> lock(unitMutex);
        return unit;
    }
    virtual void increaseUnit() {
        //do nothing
    }
    bool isActive(Ulli count) {
        std::lock_guard<std::mutex> lockU(unitMutex);
        return (count % unit == 0);
    }
};

#endif // COUNTER_HPP
