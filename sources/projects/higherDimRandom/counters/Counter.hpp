#ifndef COUNTER_HPP
#define COUNTER_HPP

#include <mutex>

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

    Ulli getCounter() {
        std::lock_guard<std::mutex> lock(counterMutex);
        return counter;
    }
    void incrementCounter() {
        std::lock_guard<std::mutex> lock(counterMutex);
        counter++;
    }
    Ulli getUnit() {
        std::lock_guard<std::mutex> lock(unitMutex);
        return unit;
    }
    bool isActive() {
        std::lock_guard<std::mutex> lockC(counterMutex);
        std::lock_guard<std::mutex> lockU(unitMutex);
        return (counter % unit == 0);
    }
    virtual void setCounter(Ulli x) {
        std::lock_guard<std::mutex> lock(counterMutex);
        counter = x;
    }
};

#endif // COUNTER_HPP
