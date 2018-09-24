#ifndef COUNTER_HPP
#define COUNTER_HPP

typedef unsigned long long int Ulli;

class Counter {
private:
    Ulli counter = 0;

protected:
    Ulli unit;

public:
    Counter(Ulli un):unit(un) {
    }

    Ulli getCounter() {
        return counter;
    }
    void incrementCounter() {
        counter++;
    }
    Ulli getUnit() {
        return unit;
    }
    bool isActive() {
        return (counter % unit == 0);
    }
};

#endif // COUNTER_HPP
