#ifndef PROCESSOR_HPP
#define PROCESSOR_HPP

#include "Generator.hpp"

template<size_t N>
class Processor {
private:
    Generator<N> generator;
    Collector collector;
    Printer printer;

    void printResults();

public:
    Processor ();

    void addRoutine(Routine);
    void process();
};

#endif
