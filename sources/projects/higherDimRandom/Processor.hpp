#ifndef PROCESSOR_HPP
#define PROCESSOR_HPP

#include "Generator.hpp"

template<size_t N>
class Processor {
private:
    Generator<N> generator;
    Collector<N> collector;
    std::ostream* output;

    void processCore() {
        for (Matrix<N> m = generator.generate(); ; m = generator.generate()) {
            collector.collectRoutineOutputs(m);
            collector.updateHistories();
            if (collector.exitConidtion()) {
                break;
            }
        }
    }
    void printResults() {
        std::cout << "Printing results\n";
        collector.printResults(output);
    }

public:
    Processor (std::ostream* out):output(out) {}

    void addRoutine(RoutinePtr<N> r) {
        std::cout << "Routine added\n";
        collector.registerRoutine(r);
    }
    void process() {
        std::cout << "Process started\n";
        processCore();
        printResults();
    }
};

#endif
