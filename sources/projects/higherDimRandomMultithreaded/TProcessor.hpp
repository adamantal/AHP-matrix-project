#ifndef PROCESSOR_HPP
#define PROCESSOR_HPP

#include <thread>
#include <chrono>

#include "TGenerator.hpp"

template<size_t N>
class Processor;

template<size_t N>
void threadFunction(Collector<N>* c, const Matrix<N>& m, Processor<N>* p) {
    if (c->cycle(m)) {
        p->stopCycle();
    }
    p->threadExited();
}

template<size_t N>
class Processor {
private:
    Generator<N> generator;
    Collector<N> collector;
    std::ostream* output;

    bool stop = false;
    unsigned short activeThreads = 0;
    const unsigned short maxThreads = 4;//2 * std::thread::hardware_concurrency();
    std::mutex threadMutex;

    void printResults() {
        std::cout << "Printing results\n";
        collector.printResults(output);
    }

public:
    Processor (std::ostream* out):output(out) {}

    void stopCycle() {
        stop = true;
    }

    void startNewThread() {
        Matrix<N> m = generator.generate();
        activeThreads++;
        std::thread thread(&threadFunction<N>, &collector, m, this);
        thread.detach();
    }
    void threadExited() {
        std::lock_guard<std::mutex> lock(threadMutex);
        activeThreads--;
        if (!stop) {
            startNewThread();
        }
    }

    void addRoutine(RoutinePtr<N> r) {
        std::cout << "Routine added\n";
        collector.registerRoutine(r);
    }
    void process() {
        std::cout << "Process started\n";
        {
            std::lock_guard<std::mutex> lock(threadMutex);
            for (unsigned short i = 0; i < maxThreads; i++) {
                startNewThread();
            }
        }
        while (true) {
            std::this_thread::sleep_for(std::chrono::seconds(10));
            if (activeThreads == 0) {
                printResults();
                break;
            }
        }
    }
};

#endif
