/*
- will be a central class that makes the whole process (Processor)
- generate matrices randomly, and (Generator)
- DO NOT STORE them, just calculate some things (Routine)
- and accumulate (Collector) the results
- and then print out, what was found (Printer)
*/

#include <fstream>

#include "TProcessor.hpp"
#include "TRoutine.hpp"
#include "TTestRoutine.hpp"
#include "THistogramRoutine.hpp"
#include "TMultipleHistogramRoutine.hpp"
#include "TLoggingRoutine.hpp"
#include "TPartialSaveRoutine.hpp"

int main () {
    std::cout << "Case 5\n";
    try {
        const unsigned int DIM = 5;
        RoutinePtr<DIM> r1 = std::make_shared<TestRoutine<DIM>> (10000);
        //RoutinePtr<DIM> r2 = std::make_shared<HistogramRoutine<DIM>> (0.01);
        RoutinePtr<DIM> r2 = std::make_shared<MultipleHistogramRoutine<DIM>> (0.01);
        RoutinePtr<DIM> r3 = std::make_shared<LoggingRoutine<DIM>> (10000);
        RoutinePtr<DIM> r4 = std::make_shared<PartialSaveRoutine<DIM>> (r2, 100000);
        std::ostream* out = new std::ofstream("Histogram_" + std::to_string(DIM) + ".out");
        Processor<DIM> p (out);
        p.addRoutine(r1);
        p.addRoutine(r2);
        p.addRoutine(r3);
        p.addRoutine(r4);
        p.process();
        delete out;
    } catch (const char * c) {
        std::cout << c << std::endl;
        std::terminate();
    }
    /*std::cout << "Case 6\n";
    try {
        const unsigned int DIM = 6;
        RoutinePtr<DIM> r1 = std::make_shared<TestRoutine<DIM>> (10000);
        RoutinePtr<DIM> r2 = std::make_shared<MultipleHistogramRoutine<DIM>> (0.01);
        RoutinePtr<DIM> r3 = std::make_shared<LoggingRoutine<DIM>> (2000);
        RoutinePtr<DIM> r4 = std::make_shared<PartialSaveRoutine<DIM>> (r2, 100000);
        std::ostream* out = new std::ofstream("Histogram_" + std::to_string(DIM) + ".out");
        Processor<DIM> p (out);
        p.addRoutine(r1);
        p.addRoutine(r2);
        p.addRoutine(r3);
        p.addRoutine(r4);
        p.process();
        delete out;
    } catch (const char * c) {
        std::cout << c << std::endl;
        std::terminate();
    }*/
    return 0;
}
