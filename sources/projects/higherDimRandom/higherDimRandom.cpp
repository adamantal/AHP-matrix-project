/*
- will be a central class that makes the whole process (Processor)
- generate matrices randomly, and (Generator)
- DO NOT STORE them, just calculate some things (Routine)
- and accumulate (Collector) the results
- and then print out, what was found (Printer)
*/

#include "Processor.hpp"
#include "routines/Routine.hpp"
#include "routines/TestRoutine.hpp"
#include "routines/HistogramRoutine.hpp"
#include "routines/MultipleHistogramRoutine.hpp"
#include "routines/LoggingRoutine.hpp"
#include "routines/PartialSaveRoutine.hpp"

int main () {
    try {
        RoutinePtr<5> r1 = std::make_shared<TestRoutine<5>> (10000);
        RoutinePtr<5> r2 = std::make_shared<MultipleHistogramRoutine<5>> (0.1);
        RoutinePtr<5> r3 = std::make_shared<LoggingRoutine<5>> (2000);
        RoutinePtr<5> r4 = std::make_shared<PartialSaveRoutine<5>> (r2, 100000);
        Processor<5> p (&std::cout);
        p.addRoutine(r1);
        p.addRoutine(r2);
        p.addRoutine(r3);
        p.addRoutine(r4);
        p.process();
    } catch (const char * c) {
        std::cout << c << std::endl;
    }
    return 0;
}
