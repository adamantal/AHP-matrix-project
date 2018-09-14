#ifndef COLLECTOR_HPP
#define COLLECTOR_HPP

#include "Routine.hpp"

class Collector {
private:
    std::list<RoutinePtr> routines;

public:
    Collector ();

    void registerRoutine(RoutinePtr);
};

#endif
