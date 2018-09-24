#ifndef PARTIALSAVEROUTINE_HPP
#define PARTIALSAVEROUTINE_HPP

#include <fstream>

#include "Routine.hpp"

typedef unsigned long long int Ulli;

template<size_t N>
class PartialSaveRoutine : public Routine<N>, private ExponentialCounter {
private:
    RoutinePtr<N> routine;
    Ulli minRound;

public:
    PartialSaveRoutine(RoutinePtr<N> r, Ulli min):ExponentialCounter(1), routine(r), minRound(min) {
        std::system("rm -rf ../results/higherDimRandomRun");
        std::system("mkdir ../results/higherDimRandomRun");
    }
    void calculate(const Matrix<N>& /*m*/) {
        //do nothing
    }
    void updateHistory() {
        incrementCounter();
        if (isActive() && getCounter() > minRound) {
            std::time_t exactTime = std::time(0);
            std::tm* now = std::localtime(&exactTime);
            std::string fileName = "../results/higherDimRandomRun/" + std::to_string(now->tm_mday)
                + "_" + std::to_string(now->tm_hour) + "_" + std::to_string(now->tm_min) + "_"
                + std::to_string(now->tm_sec);
            std::cout << "Writing to " << fileName << std::endl;
            std::ofstream file(fileName);
            routine->printResult(&file);
            file.close();
        }
        if (getCounter() / 20 > getUnit()) {
            increaseUnit();
        }
    }
    bool testExitCondition() {
        return true;
    }
    void printResult(std::ostream*) {
        //do nothing
    }
};

#endif // PARTIALSAVEROUTINE_HPP
