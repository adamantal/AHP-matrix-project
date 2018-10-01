#ifndef PARTIALSAVEROUTINE_HPP
#define PARTIALSAVEROUTINE_HPP

#include <fstream>

#include "TRoutine.hpp"

typedef unsigned long long int Ulli;

template<size_t N>
class PartialSaveRoutine : public Routine<N> {
private:
    RoutinePtr<N> routine;
    Ulli minRound;

public:
    PartialSaveRoutine(RoutinePtr<N> r, Ulli min):routine(r), minRound(min) {
        std::system("rm -rf ../results/higherDimRandomRun");
        std::system("mkdir ../results/higherDimRandomRun");

        Routine<N>::setCounter(std::make_shared<ExponentialCounter>(1));
    }
    void calculate(Ulli /*count*/, const Matrix<N>& /*m*/) {
        //do nothing
    }
    void updateHistory(Ulli count) {
        if (Routine<N>::isActive(count) && count > minRound) {
            std::time_t exactTime = std::time(0);
            std::tm* now = std::localtime(&exactTime);
            std::string fileName = "../results/higherDimRandomRun/" + std::to_string(now->tm_mday)
                + "_" + std::to_string(now->tm_hour) + "_" + std::to_string(now->tm_min) + "_"
                + std::to_string(now->tm_sec);
            std::cout << "Writing to " << fileName << std::endl;
            std::ofstream file(fileName);
            routine->printResult(count, &file);
            file.close();
        }
        if (count / 20 > Routine<N>::getUnit()) {
            Routine<N>::increaseUnit();
        }
    }
    bool testExitCondition(const Ulli /*count*/) {
        return true;
    }
    void printResult(const Ulli /*count*/, std::ostream*) {
        //do nothing
    }
};

#endif // PARTIALSAVEROUTINE_HPP
