#ifndef LOGGINGROUTINE_HPP
#define LOGGINGROUTINE_HPP

typedef unsigned long long int Ulli;

#include "Routine.hpp"

template<size_t N>
class LoggingRoutine : public Routine<N>, private Counter {
public:
    LoggingRoutine(Ulli num):Counter(num) {}

    void calculate(const Matrix<N>& /*m*/) override {
        if (isActive()) {
            std::cout << getCounter() << " of matrix has been processed so far.\n";
        }
    }
    void updateHistory() override {
        incrementCounter();
    }
    virtual bool testExitCondition() override {
        return true;
    }
    virtual void printResult(std::ostream*) override {
    }
    virtual void loadFromFolder(std::string folder) override {
        History history;
        std::set<std::string> files = scanFolder(folder);
        for (const auto & fileName : files) {
            std::string fullFileName = folder + "/" + fileName;
            std::ifstream file = std::ifstream(fullFileName);
            if (!file.is_open()) {
                throw std::runtime_error("Could not open " + fullFileName + " file!");
            }
            unsigned short i = 0;
            std::string s;
            std::vector<Ulli> seglist;
            while (std::getline(file, s, '\n')) {
                if (i == 1) {
                    std::stringstream test(s);
                    std::string segment;
                    bool first = true;
                    while(std::getline(test, segment, '\t')) {
                        if (!first && !segment.empty()) {
                            seglist.push_back(std::stoi(segment));
                        }
                        first = false;
                    }
                    break;
                }
                i++;
            }
            Ulli sum = std::accumulate(seglist.begin(), seglist.end(), 0);
            history.insert(std::pair<Ulli, std::vector<Ulli>>(sum, seglist));
        }
        Ulli maxElem = history.rbegin()->first;
        setCounter(maxElem);
    }
};

#endif // LOGGINGROUTINE_HPP
