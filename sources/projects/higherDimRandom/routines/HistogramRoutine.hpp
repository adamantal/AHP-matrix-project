#ifndef HISTOGRAMROUTINE_HPP
#define HISTOGRAMROUTINE_HPP

#include "Routine.hpp"
#include "ExponentialCounter.hpp"
#include "../FileUtil.hpp"

typedef unsigned long long int Ulli;
typedef std::map<Ulli, std::vector<Ulli>> History;

const double EPS = 1e-5;

template<size_t N>
class HistogramRoutine : public Routine<N>, protected ExponentialCounter {
private:
    const double width;

    std::vector<Ulli> numberOfMatrices;
    History history;

protected:
    static const std::string delimiter;

private:
    static double getMaxInconsistencyEstimate() {
        double lambdaMaxEstimate = 1.0f + 0.5f * (N - 1) * (9.0f + 1.0f / 9.0f);
        return (lambdaMaxEstimate - N) / (Matrix<N>::ConsistencyIndex[N] - N);
    }
protected:
    unsigned int getMaxNumberOfGroups() {
        double incons = getMaxInconsistencyEstimate();
        return ceil(incons / width);
    }
    unsigned int getGroupFromRatio(double x) const {
        double r = floor(x / width);
        return isClose(r, 0.0f) ? 0 : r;
    }
private:
    void saveCurrentStatus() {
        history.insert(std::pair<Ulli, std::vector<Ulli>>(getCounter(), numberOfMatrices));
    }

    static bool isClose(double x, double y) {
        //std::cout << "isClose between " << x << " and " << y << " is " << (fabs(x - y) < 1e-6) << std::endl;
        return fabs(x - y) < EPS;
    }
protected:
    void updateIfActive() {
        if (isActive()) {
            saveCurrentStatus();
            for (auto it = history.begin(); it != history.end(); ) {
                if (getCounter() / 2 > it->first) {
                    history.erase(it++);
                } else {
                    ++it;
                }
            }
            //std::cout << "The number of history is " << history.size() << "\n"; DEBUG
        }
    }
    void increaseUnitIf() {
        if (getCounter() / 20 > getUnit()) {
            increaseUnit();
        }
    }

public:
    HistogramRoutine(double w):ExponentialCounter(1),width(w) {
        unsigned int groups = getMaxNumberOfGroups();
        numberOfMatrices = std::vector<Ulli> (groups);

        std::cout << getMaxInconsistencyEstimate() << " is the inconsistency estimate.\n";
        std::cout << getMaxNumberOfGroups() << " is the max number of groups.\n";
    }

    virtual void updateHistory() override {
        incrementCounter();
        updateIfActive();
        increaseUnitIf();
    }
protected:
    void calculateCore(double lambda) {
        unsigned int group = getGroupFromRatio(lambda);
        numberOfMatrices.at(group)++;
    }
public:
    virtual void calculate(const Matrix<N>& m) override { // thread-safe
        double x = m.getConsistencyRatio();
        calculateCore(x);
    }
protected:
    static bool checkSingleHistory(History& customHistory) {
        std::cout << "checkSingleHistory exit condition testing starts.\n";
        auto it = customHistory.begin();
        std::list<double> results;
        {
            Ulli n = it->first;
            for (auto listIt = it->second.begin(); listIt != it->second.end(); listIt++) {
                results.push_back(float(*listIt) / n);
            }
        }
        for (it++; it != customHistory.end() ; it++) {
            Ulli n = it->first;
            {
                auto listIt = it->second.begin();
                auto resultIt = results.begin();
                for (; listIt != it->second.end(); listIt++, resultIt++) {
                    if (!isClose(float(*resultIt), float(*listIt) / n)) {
                        std::cout << "Not exit now -- " << float(*resultIt) << " and " << float(*listIt) / n <<
                            " differ between indices " << customHistory.begin()->first << " and " << n << "\n";
                        std::cout << "The difference compared to the acceptence quota (100%) is " << fabs(float(*resultIt) - float(*listIt) / n) / EPS * 100 << "%\n";
                        return false;
                    }
                }
            }
        }
        std::cout << "checkSingleHistory exit condition TRUE.\n";
        return true;
    }
protected:
    bool testExitConditionCore() {
        if (isActive()) {
            return checkSingleHistory(history);
        } else {
            return false;
        }
    }

public:
    virtual bool testExitCondition() override {
        return testExitConditionCore();
    }
    virtual void printResult(std::ostream* out) override {
        *out << "CR" << delimiter;
        double n = getMaxNumberOfGroups();
        for (unsigned int i = 0; i < n; i++) {
            *out << i * width << delimiter;
        }
        *out << "\n";

        *out << "#ofMatrices" << delimiter;
        for (auto it = numberOfMatrices.begin(); it != numberOfMatrices.end(); it++) {
            *out << *it << delimiter;
        }
        *out << "\n";
    }
public:
    virtual void loadFromFolder(std::string folder) override {
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
            if (seglist.size() != getMaxNumberOfGroups()) {
                std::cout << seglist.size() << std::endl;
                throw std::runtime_error("The read numers don't have the expected size!");
            }
            Ulli sum = std::accumulate(seglist.begin(), seglist.end(), 0);
            history.insert(std::pair<Ulli, std::vector<Ulli>>(sum, seglist));
        }
        Ulli maxElem = history.rbegin()->first;
        numberOfMatrices = history.rbegin()->second;
        setCounter(maxElem);
    }
};

template<size_t N>
const std::string HistogramRoutine<N>::delimiter = std::string("\t");

#endif // HISTOGRAMROUTINE_HPP
