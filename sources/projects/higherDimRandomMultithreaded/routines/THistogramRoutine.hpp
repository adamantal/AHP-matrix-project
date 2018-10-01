#ifndef HISTOGRAMROUTINE_HPP
#define HISTOGRAMROUTINE_HPP

#include <mutex>

#include "TRoutine.hpp"
#include "TExponentialCounter.hpp"

typedef unsigned long long int Ulli;
typedef std::map<Ulli, std::vector<Ulli>> History;

const double EPS = 1e-4;

template<size_t N>
class HistogramRoutine : public Routine<N> {
private:
    const double width;

    std::vector<Ulli> numberOfMatrices;
    History history;

protected: // snake oil
    std::mutex dataMutex;

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
    void saveCurrentStatus(Ulli count) {
        history.insert(std::pair<Ulli, std::vector<Ulli>>(count, numberOfMatrices));
    }

    static bool isClose(double x, double y) {
        //std::cout << "isClose between " << x << " and " << y << " is " << (fabs(x - y) < 1e-6) << std::endl;
        return fabs(x - y) < EPS;
    }
protected:
    void updateIfActive(Ulli count) {
        if (Routine<N>::isActive(count)) {
            std::lock_guard<std::mutex> lock(dataMutex);
            saveCurrentStatus(count);
            for (auto it = history.begin(); it != history.end(); ) {
                if (count / 2 > it->first) {
                    history.erase(it++);
                } else {
                    ++it;
                }
            }
            //std::cout << "The number of history is " << history.size() << "\n"; DEBUG
        }
    }
    void increaseUnitIf(Ulli count) {
        if (count / 20 > Routine<N>::getUnit()) {
            Routine<N>::increaseUnit();
        }
    }

public:
    HistogramRoutine(double w):width(w) {
        Routine<N>::setCounter(std::make_shared<ExponentialCounter> (1));

        unsigned int groups = getMaxNumberOfGroups();
        numberOfMatrices = std::vector<Ulli> (groups);

        std::cout << getMaxInconsistencyEstimate() << " is the inconsistency estimate.\n";
        std::cout << getMaxNumberOfGroups() << " is the max number of groups.\n";
    }

    virtual void updateHistory(Ulli count) override {
        updateIfActive(count);
        increaseUnitIf(count);
    }
protected:
    void calculateCore(double lambda) {
        unsigned int group = getGroupFromRatio(lambda);
        std::lock_guard<std::mutex> lock(dataMutex);
        numberOfMatrices.at(group)++;
    }
public:
    virtual void calculate(Ulli /*count*/, const Matrix<N>& m) override {
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
                        /*std::cout << "Not exit now -- " << float(*resultIt) << " and " << float(*listIt) / n <<
                            " differ between indices " << customHistory.begin()->first << " and " << n << "\n"; DEBUG */
                        std::cout << "The difference compared to the acceptence quota (100%) is " << fabs(float(*resultIt) - float(*listIt) / n) / EPS * 100 << "%\n";
                        return false;
                    }
                }
            }
        }
        std::cout << "checkSingleHistory exit condition TRUE.\n";
        return true;
    }

public:
    virtual bool testExitCondition(Ulli count) override {
        std::lock_guard<std::mutex> lock(dataMutex);
        if (Routine<N>::isActive(count)) {
            return checkSingleHistory(history);
        } else {
            return false;
        }
    }
    virtual void printResult(Ulli /*count*/, std::ostream* out) override {
        *out << "CR" << delimiter;
        double n = getMaxNumberOfGroups();
        for (unsigned int i = 0; i < n; i++) {
            *out << i * width << delimiter;
        }
        *out << "\n";

        std::lock_guard<std::mutex> lock(dataMutex);
        *out << "#ofMatrices" << delimiter;
        for (auto it = numberOfMatrices.begin(); it != numberOfMatrices.end(); it++) {
            *out << *it << delimiter;
        }
        *out << "\n";
    }
};

template<size_t N>
const std::string HistogramRoutine<N>::delimiter = std::string("\t");

#endif // HISTOGRAMROUTINE_HPP
