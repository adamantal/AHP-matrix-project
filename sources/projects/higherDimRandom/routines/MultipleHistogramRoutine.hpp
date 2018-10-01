#ifndef MULTIPLEHISTOGRAMROUTINE_HPP
#define MULTIPLEHISTOGRAMROUTINE_HPP

#include "HistogramRoutine.hpp"

typedef unsigned long long int Ulli;

template<size_t N>
class MultipleHistogramRoutine : public HistogramRoutine<N> {
private:
    std::vector<std::vector<Ulli>> effs;
    std::vector<History> histories;

    static const unsigned int numOfMethods = 3;

public:
    MultipleHistogramRoutine(double w):HistogramRoutine<N>(w) {
        //init histories
        for (unsigned short i = 0; i < numOfMethods; i++) {
            std::map<Ulli, std::vector<Ulli>> m;
            histories.push_back(m);
        }

        //init vectors
        unsigned int groups = HistogramRoutine<N>::getMaxNumberOfGroups();
        for (unsigned int i = 0; i < numOfMethods; i++) {
            std::vector<Ulli> vec = std::vector<Ulli> (groups);
            effs.push_back(vec);
        }
    }

    virtual void updateHistory() override {
        HistogramRoutine<N>::incrementCounter();
        HistogramRoutine<N>::updateIfActive();
        if (Counter::isActive()) {
            for (unsigned int i = 0; i < numOfMethods; i++) {
                histories.at(i).insert(std::pair<Ulli, std::vector<Ulli>>(Counter::getCounter(), effs.at(i)));
            }
            for (auto& localHistory : histories) {
                for (auto it = localHistory.begin(); it != localHistory.end(); ) {
                    if (Counter::getCounter() / 2 > it->first) {
                        localHistory.erase(it++);
                    } else {
                        ++it;
                    }
                }
                std::cout << "The number of history in the multiple case is " << localHistory.size() << "\n";
            }
        }
        HistogramRoutine<N>::increaseUnitIf();
    }
    virtual void calculate(const Matrix<N>& m) override {
        double x = m.getConsistencyRatio();
        unsigned int group = HistogramRoutine<N>::getGroupFromRatio(x);
        HistogramRoutine<N>::calculateCore(x);

        // EigenVectorMethod
        if (m.testParetoOptimality(filterType::EigenVectorMethod)) {
            effs.at(0).at(group)++;
        }
        // AverageSpanTreeMethod
        if (m.testParetoOptimality(filterType::AverageSpanTreeMethod)) {
            effs.at(1).at(group)++;
        }
        // CosineMethod
        if (m.testParetoOptimality(filterType::CosineMethod)) {
            effs.at(2).at(group)++;
        }
    }
    virtual bool testExitCondition() override {
        if (!HistogramRoutine<N>::testExitConditionCore()) {
            return false;
        } else {
            return std::all_of(histories.begin(), histories.end(), &HistogramRoutine<N>::checkSingleHistory);
        }
    }
    virtual void printResult(std::ostream* out) override {
        HistogramRoutine<N>::printResult(out);
        *out << "#ofEigenEffs" << HistogramRoutine<N>::delimiter;
        for (auto it = effs.at(0).begin(); it != effs.at(0).end(); it++) {
            *out << *it << HistogramRoutine<N>::delimiter;
        }
        *out << "\n";

        *out << "#ofAvgSpanTreeEffs" << HistogramRoutine<N>::delimiter;
        for (auto it = effs.at(1).begin(); it != effs.at(1).end(); it++) {
            *out << *it << HistogramRoutine<N>::delimiter;
        }
        *out << "\n";

        *out << "#ofCosineEffs" << HistogramRoutine<N>::delimiter;
        for (auto it = effs.at(2).begin(); it != effs.at(2).end(); it++) {
            *out << *it << HistogramRoutine<N>::delimiter;
        }
        *out << "\n";
    }
};

#endif // MULTIPLEHISTOGRAMROUTINE_HPP
