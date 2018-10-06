#ifndef MULTIPLEHISTOGRAMROUTINE_HPP
#define MULTIPLEHISTOGRAMROUTINE_HPP

#include "THistogramRoutine.hpp"

typedef unsigned long long int Ulli;

template<size_t N>
class MultipleHistogramRoutine : public HistogramRoutine<N> {
private:
    std::mutex multipleDataMutex;
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
        for (unsigned int i = 0; i < numOfMethods/* + 1*/; i++) {
            std::vector<Ulli> vec = std::vector<Ulli> (groups);
            effs.push_back(vec);
        }
    }

    virtual void updateHistory(const Ulli count) override {
        HistogramRoutine<N>::updateIfActive(count);
        if (Routine<N>::isActive(count)) {
            std::lock_guard<std::mutex> lock(multipleDataMutex);
            for (unsigned int i = 0; i < numOfMethods; i++) {
                histories.at(i).insert(std::pair<Ulli, std::vector<Ulli>>(count, effs.at(i)));
            }
            for (auto& localHistory : histories) {
                for (auto it = localHistory.begin(); it != localHistory.end(); ) {
                    if (count / 2 > it->first) {
                        localHistory.erase(it++);
                    } else {
                        ++it;
                    }
                }
                std::cout << "The number of history in the multiple case is " << localHistory.size() << "\n"; // DEBUG
            }
        }
        HistogramRoutine<N>::increaseUnitIf(count);
    }
    virtual void calculate(const Ulli /*count*/, const Matrix<N>& m) override {
        double x = m.getConsistencyRatio();
        unsigned int group = HistogramRoutine<N>::getGroupFromRatio(x);
        HistogramRoutine<N>::calculateCore(x);

        bool eigen = m.testParetoOptimality(filterType::EigenVectorMethod);
        bool spanTree = m.testParetoOptimality(filterType::AverageSpanTreeMethod);
        bool cosine = m.testParetoOptimality(filterType::CosineMethod);

        {
            std::lock_guard<std::mutex> lock(multipleDataMutex);
            if (eigen)
                effs.at(0).at(group)++;
            if (spanTree)
                effs.at(1).at(group)++;
            if (cosine)
                effs.at(2).at(group)++;
            /*if (eigen && spanTree && cosine)
                effs.at(3).at(group)++;*/
        }
    }
    virtual bool testExitCondition(Ulli count) override {
        if (!HistogramRoutine<N>::testExitCondition(count)) {
            return false;
        } else {
            std::lock_guard<std::mutex> lock(multipleDataMutex);
            return std::all_of(histories.begin(), histories.end(), &HistogramRoutine<N>::checkSingleHistory);
        }
    }
    virtual void printResult(Ulli count, std::ostream* out) override {
        HistogramRoutine<N>::printResult(count, out);
        std::lock_guard<std::mutex> lock(multipleDataMutex);
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

        /* *out << "#allOfThree" << HistogramRoutine<N>::delimiter;
        for (auto it = effs.at(3).begin(); it != effs.at(3).end(); it++) {
            *out << *it << HistogramRoutine<N>::delimiter;
        }
        *out << "\n";*/
    }
};

#endif // MULTIPLEHISTOGRAMROUTINE_HPP
