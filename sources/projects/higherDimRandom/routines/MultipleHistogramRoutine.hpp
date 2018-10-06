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
                //std::cout << "The number of history in the multiple case is " << localHistory.size() << "\n"; DEBUG
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
            std::vector<Ulli> simpleSegList;
            std::vector<std::vector<Ulli>> seglists;
            for (unsigned short i = 0; i < numOfMethods; i++) {
                std::vector<Ulli> vec;
                seglists.push_back(vec);
            }
            while (std::getline(file, s, '\n')) {
                if (i == 1) {
                    std::stringstream test(s);
                    std::string segment;
                    bool first = true;
                    while(std::getline(test, segment, '\t')) {
                        if (!first && !segment.empty()) {
                            simpleSegList.push_back(std::stoi(segment));
                        }
                        first = false;
                    }
                }
                if (i == 2) {
                    std::stringstream test(s);
                    std::string segment;
                    bool first = true;
                    while(std::getline(test, segment, '\t')) {
                        if (!first && !segment.empty()) {
                            seglists.at(0).push_back(std::stoi(segment));
                        }
                        first = false;
                    }
                }
                if (i == 3) {
                    std::stringstream test(s);
                    std::string segment;
                    bool first = true;
                    while(std::getline(test, segment, '\t')) {
                        if (!first && !segment.empty()) {
                            seglists.at(1).push_back(std::stoi(segment));
                        }
                        first = false;
                    }
                }
                if (i == 4) {
                    std::stringstream test(s);
                    std::string segment;
                    bool first = true;
                    while(std::getline(test, segment, '\t')) {
                        if (!first && !segment.empty()) {
                            seglists.at(2).push_back(std::stoi(segment));
                        }
                        first = false;
                    }
                    break;
                }
                i++;
            }
            if (N == 5 && seglists.at(0).size() != HistogramRoutine<N>::getMaxNumberOfGroups() &&
                          seglists.at(1).size() != HistogramRoutine<N>::getMaxNumberOfGroups() &&
                          seglists.at(2).size() != HistogramRoutine<N>::getMaxNumberOfGroups()) {
                std::cout << seglists.at(0).size() << " " << seglists.at(1).size() << " " << seglists.at(2).size() << std::endl;
                throw std::runtime_error("The read numers don't have the expected size!");
            }
            Ulli sum = std::accumulate(simpleSegList.begin(), simpleSegList.end(), 0);
            histories.at(0).insert(std::pair<Ulli, std::vector<Ulli>>(sum, seglists.at(0)));
            histories.at(1).insert(std::pair<Ulli, std::vector<Ulli>>(sum, seglists.at(1)));
            histories.at(2).insert(std::pair<Ulli, std::vector<Ulli>>(sum, seglists.at(2)));
        }
        for (unsigned int i = 0; i < numOfMethods; i++) {
            std::vector<Ulli> eff = histories.at(i).rbegin()->second;
            effs.at(i) = eff;
        }
        HistogramRoutine<N>::loadFromFolder(folder);
    }
};

#endif // MULTIPLEHISTOGRAMROUTINE_HPP
