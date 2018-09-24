#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include <random>
#include <chrono>

#include "Matrix.hpp"

#include "Collector.hpp"

//use const std::vector<double> Matrix<N>::elem !
typedef unsigned short Ush;

enum GeneratorOption {
    DEFAULT,
    MINIMAL_PERMUATED,
    SAATY_CONSISTENT
};

template<size_t N>
class Generator {
private:
    std::default_random_engine eng;
    GeneratorOption opt;

    std::vector<Ush> generateRandomVector() {
        std::vector<Ush> ret;
        std::uniform_int_distribution<Ush> distr (0, Matrix<N>::elem.size() - 1);
        for (Ush i = 0; i < N * (N - 1) / 2; i++) {
            ret.push_back(distr(eng));
        }
        return ret;
    }

public:
    Generator (GeneratorOption option = DEFAULT):opt(option)  {
        eng.seed(std::chrono::system_clock::now().time_since_epoch().count());
    }

    Matrix<N> generate() {
        std::vector<Ush> ret = generateRandomVector();

        if (opt == MINIMAL_PERMUATED) {
            return Matrix<N>(ret).getItsMinimalPermutate();
        } else if (opt == SAATY_CONSISTENT) {
            if (Matrix<N>(ret).testParetoOptimality(filterType::Consistency)) {
                return ret;
            } else {
                return generate();
            }
        } else {
            return ret;
        }
    }
};

#endif
