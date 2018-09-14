#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include "Matrix.hpp"

//use const std::vector<double> Matrix<N>::elem !

template<size_t N>
class Generator {
public:
    Generator ();
    Matrix<N> generate();
};

#endif
