#ifndef ROUTINE_HPP
#define ROUTINE_HPP

template<size_t N>
class Routine {
protected:
    Routine();

public:
    virtual void calculate(Matrix<N> m) = 0;
};

#endif
