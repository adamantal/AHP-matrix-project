#include "File.hpp"

#include <iostream>

File::File(std::string name, Ulli min, Ulli max):fileName(name),min(min),max(max) {
    if (max <= min)
        throw "Error in File creation: min is greater than max!";
    Ulli length = max - min;
    if (length % 8 != 0)
        throw "Error in File creation: length is not divisble by 8!";
    std::ofstream f(fileName);
    char zero;
    zero = zero & 0;
    for (Ulli i = 0; i < length / 8; i++) {
        f << zero;
    }
}

Ulli File::truncateIndex(Ulli i) {
    if (i < this->min || this->max < i) {
        throw "Index is out of bounds!";
    }
    return i - min;
}

bool File::get(Ulli i) {
    i = truncateIndex(i);
    Ulli section = i / 8;
    std::ifstream ifs(fileName);
    char c;
    for (Ulli i = 0; i <= section; i++) {
        ifs >> c;
    }
    return (c >> (i % 8)) & 1;
}

void File::setToTrue(Ulli i) {
    i = truncateIndex(i);
    Ulli section = i / 8;
    std::fstream ifs(fileName, std::ios::in | std::ios::out);
    char c;
    for (Ulli i = 0; i < section; i++) {
        ifs >> c;
    }
    std::streampos oldpos = ifs.tellg();
    ifs >> c;
    ifs.seekg(oldpos);
    ifs << (c |= 1 << (i % 8));
}
