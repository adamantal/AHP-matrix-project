#ifndef FILE_HPP
#define FILE_HPP

#include <string>
#include <fstream>
#include <memory>

typedef unsigned long long int Ulli;

class File;
typedef std::shared_ptr<File> FilePtr;

class File {
private:
    std::string fileName;
    std::ifstream file;
    Ulli min, max;

    Ulli truncateIndex(Ulli);

public:
    File(std::string name, Ulli min, Ulli max);

    bool get(Ulli i);
    void setToTrue(Ulli i);
};

#endif // FILE_HPP
