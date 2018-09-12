#ifndef DISTRIBUTEDFILE_HPP
#define DISTRIBUTEDFILE_HPP

#include <vector>

#include <mutex>
#include "mingw.mutex.h"

#include "LRUCache.hpp"
#include "File.hpp"

typedef unsigned long long int Ulli;

class DistributedFile {
private:
    static const Ulli ALL;
    static const unsigned long CACHE_CAPACITY;

    unsigned int numOfPieces;
    LRUCache<Ulli> cache;

    std::vector<FilePtr> files;
    std::vector<std::unique_ptr<std::mutex>> fileLocks;

    void createFiles();

public:
    DistributedFile(unsigned int p);
    // destructor for file deletion?

    bool get(Ulli i);
    void setToTrue(Ulli i);
};

#endif //DISTRIBUTEDFILE_HPP
