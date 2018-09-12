#include "DistributedFile.hpp"

#include <iostream>

const Ulli DistributedFile::ALL = 2015993900449 + 7;
const unsigned long DistributedFile::CACHE_CAPACITY = ALL / 1000;

DistributedFile::DistributedFile(unsigned int p):numOfPieces(p), cache(CACHE_CAPACITY) {
    createFiles();
}

void DistributedFile::createFiles() {
    std::cout << "Creating files...\n";
    int res = std::system("mkdir tmp/");
    (void)res;
    for (unsigned int i = 0; i < numOfPieces; i++) {
        if (i % (numOfPieces / 100) == 0) {
            std::cout << "Percentage is " << i / (numOfPieces / 100) << std::endl;
        }
        std::string s = std::string("tmp/file_") + std::to_string(i);
        FilePtr f = std::make_shared<File> (s, ALL / numOfPieces * i, ALL / numOfPieces * (i + 1));
        files.push_back(f);

        fileLocks.push_back(std::make_unique<std::mutex> ());
    }
}

bool DistributedFile::get(Ulli i) {
    if (cache.contains(i)) {
        return true;
    }
    std::lock_guard<std::mutex> lock_guard(*fileLocks[i / numOfPieces]);
    return files[i / numOfPieces]->get(i);
}

void DistributedFile::setToTrue(Ulli i) {
    if (!cache.contains(i)) {
        std::lock_guard<std::mutex> lock_guard(*fileLocks[i / numOfPieces]);
        files[i / numOfPieces]->setToTrue(i);
        cache.insert(i);
    }
}
