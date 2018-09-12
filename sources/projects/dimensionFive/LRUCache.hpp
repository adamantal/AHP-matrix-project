#ifndef LRUCACHE_HPP
#define LRUCACHE_HPP

#include <list>
#include <map>

template<class T>
class LRUCache {
private:
    unsigned int capacity;
    std::list<T> list;
    std::map<T, typename std::list<T>::iterator> links;

public:
    LRUCache(unsigned int cap);

    void insert(T t);
    bool contains(T t);
};

#include "LRUCache.tpp"

#endif //LRUCACHE_HPP
