#include "LRUCache.hpp"

template<class T>
LRUCache<T>::LRUCache(unsigned int cap):capacity(cap) {
}

template<class T>
void LRUCache<T>::insert(T t) {
    if (links.find(t) == links.end()) {
        if (list.size() == capacity) {
            T last = list.back();
            list.pop_back();
            links.erase(last);
        }
    } else {
        list.erase(links[t]);
    }
    list.push_front(t);
    links[t] = list.begin();
}

template<class T>
bool LRUCache<T>::contains(T t) {
    return links.find(t) != links.end();
}
