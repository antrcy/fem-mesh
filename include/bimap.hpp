#ifndef BIMAP_HPP
#define BIMAP_HPP

#include <unordered_map>
#include <stdexcept>

template<typename key1, typename key2>
class biMapIterator;

/**
 * @brief Provides a double key map similar to what boost::bimap does.
 */
template<typename key1, typename key2>
class biMap{

    using iterator = typename std::unordered_map<key1, key2>::const_iterator;

// Attributes
public:
    int size;
    std::unordered_map<key1, key2> left;
    std::unordered_map<key2, key1> right;

// Methods and constructors
    biMap(int s): size(0) {
        left.reserve(s); right.reserve(s);
    }

    biMap(): size(0) {
        left.reserve(1); right.reserve(1);
    }

    const iterator begin() const {
        return left.begin();
    }

    void insert(const std::pair<key1, key2>& p) {
        left.insert({p.first, p.second});
        right.insert({p.second, p.first});
        size ++;
    }

    int getSize() const {
        return size;
    }

    key1 at_second(key2 key) const {
        return right.at(key);
    }

    key2 at_first(key1 key) const {
        return left.at(key);
    }

    const biMapIterator<key1, key2> get_iterator() const {
        return biMapIterator<key1, key2>(this);
    }
};

/**
 * @brief biMap constant iterator - self explanatory 
 */
template<typename key1, typename key2>
class biMapIterator{

public:
    const biMap<key1, key2>* map;
    typename std::unordered_map<key1, key2>::const_iterator pos;

    biMapIterator(const biMap<key1, key2>* b): map(b) {
        if (map == nullptr){throw std::runtime_error{"biMapIterator was initialized with nullptr"};} 
        pos = map->begin();
    }
    
    key1 get_left() const {
        return pos->first;
    }

    key2 get_right() const {
        return pos->second;
    }
    
    bool end() const {
        return pos == map->left.end();
    }

    void operator++() {
        pos ++;
    }
};

#endif