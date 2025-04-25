#ifndef BIMAP_HPP
#define BIMAP_HPP

#include <unordered_map>
#include <stdexcept>

template<typename rightKey, typename leftKey>
class biMapIterator;

/**
 * @brief Original implementation of a double key map similar to what boost::bimap does.
 */
template<typename rightKey, typename leftKey>
class biMap {

    using iterator = typename std::unordered_map<rightKey, leftKey>::const_iterator;

// Attributes
public:
    int size;
    std::unordered_map<rightKey, leftKey> left;
    std::unordered_map<leftKey, rightKey> right;

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

    void insert(const std::pair<rightKey, leftKey>& p) {
        left.insert({p.first, p.second});
        right.insert({p.second, p.first});
        size ++;
    }

    int getSize() const {
        return size;
    }

    rightKey at_right(leftKey key) const {
        return right.at(key);
    }

    leftKey at_left(rightKey key) const {
        return left.at(key);
    }

    const biMapIterator<rightKey, leftKey> get_iterator() const {
        return biMapIterator<rightKey, leftKey>(this);
    }
};

/**
 * @brief biMap constant iterator - self explanatory 
 */
template<typename rightKey, typename leftKey>
class biMapIterator{

public:
    const biMap<rightKey, leftKey>* map;
    typename std::unordered_map<rightKey, leftKey>::const_iterator pos;

    biMapIterator(const biMap<rightKey, leftKey>* b): map(b) {
        if (map == nullptr){throw std::runtime_error{"biMapIterator was initialized with nullptr"};} 
        pos = map->begin();
    }
    
    rightKey get_left() const {
        return pos->first;
    }

    leftKey get_right() const {
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