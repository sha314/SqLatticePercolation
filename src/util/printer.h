//
// Created by shahnoor on 10/18/2017.
//

#ifndef PERCOLATION_PRINTER_H
#define PERCOLATION_PRINTER_H

#include <ostream>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <vector>


template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> & vec){
    os << '{';
    for(auto a: vec){
        os << a << ',';
    }
    return os << '}';
}


template <typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T> & vec){
    os << '{';
    for(auto a: vec){
        os << a << ',';
    }
    return os << '}';
}


template <typename T>
std::ostream& operator<<(std::ostream& os, const std::unordered_set<T> & vec){
    os << '{';
    for(auto a: vec){
        os << a << ',';
    }
    return os << '}';
}


template <typename K, typename V>
std::ostream& operator<<(std::ostream& os, const std::map<K, V> & m){
    os << '{';
    for(auto a: m){
        os << '(' << a.first << "->" << a.second << "),";
    }
    return os << '}';
};


template <typename K, typename V>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<K, V> & m){
    os << '{';
    for(auto a: m){
        os << '(' << a.first << "->" << a.second << "),";
    }
    return os << '}';
};


#endif //PERCOLATION_PRINTER_H

