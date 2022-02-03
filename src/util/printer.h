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
#include <map>
#include <initializer_list>
#include <iostream>


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

/**
 * Prints a horizontal barrier in the console.
 * @param n             : number of character on the barrier.
 * @param common        : common character of the barrier.
 * @param special       : special chatater of the barrier.
 * @param special_pos   : where the special character should be placed
 */
void print_h_barrier(size_t n, char common, char special, const std::initializer_list<size_t>& special_pos);

/**
 * Prints a horizontal barrier in the console.
 * @param n             : number of character on the barrier.
 * @param common        : common character of the barrier.
 * @param special       : special chatater of the barrier.
 * @param initial_pos   : Initial posiiton of the special character
 * @param period        : periodic position of the special character.
 */
void print_h_barrier(size_t n, char common, char special, size_t initial_pos, size_t period);

/**
 *
 * Prints a horizontal barrier in the console.
 * @param n             : how many time the middle string is repeated.
 * @param initial       : string that is printed initially.
 * @param middles       : middle string.
 * @param end           : string that is printed at the end.
 */
void print_h_barrier(size_t n, const std::string& initial, const std::string& middles, const std::string& end="\n");

template <class T>
void print_list(std::vector<T> list){
    std::cout << "{";
    for(auto ll: list){
        std::cout << ll << ", ";
    }
    std::cout << "}" << std::endl;
}
#endif //PERCOLATION_PRINTER_H

