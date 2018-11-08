//
// Created by shahnoor on 8/16/18.
//

#include "printer.h"
#include <iostream>

using namespace std;


void print_h_barrier(size_t n, char common, char special, const initializer_list<size_t>& special_pos){
    size_t k{};
    vector<size_t> v(special_pos);
    for(size_t i{}; i < n ; ++i){
        if(k < special_pos.size() && i == v[k]){
            ++k;
            cout << special;
            continue;
        }
        cout << common;
    }
    cout << endl; // end of barrier
}

void print_h_barrier(size_t n, char common, char special, size_t initial_pos, size_t period){
    size_t k{};

    for(size_t i{}; i < n ; ++i){
        if(i == (initial_pos + period * k)){
            ++k;
            cout << special;
            continue;
        }
        cout << common;
    }
    cout << endl; // end of barrier
}

void print_h_barrier(size_t n, const string& initial, const string& middles, const string& end){
    cout << initial;
    for(size_t i{}; i < n ; ++i){
        cout << middles;
    }
    cout << end; // end of barrier
}
