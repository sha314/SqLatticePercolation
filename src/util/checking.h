//
// Created by shahnoor on 9/4/18.
//

#ifndef SQLATTICEPERCOLATION_CHECKING_H
#define SQLATTICEPERCOLATION_CHECKING_H


#include "../types.h"
#include <vector>
#include <set>
#include <iostream>


template <class value_type>
bool ispresent(const std::vector<value_type> &v, value_type a){
    for(auto x: v){
        if(x == a){
            return true;
        }
    }
    return false;
}

template <class value_type>
bool ispresent(const std::set<value_type> &vset, value_type a){
    // auto pp = vset.find(a) == vset.end();
    // if(pp){
    //     std::cout << "could not find the element in set" << std::endl;
    // }
    // return pp;
    for(auto x: vset){
        if(x == a){
            return true;
        }
    }
    return false;
}

#endif //SQLATTICEPERCOLATION_CHECKING_H
