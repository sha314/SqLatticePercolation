//
// Created by shahnoor on 9/4/18.
//

#include "checking.h"

bool ispresent(const std::vector<value_type> &v, value_type a){
    for(auto x: v){
        if(x == a){
            return true;
        }
    }
    return false;
}