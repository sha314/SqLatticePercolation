//
// Created by shahnoor on 10/4/2017.
//

#include "id.h"

template <std::size_t N>
bool operator==(ID<N> id1, ID<N> id2){
    for(std::size_t i {} ; i != N ; ++ i){
        if(id1._id[i] != id2._id[i]){
            return false;
        }
    }
    return true;
}
