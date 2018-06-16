//
// Created by shahnoor on 10/4/2017.
//

#ifndef SITEPERCOLATION_ID_H
#define SITEPERCOLATION_ID_H

#include <vector>
#include <array>
#include <iostream>

template <std::size_t N>
struct ID{
    std::array<int, N> _id;

    ~ID() = default;
    ID() = default;
    ID(ID&) = default;
    ID(ID&&)= default;
    ID& operator=(ID&)= default;
    ID& operator=(ID&&)= default;

    ID(std::vector<int> id){
        if(id.size() != N){
            std::cout << "Cannot construct ID : numberOfBonds mismatched : line " << __LINE__ << std::endl;
        }
        for(std::size_t i {} ; i != N ; ++ i){
            _id[i] = id[i];
        }
    }
};

template <std::size_t N>
bool operator==(ID<N> id1, ID<N> id2);


#endif //SITEPERCOLATION_ID_H
