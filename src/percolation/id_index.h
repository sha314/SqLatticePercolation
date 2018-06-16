//
// Created by shahnoor on 11/7/2017.
//

#ifndef PERCOLATION_ID_INDEX_H
#define PERCOLATION_ID_INDEX_H

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include "../types.h"


/**
 * To track Id and Index of a cluster
 */
class IdIndexTracker{
    std::vector<int> _id;
    std::vector<value_type> _index;
public:
    ~IdIndexTracker() = default;
    IdIndexTracker() = default;
    IdIndexTracker(value_type initial_capacity){
        _id.reserve(initial_capacity);
        _index.reserve(initial_capacity);
    }

    void add(int id, value_type index){
        _id.push_back(id);
        _index.push_back(index);
    }

    void remove(int id){
        if(_id.size() != _index.size()){
            std::cout << "Warning : " << __LINE__ << std::endl;
        }
        for(value_type i{} ; i != _id.size(); ++i){
            if(id == _id[i]){
                _id.erase(_id.begin(), _id.begin()+i);
                _index.erase(_index.begin(), _index.begin()+i);
            }
        }
        resetIndices();
    }

    void remove(value_type index){
        if(_id.size() != _index.size()){
            std::cout << "Warning : " << __LINE__ << std::endl;
        }
        for(value_type i{} ; i != _index.size(); ++i){
            if(index == _index[i]){
                _id.erase(_id.begin(), _id.begin()+i);
                _index.erase(_index.begin(), _index.begin()+i);
            }
        }
        resetIndices();
    }

    void resetIndices(){
        for(value_type i{} ; i != _index.size(); ++i){
            _index[i] = i;
        }
    }



};

#endif //PERCOLATION_ID_INDEX_H
