//
// Created by shahnoor on 11/9/2017.
//

#ifndef PERCOLATION_CONNECTOR_H
#define PERCOLATION_CONNECTOR_H


#include <vector>

#include "../index/index.h"

// todo delete comment after February 2018

//
///**
// * Takes index and gives connecting bonds and nearest neighbors with and without periodicity
// */
//class Connector{
//    Index _site;
//    value_type _length;
//    value_type min_index {0}, max_index{_length - 1};
//    bool _periodicity;
//    std::vector<BondIndex> _bond_index;
//    std::vector<Index> _neighbor_index;
//
//public:
//    ~Connector() = default;
//    Connector(Index site, value_type length, bool periodicity=true);
//
//    void calculate();
//    std::vector<BondIndex>  connectingBonds(bool peridicity);
//    std::vector<Index>      nearestNeighbors(bool peridicity);
//
//    std::vector<BondIndex>  connectingBonds() const { return _bond_index;};
//    std::vector<Index>      nearestNeighbors() const { return _neighbor_index;};
//};



#endif //PERCOLATION_CONNECTOR_H
