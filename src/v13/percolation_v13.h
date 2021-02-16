//
// Created by shahnoor on 2/15/21.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_H
#define SQLATTICEPERCOLATION_PERCOLATION_H


#include "lattice_v13.h"

class Percolation_v13{
private:
    int _length = 0;
    unsigned _seed = 0;
    Lattice_v13 _lattice;
    ClusterPool_v13 _cluster_pool;

public:
    Percolation_v13(){

    }
};
#endif //SQLATTICEPERCOLATION_PERCOLATION_H
