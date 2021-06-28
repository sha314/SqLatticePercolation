//
// Created by shahnoor on 2/15/21.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_RSBD_H
#define SQLATTICEPERCOLATION_PERCOLATION_RSBD_H


#include "lattice_v13.h"
#include "cluster.h"
#include "percolation_v13.h"
#include "../types.h"
#include <set>
#include <algorithm>
#include <random>



class SitePercolationL1_v13: public SitePercolation_v13{
    int x_occupied = 0;
public:

    SitePercolationL1_v13(int length, value_type seed=0, bool generate_seed=true);
    SitePercolationL1_v13(SitePercolationL1_v13&&) = default;
    std::string get_signature() override { return SitePercolation_v13::get_signature() + "L1_";};


    P_STATUS select_site() override;


};



class SitePercolationL2_v13: public SitePercolation_v13{

public:

    SitePercolationL2_v13(int length, value_type seed=0, bool generate_seed=true);
    SitePercolationL2_v13(SitePercolationL2_v13&&) = default;
    std::string get_signature() override { return SitePercolation_v13::get_signature() + "L2_";};


    P_STATUS select_site();

};



#endif //SQLATTICEPERCOLATION_PERCOLATION_H
