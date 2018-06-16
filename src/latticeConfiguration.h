//
// Created by shahnoor on 10/25/17.
//

#ifndef PERCOLATION_LATTICECONFIGURATION_H
#define PERCOLATION_LATTICECONFIGURATION_H


#include <vector>
#include "index/index.h"


/****
 * functions return the active site indices
 */
namespace LatticeSiteConfiguration{
    // for 6x6 lattice
    std::vector<Index> configuration_6x6(int x);

    std::vector<Index> configuration_6x6_0();
    std::vector<Index> configuration_6x6_1();
    std::vector<Index> configuration_6x6_2();
    std::vector<Index> configuration_6x6_3();
    std::vector<Index> configuration_6x6_4();
    std::vector<Index> configuration_6x6_5();


    // for 12x12 lattice
    std::vector<Index> configuration_12x12(int x);

    std::vector<Index> configuration_12x12_1();
    std::vector<Index> configuration_12x12_2();
    std::vector<Index> configuration_12x12_3();

}

namespace LatticeBondConfiguration{
    // for 6x6 lattice
    std::vector<BondIndex> configuration_6x6(int x);

    std::vector<BondIndex> configuration_6x6_0();
    std::vector<BondIndex> configuration_6x6_1();
//    std::vector<BondIndex> configuration_6x6_2();
//    std::vector<BondIndex> configuration_6x6_3();
//    std::vector<BondIndex> configuration_6x6_4();
//    std::vector<BondIndex> configuration_6x6_5();


    // for 12x12 lattice
//    std::vector<Index> configuration_12x12(int x);
//
//    std::vector<Index> configuration_12x12_1();
//    std::vector<Index> configuration_12x12_2();
//    std::vector<Index> configuration_12x12_3();

}

#endif //PERCOLATION_LATTICECONFIGURATION_H
