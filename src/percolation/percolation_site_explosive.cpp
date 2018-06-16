//
// Created by shahnoor on 02-Apr-18.
//

#include "percolation.h"

SitePercolationExplosive::SitePercolationExplosive(value_type length)
        : SitePercolation_ps_v8(length)
{
    std::cout << "Constructing SitePercolationExplosive object : line " << __LINE__ << endl;
}