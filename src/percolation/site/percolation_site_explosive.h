//
// Created by shahnoor on 6/28/19.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_SITE_EXPLOSIVE_H
#define SQLATTICEPERCOLATION_PERCOLATION_SITE_EXPLOSIVE_H

#include "../percolation.h"
#include "percolation_site_v10.h"


/******************************************************************************
 * Explosive site percolation in square lattice with sum rule and product rule
 */
class SitePercolationExplosive: public SitePercolation_ps_v9{

public:
    ~SitePercolationExplosive() = default;
    SitePercolationExplosive(value_type length);
    std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_explosive_";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        s += std::to_string(length());
        return s;
    }

};




#endif //SQLATTICEPERCOLATION_PERCOLATION_SITE_EXPLOSIVE_H
