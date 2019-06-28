//
// Created by shahnoor on 6/28/19.
//

#include "test_v11.h"
#include "../lattice/lattice.h"
#include "../percolation/site/percolation_site_v11.h"

void test_sq_lattice_v2(size_t length) {
    SqLattice_v2 lattice(length);
    lattice.init(true, true, true, true);
    lattice.view_sites_by_id();
    lattice.view_sites_by_relative_index();
}

void test_percolation(size_t length){
    SitePercolation_ps_v11 sp(length);
    sp.setRandomState(0, false);
    sp.init(false);

    sp.viewLattice();
}
