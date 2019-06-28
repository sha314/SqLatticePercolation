//
// Created by shahnoor on 6/28/19.
//

#include "test_lattice.h"
#include "../lattice/lattice.h"

void test_sq_lattice_v2(size_t length) {
    SqLattice_v2 lattice(length);
    lattice.init(true, true, true, true);
    lattice.view_sites_by_id();
    lattice.view_sites_by_relative_index();
}
