//
// Created by shahnoor on 2/15/21.
//

#include "test_v14.h"
#include "../v13/lattice_v13.h"
#include "../v13/percolation_v13.h"
#include "../v13/run_v13.h"
#include "../v13/percolation_rsbd_v13.h"
#include "../v13/lattice_v14.h"
#include <iostream>

void test_reset();

using namespace std;

void percolation_seed_length_pair(int length, value_type seed);

void test_lattice_v14() {
    int length = 4;
    auto lattice = Lattice_v14(length);
    lattice.view(0);
    lattice.view(1);
    lattice.view(2);
    lattice.view_relative_index();
    lattice.view_site_gids();
    lattice.test_id_index_conversion();
//    lattice.print_site_bond_list();
}


void test_v14(int argc, char **argv) {
    test_lattice_v14();
//    test_detect_wrapping();
//    test_percolation_L0();
//    test_percolation_L1();
//    percolation_seed_length_pairL1(6, 455251785);
//    test_reset();
//    run_ensemble_v13(argc, argv);
}

