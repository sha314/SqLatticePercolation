//
// Created by shahnoor on 2/15/21.
//

#include "test_v13.h"
#include "../v13/lattice_v13.h"
#include "../v13/percolation_v13.h"
#include <iostream>

using namespace std;


void test_lattice() {
    int length = 6;
    auto lattice = Lattice_v13(length);
    lattice.view(0);
    lattice.view(1);
    lattice.view(2);
    lattice.view_relative_index();
    lattice.view_site_gids();
}

void test_detect_wrapping() {
//# take arguments from commandline
    auto sq_lattice_p = SitePercolationL0_v13(5, 18);

//# sq_lattice_p.viewLattice(3)
//# sq_lattice_p.viewCluster()
    int i = 0;
    while (sq_lattice_p.place_one_site()) {
        cout << "p= " << sq_lattice_p.occupation_prob() <<
             " entropy_v1 " << sq_lattice_p.entropy_v1() <<
             " entropy_v2 " << sq_lattice_p.entropy_v2() <<
             " order " << sq_lattice_p.order_param_wrapping() << endl;
        sq_lattice_p.viewLattice(3);
//        sq_lattice_p.viewLattice(4);
//        sq_lattice_p.lattice_ref.print_bonds();
//        sq_lattice_p.viewCluster(1);
        i += 1;
        if (sq_lattice_p.detect_wrapping()) {
            cout << "p= " << sq_lattice_p.occupation_prob() <<
                 " entropy_v1 " << sq_lattice_p.entropy_v1() <<
                 " entropy_v2 " << sq_lattice_p.entropy_v2() <<
                 " order " << sq_lattice_p.order_param_wrapping() << endl;
            cout << "Wrapping detected ***************** <<<" << endl;

        }

//        if (i > 2) break;

    }
//    sq_lattice_p.viewLattice(3);
//    sq_lattice_p.viewLattice(4);
//    sq_lattice_p.viewLattice(1);
//    sq_lattice_p.viewCluster();
}

void test_v13(int argc, char **argv) {
//    test_lattice();
    test_detect_wrapping();
}
