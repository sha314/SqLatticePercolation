//
// Created by shahnoor on 2/15/21.
//

#include "test_v13.h"
#include "../v13/lattice_v13.h"
#include "../v13/percolation_v13.h"
#include "../v13/run_v13.h"
#include <iostream>

void test_reset();

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
    auto sq_lattice_p = SitePercolationL0_v13(6, 18, false);
//    sq_lattice_p.setRandomState(0, true);

//# sq_lattice_p.viewLattice(3)
//# sq_lattice_p.viewCluster()
    int i = 0;
    while (sq_lattice_p.place_one_site()) {
        cout << "p= " << sq_lattice_p.occupation_prob() <<
             " entropy_v1 " << sq_lattice_p.entropy_v1() <<
             " entropy_v2 " << sq_lattice_p.entropy_v2() <<
             " order " << sq_lattice_p.order_param_wrapping() << endl;
//        sq_lattice_p.viewLattice(3);
//        sq_lattice_p.viewLattice(4);
//        sq_lattice_p.lattice_ref.print_bonds();
        sq_lattice_p.viewCluster(1);
        i += 1;
        if (sq_lattice_p.detect_wrapping()) {
            cout << "p= " << sq_lattice_p.occupation_prob() <<
                 " entropy_v1 " << sq_lattice_p.entropy_v1() <<
                 " entropy_v2 " << sq_lattice_p.entropy_v2() <<
                 " order " << sq_lattice_p.order_param_wrapping() << endl;
            cout << "Wrapping detected ***************** <<<" << endl;
            break;
        }

//        if (i > 8) break;

    }
    sq_lattice_p.viewLattice(3);
//    sq_lattice_p.viewLattice(4);
//    sq_lattice_p.viewLattice(1);
//    sq_lattice_p.viewCluster();
}

void test_v13(int argc, char **argv) {
//    test_lattice();
//    test_detect_wrapping();
//    test_reset();
    run_ensemble_v13(argc, argv);
}

void test_reset() {
//# take arguments from commandline
    auto sq_lattice_p = SitePercolationL0_v13(6, 18, false);
//    sq_lattice_p.setRandomState(0, true);

//    sq_lattice_p.viewLattice(3);
//    sq_lattice_p.viewCluster(1);
    sq_lattice_p.reset();
//    sq_lattice_p.viewLattice(3);
//    sq_lattice_p.viewCluster(1);

    int i = 0;
    while (sq_lattice_p.place_one_site()) {
        cout << "p= " << sq_lattice_p.occupation_prob() <<
             " entropy_v1 " << sq_lattice_p.entropy_v1() <<
             " entropy_v2 " << sq_lattice_p.entropy_v2() <<
             " order_largest " << sq_lattice_p.order_param_largest_clstr() <<
             " order_wrapping " << sq_lattice_p.order_param_wrapping() << endl;
//        sq_lattice_p.viewLattice(3);
//        sq_lattice_p.viewLattice(4);
//        sq_lattice_p.lattice_ref.print_bonds();
//        sq_lattice_p.viewCluster(1);
        i += 1;
        if (sq_lattice_p.detect_wrapping()) {
//            cout << "p= " << sq_lattice_p.occupation_prob() <<
//                 " entropy_v1 " << sq_lattice_p.entropy_v1() <<
//                 " entropy_v2 " << sq_lattice_p.entropy_v2() <<
//                 " order_largest " << sq_lattice_p.order_param_largest_clstr() <<
//                 " order_wrapping " << sq_lattice_p.order_param_wrapping() << endl;
//            cout << "Wrapping detected ***************** <<<" << endl;
//            break;
        }

//        if (i > 8) break;

    }
    sq_lattice_p.viewLattice(3);
//    sq_lattice_p.viewLattice(4);
//    sq_lattice_p.viewLattice(1);
//    sq_lattice_p.viewCluster();
}

void run_ensemble_v13(int argc, char **argv){
    cout << "length        = argv[1]" << endl;
    cout << "ensemble size = argv[2]" << endl;
//    cout << "rsbd l        = argv[3]" << endl;

    int length = stoi(argv[1]);
    int ensemble_size = stoi(argv[2]);


//    test_v9_lattice(argc, argv);
//    test_v10_number_of_clusters_per_site(argc, argv);
//    test_v10_lattice(argc, argv);

//    run_v10_regular(length, ensemble_size);

//    int rsbd_l = stoi(argv[3]);
    run_v13_rsbd_L0(length, ensemble_size);
//    run_v13_rsbd<SitePercolationL0_v13>(length, ensemble_size);

}
