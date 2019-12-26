//
// Created by shahnoor on 12/23/19.
//

#include <iostream>
#include "test_v12.h"
#include "../lattice/lattice.h"
#include "../percolation/percolation_v12.h"

using namespace std;

void test_v12(int argc, char **argv) {
    cout << "length        = argv[1]" << endl;
    cout << "ensemble size = argv[2]" << endl;
//    cout << "rsbd l        = argv[3]" << endl;

    int length = stoi(argv[1]);
    int ensemble_size = stoi(argv[2]);


//    test_v12_lattice(length);
    test_v12_percolation(length);
}

void test_v12_lattice(int length) {
    Lattice_v12 lattice(length);

//    lattice.view_as_assigned();
    lattice.view_all();
//    lattice.view_sites();
//    lattice.view_sites_by_id();
//    lattice.view_sites_by_gid();
//    lattice.view_by_id();
//    lattice.view_by_gid();
//    lattice.view_sites_by_relative_index();
//    lattice.view_by_relative_index();

    lattice.view_sites_list();
    lattice.view_bonds_list();

//    cout << lattice.getSiteLocation(0) << endl;
//    cout << lattice.getSiteLocation(5) << endl;
//    cout << lattice.getSiteLocation(12) << endl;
//    cout << lattice.getSiteLocation(22) << endl;


//    cout << lattice.getBondLocation(11) << endl;
//    cout << lattice.getBondLocation(17) << endl;
//    cout << lattice.getBondLocation(22) << endl;
//    cout << lattice.getBondLocation(35) << endl;
//    cout << lattice.getBondLocation(51) << endl;
//
//
//    cout << lattice.getBond(0)<< endl;
//    cout << lattice.getBond(11)<< endl;
//    cout << lattice.getBond(17)<< endl;
//    cout << lattice.getBond(22)<< endl;
//    cout << lattice.getBond(35)<< endl;
}


void test_v12_percolation(int length) {
    SqLatticeRegularSite percolation(length);

    percolation.init();

//    percolation.viewCluster();
//    percolation.viewLattice_by_id();
//    percolation.viewLattice_by_gid();

//    percolation.occupy();
    int i=0;
    while(percolation.occupy()) {
        cout << "********************** last site " << percolation.lastSite() << endl;
    percolation.viewCluster();
//    percolation.viewLattice_by_id();
//        percolation.viewLattice_by_gid();
//    percolation.viewLattice();
        percolation.viewLattice_by_relative_index();

        ++i;
        if(i >= 3) break;
    }


}
