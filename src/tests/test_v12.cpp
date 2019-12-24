//
// Created by shahnoor on 12/23/19.
//

#include <iostream>
#include "test_v12.h"
#include "../lattice/lattice.h"

using namespace std;

void test_v12(int argc, char **argv) {
    cout << "length        = argv[1]" << endl;
    cout << "ensemble size = argv[2]" << endl;
//    cout << "rsbd l        = argv[3]" << endl;

    int length = stoi(argv[1]);
    int ensemble_size = stoi(argv[2]);


    test_v12_lattice(length);
}

void test_v12_lattice(int length) {
    Lattice_v12 lattice(5);

//    lattice.view_as_assigned();
    lattice.view_all();
    lattice.view_sites();
    lattice.view_sites_by_id();

}
