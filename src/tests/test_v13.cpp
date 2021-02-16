//
// Created by shahnoor on 2/15/21.
//

#include "test_v13.h"
#include "../v13/lattice_v13.h"

void test_lattice() {
    int length = 6;
    auto lattice = Lattice_v13(length);
    lattice.view(0);
    lattice.view(1);
    lattice.view(2);
    lattice.view_relative_index();
    lattice.view_site_gids();
}

void test_v13(int argc, char **argv) {
    test_lattice();
}
