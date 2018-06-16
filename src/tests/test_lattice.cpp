//
// Created by shahnoor on 10/27/17.
//

#include "test_lattice.h"
#include "../lattice/lattice.h"


void TEST_Lattice::test_lattice(){
    SqLattice lattice(5);

//    _lattice.view_sites();
    lattice.activateAllSite();
//    _lattice.activate_site({0, 1});
//    _lattice.activate_site({2, 3});
    lattice.view_sites();

//    _lattice.viewBonds();
    lattice.activateBond(Bond({0, 0}, {0, 1}, lattice.length()).ID());
//    _lattice.activate_bond(Bond({1,0},{0,1}, _lattice.length()));
    lattice.activateBond(Bond({0, 1}, {0, 2}, lattice.length()).ID());
    lattice.activateBond(Bond({0, 2}, {0, 1}, lattice.length()).ID());
    lattice.activateBond(Bond({0, 4}, {0, 0}, lattice.length()).ID());//  these are neighbor but they are different in order
    lattice.activateBond(Bond({1, 0}, {2, 0}, lattice.length()).ID());
//    _lattice.activateAllBond();
    lattice.view_bonds();
}

