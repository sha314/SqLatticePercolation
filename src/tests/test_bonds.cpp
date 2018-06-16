//
// Created by shahnoor on 10/27/17.
//

#include "test_bonds.h"
#include "../lattice/bond.h"

using namespace std;

void TEST_Bonds::test_bond(){
    Bond bond0({0,0}, {0,3}, 5);
    Bond bond1({0,3}, {0,0}, 5);
    Bond bond2({3,0}, {0,0}, 5);
    Bond bond3({0,0}, {3,0}, 5);
    bond0.activate();
    bond1.activate();
    bond2.activate();
    bond3.activate();

    std::cout << "bond " << bond0.id() << " : " << bond0 << endl;
    std::cout << "bond " << bond1.id() << " : " << bond1 << endl;
    std::cout << "bond " << bond2.id() << " : " << bond2 << endl;
    std::cout << "bond " << bond3.id() << " : " << bond3 << endl;
}



