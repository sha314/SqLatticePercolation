//
// Created by shahnoor on 2/15/21.
//

#include <assert.h>
#include "site_v13.h"
#include "lattice_v13.h"

using namespace std;


Index_v13 Site_v14::get_index() {



    return index;

}

Index_v13 Site_v14::get_index_v14(int length) {

    int row = get_id() / length;
    int col = get_id() % length;
#ifdef UNIT_TEST
    if (index.row() != row){
        std::cout << "row not equal" << std::endl;
    }
    if (index.col() != col){
        std::cout << "column not equal" << std::endl;
    }


#endif

    return {row, col};

}



std::vector<int> Site_v14::get_connecting_bonds(int length) {
    int right_bond = right_bond_of_site(length);
//# bottom_bond = s0_index + self.site_count
    int bottom_bond = bottom_bond_of_site(length);
//# left_bond = (s0_index + self.length - 1) % self.length
    int left_bond = left_bond_of_site(length);
//# top_bond = (s0_index+self.bond_count-self.length)%self.site_count + self.site_count
    int top_bond = top_bond_of_site(length);
    return {right_bond, bottom_bond, left_bond, top_bond};

}


int Site_v14::bottom_bond_of_site(int length) {
    return get_id() + length*length;
}

/**
 * top index of a site 'm' is the bottom index of some other site 'n'.
            n = (m + 2*L*L - L) (mod L*L)

 * @param s0_index
 * @return
 */
int Site_v14::top_bond_of_site(int length) {
    int site_count = length*length;
    int bond_count = site_count*2;
    int of_bottom_site = (get_id()+bond_count-length)%site_count;
    return bottom_bond_of_site(of_bottom_site);
}

int Site_v14::right_bond_of_site(int s0_index) {
    return s0_index;
}

/***
left index of a site 'm' is the right index of some other site 'n'.
n =
*/
int Site_v14::left_bond_of_site(int length) {

    int row = get_id() / length;
    int col = get_id() % length;

    int left_site = row*length + (col + length - 1) % length;
//# print("left of ", s0_index, " is the right of ", left_site)
    return right_bond_of_site(left_site);
}
