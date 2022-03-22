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
    int row = get_id() / length;
    int col = get_id() % length;

    int site_count = length*length;
    int bond_count = site_count*2;

    int right_bond = get_id();
    int bottom_bond = right_bond + site_count;

    // left_bond is the right_bond of left_site
    int left_site_id = get_id() - 1;
    if (get_id()%length == 0){ // for first column elements
        left_site_id += length;
    }
    int left_bond = left_site_id;

    // top_bond is the bottom_bond of top_site
    int top_site_id = get_id() - length;
    if (get_id() < length){
        // top row elements
//        top_site_id = get_id() - length + site_count;
        top_site_id += site_count;
    }
//# top_bond = (s0_index+self.bond_count-self.length)%self.site_count + self.site_count
    int top_bond = top_site_id + site_count;
    return {right_bond, bottom_bond, left_bond, top_bond};

}

std::vector<int> Site_v14::get_connecting_bonds_v2(int length) {
    int row = get_id() / length;
    int col = get_id() % length;

    int site_count = length*length;
    int bond_count = site_count*2;

    int right_bond = get_id();
    int bottom_bond = right_bond + site_count;

    // left_bond is the right_bond of left_site
    int CX_L = (length + col - 1) % length;
    int left_bond = CX_L + length*row;

    // top_bond is the bottom_bond of top_site
    int RX_T = (length + row - 1) % length;
    int top_bond = col + length*RX_T + site_count;
    return {right_bond, bottom_bond, left_bond, top_bond};

}

std::vector<int> Site_v14::get_neighbor_sites(int length) {
    int row = get_id() / length;
    int col = get_id() % length;
    int CX_R = (col + 1) % length; // column of right site
    int CX_L = (length + col - 1) % length; // column of left site
    int RX_B = (row + 1) % length; // row of bottom site
    int RX_T = (length + row - 1) % length; // row of top site

    return {CX_R + length*row,
            CX_L + length*row,
            col + length*RX_B,
            col + length*RX_T};

}

//std::vector<int> Site_v14::get_connecting_bonds(int length) {
//    int right_bond = right_bond_of_site(length);
////# bottom_bond = s0_index + self.site_count
//    int bottom_bond = bottom_bond_of_site(length);
////# left_bond = (s0_index + self.length - 1) % self.length
//    int left_bond = left_bond_of_site(length);
////# top_bond = (s0_index+self.bond_count-self.length)%self.site_count + self.site_count
//    int top_bond = top_bond_of_site(length);
//    return {right_bond, bottom_bond, left_bond, top_bond};
//
//}


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
    int x = bottom_bond_of_site(of_bottom_site);
    cout << "top bond of site " << get_id() << " is " << x << endl;
    return x;
}

int Site_v14::right_bond_of_site(int length) {
    return get_id();
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
