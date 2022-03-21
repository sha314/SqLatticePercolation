//
// Created by shahnoor on 2/15/21.
//

#include "bond_v13.h"

using namespace std;

std::vector<int> Bond_v14::get_connected_sites(int length) {
    cout << "connected sites of bond " << get_id() << " is ";
    int id = get_id();
    int v_bond = id % (2*length*length);  // 1 if vertical bond
    int site_0=-1, site_1=-1;
    if (v_bond == 0){
        // horizontal
        site_0 = id;
        site_1 = id + 1;
    }else{
        site_0 = id - (length*length);
        site_1 = site_0 + length;
    }

    cout << site_0 << ", " << site_1 << endl;

    return {site_0, site_1};
}
