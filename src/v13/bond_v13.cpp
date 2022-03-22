//
// Created by shahnoor on 2/15/21.
//

#include "bond_v13.h"

using namespace std;

std::vector<int> Bond_v14::get_connected_sites(int length) {

    int id = get_id();
    int row = get_id() / length;
    int col = get_id() % length;
    int CX_R = (col + 1) % length; // column of right site

    int site_count = length * length;
    int v_bond = id / site_count;  // 1 if vertical bond
//    cout << "v bond " << v_bond << endl;
    int site_0=-1, site_1=-1;
    // horizontal
    site_0 = id;
    site_1 = CX_R + row*length;
    if (v_bond == 1){
//        cout << "Vertical bond" << endl;
        site_0 = id - site_count;
        site_1 = (site_0 + length) % site_count;
    }


#ifdef DEBUG_FLAG
    cout << "connected sites of bond " << get_id() << " is ";
    cout << site_0 << ", " << site_1 << endl;
    if (site_0 < 0 || site_1 << 0){
        cout << "They cannot be negative" << endl;
    }
#endif

    return {site_0, site_1};
}

int Bond_v14::get_row(int length) {
    return get_id() / length;
}

int Bond_v14::get_type(int length) {
    return get_id() / (length*length);
}

int Bond_v14::get_col(int length) {
    return get_id() % length;
}


