//
// Created by shahnoor on 6/28/21.
//

#include "percolation_rsbd_v13.h"
#include "percolation_v13.h"
#include "status.h"
#include "../index/index.h"

using namespace std;

/**
 * return :
 *  0  -> successfully chosen an empty site
 *  1  -> sites are remaining but current site is not empty
 *  -1 -> no remaining empty sites
 * @return
 */
P_STATUS SitePercolationL1_v13::select_site() {
//    cout << "Entry SitePercolationL1_v13::select_site " << endl;

    if (current_idx >= lattice_ref.get_site_count()) {
//        cout << "No sites to occupy" << endl;
        return P_STATUS::EMPTY_SITE_LIST;
    }
    value_type rnd = random_int(current_idx, site_ids_indices.size());
#ifdef UNIT_TEST
    if(rnd >= site_ids_indices.size()){
        cout << "index out of range " << __LINE__ << endl;
        cout << rnd << endl;
        cout << site_ids_indices.size() << endl;
    }
#endif
    auto central_X = site_ids_indices[rnd];
//    cout << "central_X " << central_X << endl;

//    int gid = lattice_ref.get_site_by_id(central_X).get_gid();
//    cout << "gid " << gid << endl;
    if(lattice_ref.get_site_by_id(central_X).is_occupied()){
//        viewLattice(0);
//        viewLattice(1);
//        viewLattice(2);
//        cout << "X is occupied" << endl;
//        exit(-1);
        ++x_occupied;

        auto sites = lattice_ref.get_all_neighbor_sites(central_X);


        auto central2 = sites[_random_engine()% sites.size()];

        if (lattice_ref.get_site_by_id(central2).is_occupied()){
            lattice_ref.get_site_by_id(central_X).reduce_1st_nn();// remove
            if (lattice_ref.get_site_by_id(central_X).is_removable(1)){ // remove
//            if(lattice_ref.isRemovable(central_X, 1)){
//# print("is_removable")
//# print("self.site_ids_indices before ", self.site_ids_indices)
//# print("rnd ", rnd, " self.current_idx ", self.current_idx)

                site_ids_indices[rnd] = site_ids_indices[current_idx];

                current_idx += 1;
            }
            return P_STATUS::CURRENT_SITE_NOT_EMPTY;
        }
        else {
            selected_id = central2;
        }

    }
    else {
        selected_id = central_X;
    }


//    current_site = lattice_ref.get_site_by_id(selected_id);
//# print("selected id ", self.selected_id)
//    occupied_site_count += 1;

//    current_site = lattice_ref.get_site_by_id(selected_id);
//    cout << "selected id " << selected_id << " site " << lattice_ref.get_site_by_id(selected_id).get_str() << endl;
//    cout << "selected id " << selected_id << " site " << get_current_site().get_str() << endl;
    return P_STATUS::SUCESS;





}

SitePercolationL1_v13::SitePercolationL1_v13(int length, value_type seed, bool generate_seed)
        :SitePercolationL0_v13(length, seed, generate_seed) {
//    setRandomState(seed, generate_seed);
    cout << "class : SitePercolationL1_v13" << endl;
}

P_STATUS SitePercolationL2_v13::select_site() {
    //    cout << "Entry SitePercolationL1_v13::select_site " << endl;

    if (current_idx >= lattice_ref.get_site_count()) {
//        cout << "No sites to occupy" << endl;
        return P_STATUS::EMPTY_SITE_LIST;
    }
    value_type rnd = random_int(current_idx, site_ids_indices.size());
#ifdef UNIT_TEST
    if(rnd >= site_ids_indices.size()){
        cout << "index out of range " << __LINE__ << endl;
        cout << rnd << endl;
        cout << site_ids_indices.size() << endl;
    }
#endif
    auto central_X = site_ids_indices[rnd];
//    cout << "central_X " << central_X << endl;

//    int gid = lattice_ref.get_site_by_id(central_X).get_gid();
//    cout << "gid " << gid << endl;
    if(lattice_ref.get_site_by_id(central_X).is_occupied()){
//        viewLattice(0);
//        viewLattice(1);
//        viewLattice(2);
//        cout << "X is occupied" << endl;
//        exit(-1);
//        ++x_occupied;

        auto sites = lattice_ref.get_all_neighbor_sites(central_X);


        auto nn_1_id = sites[_random_engine()% sites.size()];

        if (lattice_ref.get_site_by_id(nn_1_id).is_occupied()){
            // get site id in the direction of 1st nn
            auto nn_2_id = get_2nn_in_1nn_direction(central_X, nn_1_id);

            // if it's empty then select it

            // otherwise occupy it
            if (lattice_ref.get_site_by_id(nn_2_id).is_occupied()){
                lattice_ref.get_site_by_id(central_X).reduce_2nd_nn();
                if (lattice_ref.get_site_by_id(central_X).is_removable(2)){
//            if(lattice_ref.isRemovable(central_X, 1)){
//# print("is_removable")
//# print("self.site_ids_indices before ", self.site_ids_indices)
//# print("rnd ", rnd, " self.current_idx ", self.current_idx)

                    site_ids_indices[rnd] = site_ids_indices[current_idx];

                    current_idx += 1;

                }
                return P_STATUS::CURRENT_SITE_NOT_EMPTY;
            } else{
                selected_id = nn_2_id;
            }

        }
        else {
            selected_id = nn_1_id;
        }

    }
    else {
        selected_id = central_X;
    }


//    current_site = lattice_ref.get_site_by_id(selected_id);
//# print("selected id ", self.selected_id)
//    occupied_site_count += 1;

//    current_site = lattice_ref.get_site_by_id(selected_id);
//    cout << "selected id " << selected_id << " site " << lattice_ref.get_site_by_id(selected_id).get_str() << endl;
//    cout << "selected id " << selected_id << " site " << get_current_site().get_str() << endl;
    return P_STATUS::SUCESS;


}

SitePercolationL2_v13::SitePercolationL2_v13(int length, value_type seed, bool generate_seed)
        :SitePercolationL0_v13(length, seed, generate_seed) {
    cout << "class : SitePercolationL2_v13" << endl;
}

int SitePercolationL2_v13::get_2nn_in_1nn_direction(int central_X_id, int nn_1_id) {

    auto center = lattice_ref.get_site_by_id(central_X_id).get_index(_length);
    auto nn_1 = lattice_ref.get_site_by_id(nn_1_id).get_index(_length);

    int delta_c = nn_1.col() - center.col();
    int delta_r = nn_1.row() - center.row();
    if (delta_c == 0 && delta_r == 0){
        cout << "Both indices are same : line " << __LINE__ << endl;
    }
    else if(delta_c > 1 || delta_r > 1){
        // meaning, the site_index_sequence are on the opposite edges
//        cout << "2nd index is not the First nearest neighbor : line " << __LINE__ << " : file " << __FILE__ << endl;
    }
    int rr = (nn_1.row() + delta_r + _length) % _length;
    int cc = (nn_1.col() + delta_c + _length) % _length;
    return rr*_length + cc;

}
