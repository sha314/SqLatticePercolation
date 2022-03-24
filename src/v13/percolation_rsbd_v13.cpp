//
// Created by shahnoor on 6/28/21.
//

#include "percolation_rsbd_v13.h"
#include "percolation_v13.h"
#include "status.h"

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

SitePercolationL1_v13::SitePercolationL1_v13(int length, value_type seed, bool generate_seed):SitePercolationL0_v13(length, seed, generate_seed) {
//    setRandomState(seed, generate_seed);
    cout << "class : SitePercolationL1_v13" << endl;
}

