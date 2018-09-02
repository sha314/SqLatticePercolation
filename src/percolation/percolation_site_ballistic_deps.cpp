//
// Created by shahnoor on 31-Mar-18.
//

#include <cstdlib>
#include <climits>

#include "percolation.h"

/**************************
 * Ballistic Deposition
 */

/**
 *
 * @param length
 */
SitePercolationBallisticDeposition::SitePercolationBallisticDeposition(value_type length, bool periodicity)
        : SitePercolation_ps_v8(length, periodicity )
{

    std::cout << "Constructing SitePercolationBallisticDeposition object : line " << __LINE__ << endl;

    initialize_indices();
    indices_tmp = indices;
//    randomize_index();
}

/**
 *
 */
void SitePercolationBallisticDeposition::reset() {
    SitePercolation_ps_v8::reset();
    indices_tmp = indices;
}

/**
 * Called only once when the object is constructed for the first time
 */
void SitePercolationBallisticDeposition::initialize_indices() {
    indices = vector<value_type>(maxSites());
    for(value_type i{}; i != indices.size(); ++i){
        indices[i] = i; // assign index first
    }
}


/***************************************
 * occupy upto 1st nearset neighbor.
 * */

value_type SitePercolationBallisticDeposition::placeSite_nn_v0(int n) {
    cout << "Under Devlopment : line " << __LINE__ << endl;
    vector<BondIndex> bonds;
    vector<Index>     sites;

    try {
        if (n == 1) { // select upto 1st nearest neighbor
            _last_placed_site = select_site_upto_1nn(sites, bonds);

        } else if (n == 2) { // select upto 2nd nearest neighbor
            _last_placed_site = select_site_upto_2nn(sites, bonds);
        } else {
            _last_placed_site = select_site(sites, bonds);
        }
    }catch (OccupiedNeighbor& on){
        on.what();
        cout << "line : " << __LINE__ << endl;
        throw; // rethrowing
    }

    // now use the selected site to occupy the lattice

    _lattice.activate_site(_last_placed_site);

    ++_number_of_occupied_sites;
//    ++_index_sequence_position;

    // find one of hv_bonds in _clusters and add ever other value to that place. then erase other position
    set<value_type> found_index_set = find_index_for_placing_new_bonds_v3(sites);

//    cout << "Found indices " << found_index_set << endl;

    subtract_entropy_for_bond(found_index_set);  // tracking entropy change
    value_type merged_cluster_index = manage_clusters_v10(
            found_index_set, bonds, _last_placed_site
    );
    add_entropy_for_bond(merged_cluster_index); // tracking entropy change

    // running tracker
    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster

    return merged_cluster_index;
}


/*******************************************
 * Site selection methods
 */

/**
 *
 * @param sites
 * @param bonds
 * @return
 */
Index SitePercolationBallisticDeposition::select_site(vector<Index> &sites, vector<BondIndex> &bonds) {
    // randomly choose a site
    value_type r = std::rand() % (indices_tmp.size());

    Index current_site = index_sequence[indices_tmp[r]];
    cout << "current site " << current_site << endl;
    // find the bonds for this site

    if (_lattice.getSite(current_site).isActive()){
        indices_tmp.erase(indices_tmp.begin()+r);

        throw OccupiedNeighbor{"all of the 1nd neighbors are occupied : line " + std::to_string(__LINE__)};
    }

    cout << "choosing " << current_site << " out of the neighbors : line " << __LINE__ << endl;
    sites.clear();
    bonds.clear();
    connection_v2(current_site, sites, bonds);
    return current_site;
}

/**
 *
 * @param sites
 * @param bonds
 * @return
 */
Index SitePercolationBallisticDeposition::select_site_upto_1nn(
        vector<Index> &sites, vector<BondIndex> &bonds
) {
    // randomly choose a site
    value_type r = std::rand() % (indices_tmp.size());

    Index current_site = index_sequence[indices_tmp[r]];
//    cout << "current site " << current_site << endl;
    // find the bonds for this site

//    connection_v1(current_site, sites, bonds);
    connection_v2(current_site, sites, bonds);

    if (_lattice.getSite(current_site).isActive()){ // if the current site is occupied or active
        value_type r2 = std::rand() % (sites.size());
        current_site = sites[r2]; // select one of the neighbor randomly

        if(_lattice.getSite(current_site).isActive()){
            // if the neighbor is also occupied cancel current step
            bool flag = true;
//            cout << "if one of the neighbor is inactive. it's engouh to go on" << endl;
            for(auto s : sites){
//                cout << s << "->";
                if(!_lattice.getSite(s).isActive()){
                    // if one of the neighber is unoccupied then
                    flag = false;
//                    cout << " inactive" << endl;
                    break;
                }
//                cout << " active"<< endl;
            }
            if(flag){
                // erase the index, since its four neighbors are occupied
                indices_tmp.erase(indices_tmp.begin()+r);
                throw OccupiedNeighbor{"all of the 1nd neighbors are occupied : line " + std::to_string(__LINE__)};
            }
            throw OccupiedNeighbor{"selected 1st neighbor is occupied : line " + std::to_string(__LINE__)};
        }
    }

//    cout << "choosing " << current_site << " out of the neighbors : line " << __LINE__ << endl;
    sites.clear();
    bonds.clear();
    connection_v2(current_site, sites, bonds);
    return current_site;
}



/**
 * Select neighbor upto 2nd nearest neighbor
 * uses direcion of motion when selecting 2nd nearest neighbor
 * @param r : index of sites in the randomized array
 * @param sites
 * @param bonds
 * @return
 */
Index SitePercolationBallisticDeposition::select_site_upto_2nn(
        vector<Index> &sites, vector<BondIndex> &bonds
){
    value_type r = std::rand() % (indices_tmp.size());

    Index central_site = index_sequence[indices_tmp[r]];
    Index selected_site;
//    cout << "current site " << central_site << endl;
    // find the bonds for this site


//    connection_v1(central_site, sites, bonds);
    connection_v2(central_site, sites, bonds);

    if (_lattice.getSite(central_site).isActive()){
//        cout << "selected " << central_site << " : line " << __LINE__ << endl;

        bool flag_nn1 = true; // true means all 1st nearest neighbors are occupied
        bool flag_nn2 = true; // true means all 2nd nearest neighbors are occupied
//        cout << "if one of the neighbor is inactive. it's engouh to go on" << endl;
        for(auto s : sites){
//            cout << s << "->";
            if(!_lattice.getSite(s).isActive()){
                // if one of the neighber is unoccupied then
                flag_nn1 = false;
//                cout << " inactive" << endl;
                break;
            }
//            cout << " active"<< endl;
        }

        value_type r2 = std::rand() % (sites.size());
        Index nn1 = sites[r2]; // select one of the neighbor randomly
//        cout << "nn1 " << nn1 << " : line " << __LINE__ <<endl;
        Index nn2;
        if(_lattice.getSite(nn1).isActive()){
            // if the neighbor is also occupied then choose the 2nd nearest neighbor in the direction of motion
            nn2 = get_2nn_in_1nn_direction(central_site, nn1, length());
            if(!_periodicity){
                // if periodic boundary condition is not enabled then sites on the opposite edges will not contribute
                vector<Index> tmp_sites;
                vector<BondIndex> tmp_bonds;
                // will find all possible neighbors of the selected first nearest neighbor
                connection_v2(nn1, tmp_sites, tmp_bonds);
                bool valid{false};
                for(auto s: tmp_sites){
                    if(nn2 == s){
//                        cout << "valid 2nd nearest neighbor : line " << __LINE__ << endl;
                        valid = true;
                        break;
                    }
                }
                if(!valid){
                    throw InvalidNeighbor{"invalid 2nd nearest neighbor : line " + std::to_string(__LINE__)};
                }
            }
//            cout << "nn2 " << nn2 << " : line " << __LINE__ <<endl;
            // if it is also occupied the skip the step
            if(_lattice.getSite(nn2).isActive()) {
                flag_nn2 = true;

                vector<Index> nn2_sites = get_2nn_s_in_1nn_s_direction(central_site, sites, length());
                for(auto x: nn2_sites){
                    if(!_lattice.getSite(x).isActive()){
                        flag_nn2 = false;
//                        cout << "inactive";
                        break;
                    }
                }

                if(flag_nn1 && flag_nn2){
                    // erase the index, since its 1st nearest neighbors are occupied
                    // and 2nd nearest neighbors are also occupied
//                    cout << "ignoring this step : line " << __LINE__ << endl;
                    indices_tmp.erase(indices_tmp.begin()+r);
                }

                throw OccupiedNeighbor{"2nd neighbor is also occupied : line " + std::to_string(__LINE__)};
            }else{
                selected_site = nn2;
            }
        }else {
            selected_site = nn1;
        }

//        cout << "choosing " << selected_site << " out of the neighbors : line " << __LINE__ << endl;
        sites.clear();
        bonds.clear();

        connection_v2(selected_site, sites, bonds);
    }else{
        selected_site = central_site;
    }
    return selected_site;
}





///**
// *
// */
//void SitePercolationBallisticDeposition::randomize_index(){
//
//    indices_random = vector<value_type>(_length_squared);
//
//    value_type r{};
//    for(value_type i{}; i != indices_random.size(); ++i){
//        r = i + std::rand() % (_length_squared - i);
//        indices_random[i] = indices[r];
//        indices_random[r] = indices[i];
//    }
//}



/********************************************************
 * SitePercolationBallisticDeposition
 * select upto 1st nearest neighbor
 */

/**
 *
 * @return
 */
bool SitePercolationBallisticDeposition::occupy() {
    // if no site is available then return false

    if(_number_of_occupied_sites == maxSites()){
        return false;
    }

    try {
//        value_type v = placeSite_1nn_v0(); // debugging version
        value_type v = placeSite_1nn_v1();
        _occuption_probability = occupationProbability(); // for super class


        return v != ULLONG_MAX;
    }catch (OccupiedNeighbor& on){
//        on.what();
//        cout << "line : " << __LINE__ << endl;
        return false;
    }

}


/**
* debugging version
*
* 1. Randomly select a site from all sites
* 2. If it is not occupied occupy it.
* 3. If it is occupied select one of the 4 neighbor to occupy
* 4. If the selected neighbor is also occupied cancel current step
* 4. form cluster and track all informations
* 5. go to step 1
* 6. untill spanning cluster appears or no unoccupied site
*/
value_type SitePercolationBallisticDeposition::placeSite_1nn_v0() {
    vector<BondIndex> bonds;
    vector<Index>     sites;

    _last_placed_site = select_site_upto_1nn(sites, bonds);

    _lattice.activate_site(_last_placed_site);

    ++_number_of_occupied_sites;
//    ++_index_sequence_position;

    // find one of hv_bonds in _clusters and add ever other value to that place. then erase other position
    set<value_type> found_index_set = find_index_for_placing_new_bonds_v3(sites);

//    cout << "Found indices " << found_index_set << endl;

    subtract_entropy_for_bond(found_index_set);  // tracking entropy change
    value_type merged_cluster_index = manage_clusters_v10(
            found_index_set, bonds, _last_placed_site
    );
    add_entropy_for_bond(merged_cluster_index); // tracking entropy change

    // running tracker
    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster

    return merged_cluster_index;
}


/**
 *
 * 1. Randomly select a site from all sites
 * 2. If it is not occupied occupy it.
 * 3. If it is occupied select one of the 4 neighbor to occupy
 * 4. If the selected neighbor is also occupied cancel current step
 * 4. form cluster and track all informations
 * 5. go to step 1
 * 6. untill spanning cluster appears or no unoccupied site
 */
value_type SitePercolationBallisticDeposition::placeSite_1nn_v1() {

    vector<BondIndex> bonds;
    vector<Index>     sites;

    // randomly choose a site
    value_type r = rand() % (indices_tmp.size());

    Index current_site = index_sequence[indices_tmp[r]];
    // find the bonds for this site

//    connection_v1(current_site, sites, bonds);
    connection_v2(current_site, sites, bonds);

    if (_lattice.getSite(current_site).isActive()){

        bool flag = true;
        for(auto s : sites){
            if(!_lattice.getSite(s).isActive()){
                // if one of the neighber is unoccupied then
                flag = false;
                break;
            }
        }
        if(flag){
            // erase the index, since its four neighbors are occupied
            indices_tmp.erase(indices_tmp.begin()+r);
            return ULLONG_MAX;
        }
        value_type r2 = rand() % (sites.size());
        current_site = sites[r2]; // select one of the neighbor randomly
        if(_lattice.getSite(current_site).isActive()){
            // if the neighbor is also occupied cancel current step
            return ULONG_MAX;// unsigned long int maximum value
        }
        sites.clear();
        bonds.clear();
        connection_v2(current_site, sites, bonds);
    }

    _last_placed_site = current_site;

    _lattice.activate_site(current_site);

    ++_number_of_occupied_sites;
//    ++_index_sequence_position;

    // find one of hv_bonds in _clusters and add ever other value to that place. then erase other position
    set<value_type> found_index_set = find_index_for_placing_new_bonds_v3(sites);

//    cout << "Found indices " << found_index_set << endl;

    subtract_entropy_for_bond(found_index_set);  // tracking entropy change
    value_type merged_cluster_index = manage_clusters_v10(
            found_index_set, bonds, current_site
    );
    add_entropy_for_bond(merged_cluster_index); // tracking entropy change

    // running tracker
    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster

    return merged_cluster_index;
}



/**
 *
 * 1. Randomly select a site from all sites
 * 2. If it is not occupied occupy it.
 * 3. If it is occupied select one of the 4 neighbor to occupy
 * 4. If the selected neighbor is also occupied cancel current step
 * 4. form cluster and track all informations
 * 5. go to step 1
 * 6. untill spanning cluster appears or no unoccupied site
 */
value_type SitePercolationBallisticDeposition::placeSite_1nn_v2() {

    vector<BondIndex> bonds;
    vector<Index>     sites;

    _last_placed_site = select_site_upto_1nn(sites, bonds);

    return placeSite_v12(_last_placed_site, sites, bonds);
}




/**
 *
 * @return
 */
value_type SitePercolationBallisticDeposition::placeSite_2nn_v0() {
    vector<BondIndex> bonds;
    vector<Index>     sites;

    _last_placed_site = select_site_upto_2nn(sites, bonds);

    _lattice.activate_site(_last_placed_site);

    ++_number_of_occupied_sites;

    // find one of hv_bonds in _clusters and add ever other value to that place. then erase other position
    set<value_type> found_index_set = find_index_for_placing_new_bonds_v3(sites);

//    cout << "Found indices " << found_index_set << endl;

    subtract_entropy_for_bond(found_index_set);  // tracking entropy change
    value_type merged_cluster_index = manage_clusters_v10(
            found_index_set, bonds, _last_placed_site
    );
    add_entropy_for_bond(merged_cluster_index); // tracking entropy change

    // running tracker
    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster

    return merged_cluster_index;
}



/**
 *
 * @return
 */
value_type SitePercolationBallisticDeposition::placeSite_2nn_v1() {
    vector<BondIndex> bonds;
    vector<Index>     sites;

    try {
        _last_placed_site = select_site_upto_2nn(sites, bonds);
        return placeSite_v12(_last_placed_site, sites, bonds);
        //    return placeSite_v11(_last_placed_site);
    }catch (OccupiedNeighbor& e){
//        cout << "Exception !!!!!!!!!!!!!!!!!!" << endl;
//        e.what();
    }catch (InvalidNeighbor& b){
//        cout << "Exception !!!!!!!!!!!!!!!!!!" << endl;
//        b.what();
    }
    return ULONG_MAX;
}

