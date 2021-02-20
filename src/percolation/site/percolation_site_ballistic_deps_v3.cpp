//
// Created by shahnoor on 8 Sept 2018
//

#include <cstdlib>
#include <climits>
#include <iterator>

#include "../percolation.h"
#include "percolation_site_v9.h"
#include "percolation_site_explosive.h"
#include "percolation_site_ballistic.h"

using namespace std;

/**************************
 * Ballistic Deposition
 */

/**
 *
 * @param length
 */
SitePercolationBallisticDeposition_v3::SitePercolationBallisticDeposition_v3(value_type length, bool periodicity)
        : SitePercolation_ps_v10(length, periodicity )
{

    std::cout << "Constructing SitePercolationBallisticDeposition_v3 object : line " << __LINE__ << endl;

    initialize_indices();
    indices_tmp = indices;
//    randomize_index();
}

/**
 *
 */
void SitePercolationBallisticDeposition_v3::reset() {
    SitePercolation_ps_v10::reset();
    indices_tmp = indices;
    _search_position = 0;
}

/**
 * Called only once when the object is constructed for the first time
 */
void SitePercolationBallisticDeposition_v3::initialize_indices() {
    indices = vector<value_type>(maxSites());
    for(value_type i{}; i != indices.size(); ++i){
        indices[i] = i; // assign index first
    }
}


/***************************************
 * occupy upto 1st nearset neighbor.
 * */


/*******************************************
 * Site selection methods
 */

/**
 *
 * @param sites
 * @param bonds
 * @return
 */
Index SitePercolationBallisticDeposition_v3::select_site(vector<Index> &sites, vector<BondIndex> &bonds) {
    cout << "TODO : select_site __LINE__ " << __LINE__ << endl;
    // randomly choose a site
    value_type r = std::rand() % (indices_tmp.size());

//    Index current_site = site_indices[indices_tmp[r]];
    Index current_site;
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
 * Instead of erasing the indices, place them to an inaccesible area by swapping or replacing.
 * If erasing is not performed then the performance will be better.
 *
 * L    time(v1)    time(v2)
 * 200  1.12 sec    0.75 sec
 * 400  7.8 sec     3.4 sec
 *
 * @param sites
 * @param bonds
 * @return
 */
Index SitePercolationBallisticDeposition_v3::select_site_upto_1nn_v2(
        vector<Index> &sites, vector<BondIndex> &bonds
) {
    // randomly choose a site
    uniform_int_distribution<size_t> dist(_search_position, randomized_index.size());
//    value_type  r = _search_position + (_random() % (randomized_index.size() - _search_position));
    value_type  r = dist(_random_engine);
    value_type index = randomized_index[r];
    Index current_site = site_indices[index];
//    cout << "current site " << current_site << endl;
    // find the bonds for this site

//    connection_v1(current_site, site_index_sequence, bonds);
//    connection_v2(current_site, sites, bonds);
    _lattice.get_neighbors(current_site, sites, bonds);
//    cout << "current site " << current_site << " active ? " << lattice_ref.getSite(current_site).isActive() << endl;
    if (_lattice.getSite(current_site).isActive()){ // if the current site is occupied or active
#ifdef DEBUG_FLAG
        cout << "1st neighbors of " << current_site << ":{";
        copy(sites.begin(), sites.end(), ostream_iterator<Index>(cout, ","));
        cout << "}";
#endif
        value_type r2 = _random_engine() % (sites.size());
        current_site = sites[r2]; // select one of the neighbor randomly

        if(_lattice.getSite(current_site).isActive()){
            // if the neighbor is also occupied cancel current step
            bool flag = true;
#ifdef DEBUG_FLAG
            cout << "if one of the neighbor is inactive. it's engouh to go on" << endl;
#endif
            for(auto s : sites){
#ifdef DEBUG_FLAG
                cout << s << "->";
#endif
                if(!_lattice.getSite(s).isActive()){
                    // if one of the neighber is unoccupied then
                    flag = false;
#ifdef DEBUG_FLAG
                    cout << " inactive" << endl;
#endif
                    break;
                }
#ifdef DEBUG_FLAG
                cout << " active"<< endl;
#endif
            }
            if(flag){
                /*
                 * Instead of erasing the indices, place them to an inaccesible area by swapping or replacing.
                 * If erasing is not performed then the performance will be better.
                 */
                randomized_index[r] = randomized_index[_search_position]; // replace r-th value with the unused value
                ++_search_position;// increment search position so that this value is never selected again

                throw OccupiedNeighbor{"all of the 1nd neighbors are occupied : line " + std::to_string(__LINE__)};
            }
            throw OccupiedNeighbor{"selected 1st neighbor is occupied : line " + std::to_string(__LINE__)};
        }
    }

//    cout << "choosing " << current_site << " out of the neighbors : line " << __LINE__ << endl;
    sites.clear();
    bonds.clear();
//    connection_v2(current_site, sites, bonds);
    _lattice.get_neighbors(current_site, sites, bonds); // since neighbor is changed now
    return current_site;
}



/**
 * Select neighbor upto 2nd nearest neighbor
 * uses direcion of motion when selecting 2nd nearest neighbor
 * @param r : index of sites in the randomized array
 *
 * L    time(v1)    time(v2)
 * 400  8.4 sec     3.4 sec
 *
 *
 * @param sites
 * @param bonds
 * @return
 */
Index SitePercolationBallisticDeposition_v3::select_site_upto_2nn_v2(
        vector<Index> &sites, vector<BondIndex> &bonds
){
// randomly choose a site
    value_type  r = _search_position + (_random_engine() % (randomized_index.size() - _search_position));
    value_type index = randomized_index[r];
    Index central_site = site_indices[index];
    Index selected_site{};
//    cout << "current site " << current_site << endl;
    // find the bonds for this site

//    connection_v1(current_site, site_index_sequence, bonds);
//    connection_v2(central_site, sites, bonds);
    _lattice.get_neighbors(central_site, sites, bonds);
//    cout << "current site " << central_site << " active ? " << lattice_ref.getSite(central_site).isActive() << endl;
    if (_lattice.getSite(central_site).isActive()){ // if the current site is occupied or active
//        cout << "1st neighbors :{";
//        copy(sites.begin(), sites.end(), ostream_iterator<Index>(cout, ","));
//        cout << "}";

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

        value_type r2 = _random_engine() % (sites.size());
        Index nn1 = sites[r2]; // select one of the neighbor randomly
//        cout << "nn1 " << nn1 << " : line " << __LINE__ <<endl;
        Index nn2;
        if(_lattice.getSite(nn1).isActive()){
            // if the neighbor is also occupied then choose the 2nd nearest neighbor in the direction of motion
            nn2 = get_2nn_in_1nn_direction(central_site, nn1, length());

//            cout << "nn2 " << nn2 << " : line " << __LINE__ <<endl;
            // if it is also occupied the skip the step
            if(_lattice.getSite(nn2).isActive()) {
                flag_nn2 = true;

                vector<Index> nn2_sites = get_all_2nn_in_1nn_s_direction(central_site, sites, length());
                for(auto x: nn2_sites){
                    if(!_lattice.getSite(x).isActive()){
                        flag_nn2 = false;
//                        cout << "inactive";
                        break;
                    }
                }

                if(flag_nn1 && flag_nn2){
                    // since its 1st nearest neighbors are occupied
                    // and 2nd nearest neighbors are also occupied
                /*
                 * Instead of erasing the indices, place them to an inaccesible area by swapping or replacing.
                 * If erasing is not performed then the performance will be better.
                 */
                    randomized_index[r] = randomized_index[_search_position]; // replace r-th value with the unused value
                    ++_search_position;// increment search position so that this value is never selected again
                    throw OccupiedNeighbor{"all 1st and 2nd neighbors are occupied : line " + std::to_string(__LINE__)};
                }

                throw OccupiedNeighbor{"2nd neighbor is also occupied : line " + std::to_string(__LINE__)};
            }else{
                selected_site = nn2;
            }
        }else {
            selected_site = nn1;
        }
//        cout << "choosing " << selected_site << " out of the neighbors : line " << __LINE__ << endl;
//        connection_v2(selected_site, sites, bonds);
        _lattice.get_neighbors(selected_site, sites, bonds);
    }else{
        selected_site = central_site;
    }
    return selected_site;
}



/********************************************************
 * SitePercolationBallisticDeposition_v3
 * select upto 1st nearest neighbor
 */

/**
 *
 * @return
 */
bool SitePercolationBallisticDeposition_v3::occupy() {
    // if no site is available then return false

    if(_number_of_occupied_sites == maxSites()){
        return false;
    }

    try {

        value_type v = placeSite_1nn_v2();
//        _occuption_probability = occupationProbability(); // for super class


        return v != ULLONG_MAX;
    }catch (OccupiedNeighbor& on){
//        on.what();
//        cout << "line : " << __LINE__ << endl;
        return false;
    }

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
value_type SitePercolationBallisticDeposition_v3::placeSite_1nn_v2() {

    vector<BondIndex> bonds;
    vector<Index>     sites;

//    _last_placed_site = select_site_upto_1nn(sites, bonds);
    auto current_site = select_site_upto_1nn_v2(sites, bonds);

#ifdef DEBUG_FLAG
    cout << "Neighbors of " << current_site << endl;
    cout << "{" ;
    for(auto a: sites){
        cout << a << ",";
    }
    cout << "}" << endl;
    cout << "{" ;
    for(auto b: bonds){
        cout << b << ",";
    }
    cout << "}" << endl;
#endif
    return placeSite(current_site, sites, bonds);
}




/**
 *
 * @return
 */
value_type SitePercolationBallisticDeposition_v3::placeSite_2nn_v1() {
    vector<BondIndex> bonds;
    vector<Index>     sites;

    try {
//        _last_placed_site = select_site_upto_2nn(sites, bonds);
        auto current_site = select_site_upto_2nn_v2(sites, bonds);
        return placeSite(current_site, sites, bonds);
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

void SitePercolationBallisticDeposition_v3::viewRemainingSites() {
    cout << "viewRemainingSites {";
    for(size_t i{_search_position}; i < randomized_index.size(); ++i){
        cout << site_indices[randomized_index[i]] << ",";
    }

    cout << "}" << endl;
}

