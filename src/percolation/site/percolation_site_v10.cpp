//
// Created by shahnoor on 01.01.2019
//


#include <cstdlib>
#include <climits>
#include <unordered_set>
#include <mutex>

#include "../percolation.h"
#include "percolation_site_v10.h"

#include "../../util/printer.h"
#include "../site_position.h"
#include "../../lattice/connector.h"
#include "../../index/delta.h"
#include <omp.h>
#include <thread>
#include <algorithm>
#include <iterator>

#include "../../util/checking.h"
#include "../../util/time_tracking.h"

using namespace std;



/**
 *
 * @param length       : length of the lattice
 * @param impure_sites : number of impure sites. cannot be greater than length*length
 */
SitePercolation_ps_v10::SitePercolation_ps_v10(value_type length, bool periodicity)
        :SqLatticePercolation(length)
{
    std::cout << "Constructing SitePercolation_ps_v10 object : line " << __LINE__ << endl;
    SqLatticePercolation::set_type('s');
    _periodicity = periodicity;
    _index_sequence_position = 0;

    min_index = 0;
    max_index = length - 1;
    _max_iteration_limit = maxSites();
//    site_indices.resize(maxSites());

    // to improve performence
    number_of_sites_to_span.reserve(maxSites());
    number_of_bonds_to_span.reserve(maxSites());

    _top_edge.reserve(length);
    _bottom_edge.reserve(length);
    _left_edge.reserve(length);
    _right_edge.reserve(length);

    initialize();
}


/**
 * collections of statements that must be called once while constructing
 * these statements does not require any random function
 */
void SitePercolation_ps_v10::initialize() {
    //
    randomized_index.resize(maxSites());
    for(size_t i{}; i < maxSites(); ++i){randomized_index[i]=i;}
    _lattice = SqLattice(length(), true, false, false, true);   // since it is a site percolation all bonds will be activated by default
    site_indices = _lattice.getSites();
    bond_indices = _lattice.getBonds();
}


/**
 * Called only once when the object is constructed for the first time
 */
void SitePercolation_ps_v10::initialize_index_sequence() {


//    value_type m{}, n{};
//    for (value_type i{}; i != index_sequence.size(); ++i) {
//        randomized_index[i] = i;
//        index_sequence[i] = Index(m, n);
//        ++n;
//        if (n == length()) {
//            n = 0;
//            ++m;
//        }
//    }
    //for (value_type i{}; i != bond_index_sequence.size(); ++i) {cout << bond_index_sequence[i] << endl;}
}

/**
 * Called only once when the object is constructed for the first time
 */
void SitePercolation_ps_v10::initialize_cluster() {

    _clusters.resize(bond_indices.size());
    for (int i{}; i != bond_indices.size(); ++i) {
        _clusters[i].addBondIndex(bond_indices[i]);
        _lattice.setGroupID(bond_indices[i], i);
        _clusters[i].set_ID(i);
    }
}


/**
 * Reset all calculated values and then call initiate()
 * to initiallize for reuse
 *
 * caution -> it does not erase _calculation_flags, for it will be used for calculation purposes
 */
void SitePercolation_ps_v10::reset() {
    SqLatticePercolation::reset();
    // variables
    _number_of_occupied_sites = 0;
    _index_sequence_position = 0;
//    _first_spanning_cluster_id = -1;
    _cluster_id = 0;
//    bonds_in_largest_cluster = 0;
//    sites_in_largest_cluster = 0;

    // containers
//    randomized_index_sequence.clear();    // reseted in the initialize function

//    _number_of_occupied_sites.clear();
//    _entropy_by_bond.clear();
    number_of_sites_to_span.clear();
    number_of_bonds_to_span.clear();
//    spanning_cluster_ids.clear();
    _spanning_sites.clear();
    _wrapping_sites.clear();
//    wrapping_cluster_ids.clear();

    _bonds_in_cluster_with_size_two_or_more = 0;
//    _id_largest_cluster = 0;

//    _id_last_modified_cluster = -1;

    _index_last_modified_cluster = 0;  // id of the last modified cluster
//    _index_largest_cluster = 0;
    _number_of_bonds_in_the_largest_cluster = 0;
    _number_of_sites_in_the_largest_cluster = 0;
//    _cluster_id_set.clear();

    // clearing edges
    _top_edge.clear();
    _bottom_edge.clear();
    _left_edge.clear();
    _right_edge.clear();

    _spanning_occured = false;

    time_relabel = 0;
    _total_relabeling = 0;

    init();
}


/**
 * Randomize the indices
 */
void SitePercolation_ps_v10::randomize_v2(){

    std::shuffle(randomized_index.begin(), randomized_index.end(), _random_engine);
//
//    cout << "randomized {";
//    for(size_t i{}; i < randomized_index.size(); ++i){
//        cout << randomized_index[i] << ",";
//    }
//    cout << "}" << endl;
}


/*************************************************
 * Calculation methods
 *
 ***********************************/

/*
 * Instead of calculating entropy for 1000s of cluster in every iteration
 * just keep track of entropy change, i.e.,
 * how much to subtract and how much to add.
 */
/**
 * Must be called before merging the clusters
 * @param found_index_set
 */
void SitePercolation_ps_v10::subtract_entropy_for_bond(const vector<BondIndex> bonds){
    double nob, mu_bond;
    long double H{};
    int id{-1};
#ifdef DEBUG_FLAG
    cout << "subtracting bond of id {" ;
#endif
    unordered_set<int> ids;
    for(auto b: bonds){
        id = _lattice.getGroupID(b);
        if(ids.count(id) > 0) continue; // preventing double counting
        ids.insert(id);
#ifdef DEBUG_FLAG
        cout << id << ",";
#endif
        nob = _clusters[id].numberOfBonds();
        mu_bond = nob / maxBonds();
        H += logl(mu_bond) * mu_bond;
    }
#ifdef DEBUG_FLAG
    cout << "}" << endl;
#endif
    _entropy += H; // since subtracting
}



/**
 * Must be called after merging the clusters
 * Cluster length is measured by bonds
 * @param index
 */
void SitePercolation_ps_v10::add_entropy_for_bond(value_type index){
    double nob = _clusters[index].numberOfBonds();
    double mu_bond = nob / maxBonds();
    long double H = logl(mu_bond) * mu_bond;
    _entropy += -H;
}



/**
 * Condition: must be called each time a site is placed
 */
void SitePercolation_ps_v10::track_numberOfBondsInLargestCluster() {

    // calculating number of bonds in the largest cluster // by cluster index
    // checking number of bonds
    if(_clusters[_index_last_modified_cluster].numberOfBonds() > _number_of_bonds_in_the_largest_cluster){
        _number_of_bonds_in_the_largest_cluster = _clusters[_index_last_modified_cluster].numberOfBonds();
    }

}

/**
 *
 */
void SitePercolation_ps_v10::track_numberOfSitesInLargestCluster(){

    // calculating number of bonds in the largest cluster // by cluster index
    // checking number of bonds
    if(_clusters[_index_last_modified_cluster].numberOfSites() > _number_of_sites_in_the_largest_cluster){
        _number_of_sites_in_the_largest_cluster = _clusters[_index_last_modified_cluster].numberOfSites();
    }
}




/**
 * When merging two cluster this function perform relabeling based
 * on the gived indices. If the index_1 is smaller then relabel
 * sites that belong to index_1 using index_2 then return index_1
 *
 * @param index_1
 * @param index_2
 * @return  : index of the cluster whose members are relabeled and this cluster must be erased
 */
value_type SitePercolation_ps_v10::relabel(value_type index_1, value_type index_2){
    if(_clusters[index_1].numberOfSites() < _clusters[index_2].numberOfSites()){
        // _clusters[index_1] is smaller

        relabel_sites(_clusters[index_1], _clusters[index_2].get_ID()); // for site percolation
//        relabel_bonds(_clusters[index_1], _clusters[index_2].set_ID()); // for bond percolation
        return index_1;
    }else{
        relabel_sites(_clusters[index_2], _clusters[index_1].get_ID()); // for site percolation
//        relabel_bonds(_clusters[index_2], _clusters[index_1].set_ID()); // for bond percolation
        return index_2;
    }
}



///**
// * Relative index of site_new with respect to root
// * @param root
// * @param site_new
// * @return
// */
//IndexRelative SitePercolation_ps_v10::getRelativeIndex(Index root, Index site_new){
////    cout << "Entry \"SitePercolation_ps_v10::getRelativeIndex\" : line " << __LINE__ << endl;
//    int delta_x = -int(root.column_) + int(site_new.column_); // if +1 then root is on the right ??
//    int delta_y = int(root.row_) - int(site_new.row_); // if +1 then root is on the top ??
//
//
//    // normalizing delta_x
//    if(delta_x > 1){
//        delta_x /= -delta_x;
//    }
//    else if(delta_x < -1){
//        delta_x /= delta_x;
//    }
//
//    // normalizing delta_y
//    if(delta_y > 1){
//        delta_y /= -delta_y;
//    }else if(delta_y < -1){
//        delta_y /= delta_y;
//    }
//
//    IndexRelative indexRelative_root = _lattice.getSite(root).relativeIndex();
////    cout << "Relative index of root " << indexRelative_root << endl;
////    cout << "Delta x,y " << delta_x << ", " << delta_y << endl;
//    IndexRelative r =  {indexRelative_root.x_ + delta_x, indexRelative_root.y_ + delta_y};
////    cout << "Relative index of site_new " << r << endl;
//    return r;
//}




/**
 * Take a bond index only if the corresponding site is active
 * lengthy but straight forward
 * @param site
 * @param neighbors
 * @param bonds
 */
void SitePercolation_ps_v10::connection_v1(Index site, vector<Index> &neighbors, vector<BondIndex> &bonds)
{
    clock_t t = clock();
    value_type prev_column  = (site.column_ + length() - 1) % length();
    value_type prev_row     = (site.row_ + length() - 1) % length();
    value_type next_row     = (site.row_ + 1) % length();
    value_type next_column  = (site.column_ + 1) % length();
    if(_periodicity){
        neighbors.resize(4);
        neighbors[0] = {site.row_, next_column};
        neighbors[1] = {site.row_, prev_column};
        neighbors[2] = {next_row, site.column_};
        neighbors[3] = {prev_row, site.column_};

        bonds.reserve(4);
        if(!_lattice.getSite(neighbors[0]).isActive()) {
            bonds.push_back({BondType::Horizontal, site.row_, site.column_});
        }
        if(!_lattice.getSite(neighbors[1]).isActive()){
            bonds.push_back({BondType::Horizontal, site.row_, prev_column});
        }
        if(!_lattice.getSite(neighbors[2]).isActive()){
            bonds.push_back({BondType::Vertical,    site.row_, site.column_});
        }
        if(!_lattice.getSite(neighbors[3]).isActive()){
            bonds.push_back({BondType::Vertical, prev_row, site.column_});
        }

    }
    else{
        // without periodicity
        if (site.row_ == min_index) { // top edge including corners
            if(site.column_ == min_index){
                // upper left corner

                neighbors.resize(2);
                neighbors[0] = {site.row_, next_column};
                neighbors[1] = {next_row, site.column_};

                bonds.reserve(2);
                if(!_lattice.getSite(neighbors[0]).isActive()){
                    bonds.push_back({BondType::Horizontal, site.row_, site.column_});
                }
                if(!_lattice.getSite(neighbors[1]).isActive()){
                    bonds.push_back({BondType::Vertical, site.row_, site.column_});
                }

            }
            else if(site.column_ == max_index){
                // upper right corner

                neighbors.resize(2);
                neighbors[0] = {site.row_, prev_column};
                neighbors[1] = {next_row, site.column_};

                bonds.reserve(2);
                if(!_lattice.getSite(neighbors[0]).isActive()){
                    bonds.push_back({BondType::Horizontal, site.row_, prev_column});
                }
                if(!_lattice.getSite(neighbors[1]).isActive()){
                    bonds.push_back({BondType::Vertical, site.row_, site.column_});
                }

            }
            else{
                // top edge excluding corners
                neighbors.resize(3);
                neighbors[0] = {site.row_, next_column};
                neighbors[1] = {site.row_, prev_column};
                neighbors[2] = {next_row, site.column_};

                bonds.reserve(4);
                if(!_lattice.getSite(neighbors[0]).isActive()) {
                    bonds.push_back({BondType::Horizontal, site.row_, site.column_});
                }
                if(!_lattice.getSite(neighbors[1]).isActive()){
                    bonds.push_back({BondType::Horizontal, site.row_, prev_column});
                }
                if(!_lattice.getSite(neighbors[2]).isActive()){
                    bonds.push_back({BondType::Vertical,    site.row_, site.column_});
                }

            }
        }
        else if (site.row_ == max_index) { // bottom edge including corners
            if (site.column_ == min_index) {
                // lower left corner
                neighbors.resize(2);
                neighbors[0] = {site.row_, next_column};
                neighbors[1] = {prev_row, site.column_};

                bonds.reserve(2);
                if(!_lattice.getSite(neighbors[0]).isActive()){
                    bonds.push_back({BondType::Horizontal, site.row_, site.column_});
                }
                if(!_lattice.getSite(neighbors[1]).isActive()){
                    bonds.push_back({BondType::Vertical, prev_row, site.column_});
                }

            } else if (site.column_ == max_index) {
                // lower right corner
                neighbors.resize(2);
                neighbors[0] = {site.row_, prev_column};
                neighbors[1] = {prev_row, site.column_};

                bonds.reserve(2);
                if(!_lattice.getSite(neighbors[0]).isActive()){
                    bonds.push_back({BondType::Horizontal, site.row_, prev_column});
                }
                if(!_lattice.getSite(neighbors[1]).isActive()){
                    bonds.push_back({BondType::Vertical, prev_row, site.column_});
                }

            } else {
                // bottom edge excluding corners
                //  bottom edge
                neighbors.resize(3);
                neighbors[0] = {site.row_, next_column};
                neighbors[1] = {site.row_, prev_column};
                neighbors[2] = {prev_row, site.column_};

                bonds.reserve(3);
                if(!_lattice.getSite(neighbors[0]).isActive()) {
                    bonds.push_back({BondType::Horizontal, site.row_, site.column_});
                }
                if(!_lattice.getSite(neighbors[1]).isActive()){
                    bonds.push_back({BondType::Horizontal, site.row_, prev_column});
                }
                if(!_lattice.getSite(neighbors[2]).isActive()){
                    bonds.push_back({BondType::Vertical, prev_row, site.column_});
                }
            }
        }
            /* site.x_ > min_index && site.x_ < max_index &&  is not possible anymore*/
        else if (site.column_ == min_index) { // left edge not in the corners
            neighbors.resize(3);
            neighbors[0] = {site.row_, next_column};
            neighbors[1] = {next_row, site.column_};
            neighbors[2] = {prev_row, site.column_};

            bonds.reserve(3);
            if(!_lattice.getSite(neighbors[0]).isActive()) {
                bonds.push_back({BondType::Horizontal, site.row_, site.column_});
            }
            if(!_lattice.getSite(neighbors[1]).isActive()){
                bonds.push_back({BondType::Vertical,    site.row_, site.column_});
            }
            if(!_lattice.getSite(neighbors[2]).isActive()){
                bonds.push_back({BondType::Vertical, prev_row, site.column_});
            }
        }
        else if (site.column_ == max_index) {
            // right edge no corners

            neighbors.resize(3);
            neighbors[0] = {site.row_, prev_column};
            neighbors[1] = {next_row, site.column_};
            neighbors[2] = {prev_row, site.column_};

            bonds.reserve(3);
            if(!_lattice.getSite(neighbors[0]).isActive()){
                bonds.push_back({BondType::Horizontal, site.row_, prev_column});
            }
            if(!_lattice.getSite(neighbors[1]).isActive()){
                bonds.push_back({BondType::Vertical,    site.row_, site.column_});
            }
            if(!_lattice.getSite(neighbors[2]).isActive()){
                bonds.push_back({BondType::Vertical, prev_row, site.column_});
            }
        }
        else {
            // 1 level inside the lattice
            // not in any the boundary
            neighbors.resize(4);
            neighbors[0] = {site.row_, next_column};
            neighbors[1] = {site.row_, prev_column};
            neighbors[2] = {next_row, site.column_};
            neighbors[3] = {prev_row, site.column_};

            bonds.reserve(4);
            if(!_lattice.getSite(neighbors[0]).isActive()) {
                bonds.push_back({BondType::Horizontal, site.row_, site.column_});
            }
            if(!_lattice.getSite(neighbors[1]).isActive()){
                bonds.push_back({BondType::Horizontal, site.row_, prev_column});
            }
            if(!_lattice.getSite(neighbors[2]).isActive()){
                bonds.push_back({BondType::Vertical,    site.row_, site.column_});
            }
            if(!_lattice.getSite(neighbors[3]).isActive()){
                bonds.push_back({BondType::Vertical, prev_row, site.column_});
            }
        }

    }


}



/**
 * Take a bond index only if the corresponding site is active
 * takes longer? time than version 1?, i.e.,  connection()
 * @param site
 * @param site_neighbor
 * @param bond_neighbor
 */
void SitePercolation_ps_v10::connection_v2(Index site, vector<Index> &site_neighbor, vector<BondIndex> &bond_neighbor)
{

    value_type prev_column  = (site.column_ + length() - 1) % length();
    value_type prev_row     = (site.row_ + length() - 1) % length();
    value_type next_row     = (site.row_ + 1) % length();
    value_type next_column  = (site.column_ + 1) % length();

    if(!_periodicity){
        // without periodicity
        if (site.row_ == min_index) { // top edge including corners
            if(site.column_ == min_index){
                // upper left corner

                site_neighbor.resize(2);
                site_neighbor[0] = {site.row_, next_column};
                site_neighbor[1] = {next_row, site.column_};

                bond_neighbor.reserve(2);
                if(!_lattice.getSite(site_neighbor[0]).isActive()){
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
                }
                if(!_lattice.getSite(site_neighbor[1]).isActive()){
                    bond_neighbor.push_back({BondType::Vertical, site.row_, site.column_});
                }

                return;

            }
            else if(site.column_ == max_index){
                // upper right corner

                site_neighbor.resize(2);
                site_neighbor[0] = {site.row_, prev_column};
                site_neighbor[1] = {next_row, site.column_};

                bond_neighbor.reserve(2);
                if(!_lattice.getSite(site_neighbor[0]).isActive()){
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
                }
                if(!_lattice.getSite(site_neighbor[1]).isActive()){
                    bond_neighbor.push_back({BondType::Vertical, site.row_, site.column_});
                }

                return;
            }
            else{
                // top edge excluding corners
                site_neighbor.resize(3);
                site_neighbor[0] = {site.row_, next_column};
                site_neighbor[1] = {site.row_, prev_column};
                site_neighbor[2] = {next_row, site.column_};

                bond_neighbor.reserve(4);
                if(!_lattice.getSite(site_neighbor[0]).isActive()) {
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
                }
                if(!_lattice.getSite(site_neighbor[1]).isActive()){
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
                }
                if(!_lattice.getSite(site_neighbor[2]).isActive()){
                    bond_neighbor.push_back({BondType::Vertical,    site.row_, site.column_});
                }

                return;

            }
        }
        else if (site.row_ == max_index) { // bottom edge including corners
            if (site.column_ == min_index) {
                // lower left corner
                site_neighbor.resize(2);
                site_neighbor[0] = {site.row_, next_column};
                site_neighbor[1] = {prev_row, site.column_};

                bond_neighbor.reserve(2);
                if(!_lattice.getSite(site_neighbor[0]).isActive()){
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
                }
                if(!_lattice.getSite(site_neighbor[1]).isActive()){
                    bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
                }


                return;

            } else if (site.column_ == max_index) {
                // lower right corner
                site_neighbor.resize(2);
                site_neighbor[0] = {site.row_, prev_column};
                site_neighbor[1] = {prev_row, site.column_};

                bond_neighbor.reserve(2);
                if(!_lattice.getSite(site_neighbor[0]).isActive()){
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
                }
                if(!_lattice.getSite(site_neighbor[1]).isActive()){
                    bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
                }

                return;

            } else {
                // bottom edge excluding corners
                //  bottom edge
                site_neighbor.resize(3);
                site_neighbor[0] = {site.row_, next_column};
                site_neighbor[1] = {site.row_, prev_column};
                site_neighbor[2] = {prev_row, site.column_};

                bond_neighbor.reserve(3);
                if(!_lattice.getSite(site_neighbor[0]).isActive()) {
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
                }
                if(!_lattice.getSite(site_neighbor[1]).isActive()){
                    bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
                }
                if(!_lattice.getSite(site_neighbor[2]).isActive()){
                    bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
                }

                return;
            }
        }
            /* site.x_ > min_index && site.x_ < max_index &&  is not possible anymore*/
        else if (site.column_ == min_index) { // left edge not in the corners
            site_neighbor.resize(3);
            site_neighbor[0] = {site.row_, next_column};
            site_neighbor[1] = {next_row, site.column_};
            site_neighbor[2] = {prev_row, site.column_};

            bond_neighbor.reserve(3);
            if(!_lattice.getSite(site_neighbor[0]).isActive()) {
                bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
            }
            if(!_lattice.getSite(site_neighbor[1]).isActive()){
                bond_neighbor.push_back({BondType::Vertical,    site.row_, site.column_});
            }
            if(!_lattice.getSite(site_neighbor[2]).isActive()){
                bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
            }

            return;
        }
        else if (site.column_ == max_index) {
            // right edge no corners

            site_neighbor.resize(3);
            site_neighbor[0] = {site.row_, prev_column};
            site_neighbor[1] = {next_row, site.column_};
            site_neighbor[2] = {prev_row, site.column_};

            bond_neighbor.reserve(3);
            if(!_lattice.getSite(site_neighbor[0]).isActive()){
                bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
            }
            if(!_lattice.getSite(site_neighbor[1]).isActive()){
                bond_neighbor.push_back({BondType::Vertical,    site.row_, site.column_});
            }
            if(!_lattice.getSite(site_neighbor[2]).isActive()){
                bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
            }

            return;
        }

    }
    // 1 level inside the lattice
    // not in any the boundary
    site_neighbor.resize(4);
    site_neighbor[0] = {site.row_, next_column};
    site_neighbor[1] = {site.row_, prev_column};
    site_neighbor[2] = {next_row, site.column_};
    site_neighbor[3] = {prev_row, site.column_};

    bond_neighbor.reserve(4);
    if(!_lattice.getSite(site_neighbor[0]).isActive()) {
        bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
    }
    if(!_lattice.getSite(site_neighbor[1]).isActive()){
        bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
    }
    if(!_lattice.getSite(site_neighbor[2]).isActive()){
        bond_neighbor.push_back({BondType::Vertical,    site.row_, site.column_});
    }
    if(!_lattice.getSite(site_neighbor[3]).isActive()) {
        bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
    }

}


/**
 * Take a bond index only if the corresponding site is active
 * takes longer? time than version 1?, i.e.,  connection()
 * @param site
 * @param site_neighbor
 * @param bond_neighbor
 */
void SitePercolation_ps_v10::connection_periodic(Index site, vector<Index> &site_neighbor, vector<BondIndex> &bond_neighbor)
{
    site_neighbor.clear();
    bond_neighbor.clear();

    value_type prev_column  = (site.column_ + length() - 1) % length();
    value_type prev_row     = (site.row_ + length() - 1) % length();
    value_type next_row     = (site.row_ + 1) % length();
    value_type next_column  = (site.column_ + 1) % length();

    // 1 level inside the lattice
    // not in any the boundary
    site_neighbor.resize(4);
    site_neighbor[0] = {site.row_, next_column};
    site_neighbor[1] = {site.row_, prev_column};
    site_neighbor[2] = {next_row, site.column_};
    site_neighbor[3] = {prev_row, site.column_};

    bond_neighbor.reserve(4);

//    if(!_lattice.getSite(site_neighbor[0]).isActive()) {
//        bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
//    }
//    if(!_lattice.getSite(site_neighbor[1]).isActive()){
//        bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
//    }
//    if(!_lattice.getSite(site_neighbor[2]).isActive()){
//        bond_neighbor.push_back({BondType::Vertical,    site.row_, site.column_});
//    }
//    if(!_lattice.getSite(site_neighbor[3]).isActive()) {
//        bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
//    }
    bond_neighbor.resize(4);
    bond_neighbor[0] = {BondType::Horizontal, site.row_, site.column_};
    bond_neighbor[1] = {BondType::Horizontal, site.row_, prev_column};
    bond_neighbor[2] = {BondType::Vertical,    site.row_, site.column_};
    bond_neighbor[3] = {BondType::Vertical, prev_row, site.column_};

}



/**
 * This function must be called after a site is placed in the
 * newest cluster
 * and an id is given from the lattice according to the cluster they are on.
 * @param site
 */
void SitePercolation_ps_v10::save_index_if_in_boundary_v2(const Index& site){
    cout << "site with new id found " << site << " : line " << __LINE__ << endl;
    // checking row indices for Top-Bottom boundary
    if(site.row_ == min_index){ // top index
        if(!check_if_id_matches(site, _top_edge)) {
            _top_edge.push_back(site);
        }
    }
    else if(site.row_ == max_index){
        if(!check_if_id_matches(site, _bottom_edge)) {
            _bottom_edge.push_back(site);
        }
    }

    // checking column indices for Left-Right boundary
    if(site.column_ == min_index){
        if(!check_if_id_matches(site, _left_edge)) {
            _left_edge.push_back(site);
        }
    }
    else if(site.column_ == max_index){
        if(!check_if_id_matches(site, _right_edge)) {
            _right_edge.push_back(site);
        }
    }
}

/**
 *
 * @param site
 * @param edge
 * @return
 */
bool SitePercolation_ps_v10::check_if_id_matches(Index site, const vector<Index> &edge){
    for(auto s :edge){
        if(_lattice.getSite(site).get_groupID() == _lattice.getSite(s).get_groupID()){
            // no need to put the site here
//            cout << "Site " << site << " and Id " << _lattice.getSite(site).set_groupID()
//                 << " is already in the edge : line " << __LINE__ << endl;
            return true;
        }
    }
    return false;
}


/**
 * if at least one element is erased that means we have found a match and we have erased it
 * @param site
 * @param edge
 * @return
 */
bool SitePercolation_ps_v10::check_if_id_matches_and_erase(Index site, vector<Index> &edge){
    vector<Index>::iterator it = edge.begin();
    value_type sz = edge.size();
    for(; it < edge.end(); ++it){
        if(_lattice.getSite(site).get_groupID() == _lattice.getSite(*it).get_groupID()){
            // no need to put the site here
            cout << "Site " << site << " and Id " << _lattice.getSite(site).get_groupID()
                 << " is already in the edge : line " << __LINE__ << endl;
            edge.erase(it);
        }
    }
    // if at least one element is erased that means we have found a match and we have erased it
    return edge.size() < sz;
}


/***********************************************
 *  Placing sites
 *
 *****************************************/

/**
 * All site placing method in one place
 *
 * @return true if operation is successfull
 */
bool SitePercolation_ps_v10::occupy() {
    if(_index_sequence_position >= maxSites()){
        return false;
    }

    Index site = selectSite();
//    cout << "******************current site " << site  << "***********" << endl;
    placeSite(site);
//    _occuption_probability = occupationProbability(); // for super class
    return true;
}



value_type SitePercolation_ps_v10::placeSite(Index current_site) {
    return placeSite_weighted_v2(current_site);
}

value_type SitePercolation_ps_v10::placeSite(Index current_site, const std::vector<Index>& sites, const std::vector<BondIndex> &bonds) {
    return placeSite_weighted_v3(current_site, sites, bonds);
}

/***
 * Index of the selected site must be provided with the argument
 * Date : 2019.06.24
 * @param current_site
 * @return
 */
value_type SitePercolation_ps_v10::placeSite_weighted_v2(Index current_site) {
    // randomly choose a site
    if (_number_of_occupied_sites == maxSites()) {
        return ULONG_MAX;// unsigned long int maximum value
    }

    _last_placed_site = current_site;
//    cout << "placing site " << current_site << endl;
    _lattice.activate_site(current_site);
    ++_number_of_occupied_sites;

    // find the bonds for this site
    vector<BondIndex> bonds ;
    vector<Index> sites ;
    _lattice.get_neighbors(current_site, sites, bonds);
//    cout << bonds << endl;
    subtract_entropy_for_bond(bonds);  // tracking entropy change

    value_type base = manageClusters(sites, bonds);
    add_entropy_for_bond(base); // tracking entropy change
    // running tracker
    track_largestCluster(base);
//    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster
//    track_numberOfSitesInLargestCluster();

    return base;
}

/***
 * Index of the selected site must be provided with the argument
 * Date : 2019.06.24
 * @param current_site
 * @return
 */
value_type SitePercolation_ps_v10::placeSite_weighted_v3(
        Index current_site, const std::vector<Index>& sites, const vector<BondIndex> &bonds)
{
    // randomly choose a site
    if (_number_of_occupied_sites == maxSites()) {
        return ULONG_MAX;// unsigned long int maximum value
    }

    _last_placed_site = current_site;
//    cout << "placing site " << current_site << endl;
    _lattice.activate_site(current_site);
    ++_number_of_occupied_sites;

//    cout << bonds << endl;
    subtract_entropy_for_bond(bonds);  // tracking entropy change

    value_type base = manageClusters(sites, bonds);
    add_entropy_for_bond(base); // tracking entropy change
    // running tracker
    track_largestCluster(base);
//    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster
//    track_numberOfSitesInLargestCluster();

    return base;
}

/**
 * Works fine
 *
 * @param sites
 * @param bonds
 * @return
 */
value_type SitePercolation_ps_v10::merge_cluster_v2(
        const std::vector<Index>& sites, const vector<BondIndex> &bonds) {// mergeing procedure
    // we got four bonds
    int base{}, id{};
    size_t base_sites{}, size{}, base_bonds{};
    unordered_set<int> ids; // selected cluster to merge including the base cluster id
    for(const BondIndex &b: bonds){
        // choosing the largest cluster as base cluster
        id = _lattice.getGroupID(b);
        size = _clusters[id].numberOfSites();
        if(size >= base_sites){
            // greater equal is important here.
            // since we are measuring with site and a cluster
            // can have zero sites in bond percolation
            base_sites = size;
            base = id;
        }
        ids.insert(id);
    }

    // find which of the neighbors are of id_base as the base cluster
    IndexRelative r;
    for(auto n: sites){
        if(_lattice.getSite(n).get_groupID() == base){
            // find relative index with respect to this site
            r = getRelativeIndex(n, _last_placed_site);
            break; // since first time r is set running loop is doing no good
        }
    }
    base_bonds = _clusters[base].numberOfBonds();
    _clusters[base].addSiteIndex(_last_placed_site);
    _lattice.setRelativeIndex(_last_placed_site, r);
    _lattice.setGroupID(_last_placed_site, base);
//    cout << "base " << base << " ids {";
    for(int i: ids){
        if(i == base) {
            continue;
        }
//        cout << i << ",";
        _clusters[base].insert_v2(_clusters[i]);
        relabel_sites_v5(_last_placed_site, _clusters[i]);
        relabel_bonds(_clusters[i], base);

//        relabel_cluster(_last_placed_site, _clusters[i], 0, 0);
        _clusters[i].clear(); // clear the cluster
        --number_of_clusters;
    }

//    relabel_cluster(_last_placed_site, _clusters[base], base_bonds, base_sites);
//    cout << "}" << endl;
#ifdef UNIT_TEST
    if(base < 0){
        cerr << "base cannot be negative: SitePercolation_ps_v10::merge_cluster_v2" << endl;
        exit(0);
    }
#endif
    return base;
}

/**
 * todo On test
 *
 * @param sites
 * @param bonds
 * @return
 */
value_type SitePercolation_ps_v10::merge_cluster_v3(
        const std::vector<Index>& sites, const vector<BondIndex> &bonds) {// mergeing procedure
    // we got four bonds
    int base{}, id{};
    size_t base_sites{}, size{}, base_bonds{};
    unordered_set<int> ids; // selected cluster to merge including the base cluster id
//    cout << "found ids {";
    for(const BondIndex &b: bonds){
        // choosing the largest cluster as base cluster
        id = _lattice.getGroupID(b);
//        cout << id << ",";
        if(_clusters[id].emptyBond()) continue; // empty cluster cannot be a part of what's coming next
        size = _clusters[id].numberOfSites();
        if(size >= base_sites){
            // greater equal is important here.
            // since we are measuring with site and a cluster
            // can have zero sites in bond percolation
            base_sites = size;
            base = id;
        }
        ids.insert(id);
    }
//    cout << "}" << endl;
//    cout << "non empty {";
//    copy(ids.begin(), ids.end(), ostream_iterator<int>(cout, ","));
//    cout << "}" << endl;

    // find which of the neighbors are of id base
    IndexRelative r;
    for(auto n: sites){
        if(_lattice.getGroupID(n) == base){
            // find relative index with respect to this site
            r = getRelativeIndex(n, _last_placed_site);
            break; // since first time r is set running loop is doing no good
        }
    }
    base_bonds = _clusters[base].numberOfBonds();
    _clusters[base].addSiteIndex(_last_placed_site);
    _lattice.setRelativeIndex(_last_placed_site, r);
    _lattice.setGroupID(_last_placed_site, base);
    for(int i: ids){
        if(i == base) {
            continue;
        }
        _clusters[base].insert_v2(_clusters[i]);
        relabel_cluster(_last_placed_site, _clusters[i], 0, 0);
        _clusters[i].clear(); // clear the cluster
        --number_of_clusters;
    }
#ifdef UNIT_TEST
    if(base < 0){
        cerr << "base cannot be negative: SitePercolation_ps_v10::merge_cluster_v3" << endl;
        exit(0);
    }
#endif
    return base;
}




/**
 *
 * @return
 */
Index SitePercolation_ps_v10::selectSite(){
//    Index current_site = randomized_index_sequence[_index_sequence_position]; // old
    value_type index = randomized_index[_index_sequence_position];
    Index current_site = site_indices[index]; // new process
    ++_index_sequence_position;

    return current_site;
}


/***************************************************
 * View methods
 ****************************************/


/**
 *
 */
void SitePercolation_ps_v10::spanningIndices() const {
    cout << "Spanning Index : id" << endl;
    for(auto i: _spanning_sites){
        cout << "Index " << i << " : id " << _lattice.getSite(i).get_groupID() << endl;
    }
}

void SitePercolation_ps_v10::wrappingIndices() const {
    cout << "Wrapping Index : id : relative index" << endl;
    for(auto i: _wrapping_sites){
        cout << "Index " << i << " : id "
             << _lattice.getSite(i).get_groupID()
             << " relative index : " << _lattice.getSite(i).relativeIndex() << endl;
    }
}


/****************************************
 * Spanning Detection
 ****************************************/


/**
 * todo problem
 * successful : but result is not entirely correct
 * length       time
 * 200          5.265000 sec
 * 500          2 min 29.656000
 * @param site : Check spanning for this argument
 * @return
 */
bool SitePercolation_ps_v10::detectSpanning_v5() {
//    cout << "Entry -> detectSpanning_v4() : line " << __LINE__ << endl;
    if(_periodicity) {
        cout << "Cannot detect spanning if _periodicity if ON: line " << __LINE__ << endl;
        return false;
    }
    Index site = _last_placed_site;

    // first check if the site with a cluster id is already a spanning site
    for(const Index& ss: _spanning_sites){
        if(_lattice.getSite(ss).get_groupID() == _lattice.getSite(site).get_groupID()){
            cout << "Already a spanning site : line " << __LINE__ << endl;
            return true;
        }
    }

    // only check for the newest site placed
    if(site.row_ == min_index){ // top index
        _top_edge.push_back(site);
    }
    else if(site.row_ == max_index){
        _bottom_edge.push_back(site);
    }

    // checking column indices for Left-Right boundary
    if(site.column_ == min_index){ // left edge
        _left_edge.push_back(site);
    }
    else if(site.column_ == max_index){
        _right_edge.push_back(site);
    }

    if(_index_sequence_position < length()){
//        cout << "Not enough site to span : line " << __LINE__ << endl;
        return false;
    }

    // now do the matching with top and bottom for vertical spanning
    // meaning new site is added to _spanning_sites so remove them from top and bottom edges
    vector<Index>::iterator it_top = _top_edge.begin();
    vector<Index>::iterator it_bot = _bottom_edge.begin();
    bool found_spanning_site = false;
    for(; it_top < _top_edge.end(); ++it_top){
        for(; it_bot < _bottom_edge.end(); ++it_bot){
            if(_lattice.getSite(*it_top).get_groupID() == _lattice.getSite(*it_bot).get_groupID()){
                _spanning_sites.push_back(*it_top);
//                _spanning_occured = true;
//                spanning_id = _lattice.getSite(*it_bot).set_groupID();
                found_spanning_site = true;
                _bottom_edge.erase(it_bot);
            }
        }
        if(found_spanning_site){
            found_spanning_site = false;
            _top_edge.erase(it_top);
        }
    }


    found_spanning_site = false;
    // now do the matching with left and right for horizontal spanning
    // meaning new site is added to _spanning_sites so remove them from top and bottom edges
    vector<Index>::iterator it_lft = _left_edge.begin();
    vector<Index>::iterator it_rht = _right_edge.begin();
    for(; it_lft < _left_edge.end(); ++it_lft){
        for(; it_rht < _right_edge.end(); ++it_rht){
            if(_lattice.getSite(*it_lft).get_groupID() == _lattice.getSite(*it_rht).get_groupID()){
                _spanning_sites.push_back(*it_lft);
//                _spanning_occured = true;
//                spanning_id = _lattice.getSite(*it_bot).set_groupID();
                found_spanning_site = true;
                _right_edge.erase(it_rht);
            }
        }
        if(found_spanning_site){
            found_spanning_site = false;
            _left_edge.erase(it_lft);
        }
    }


    return !_spanning_sites.empty();
}


/**
 * success : gives correct result
 * length       time
 * 200          7.859000 sec
 * 500          2 min 18.874000 sec
 * only check for the cluster id of the recently placed site
 * @param site : Check spanning for this argument
 * @return
 */
bool SitePercolation_ps_v10::detectSpanning_v6() {
//    cout << "Entry -> detectSpanning_v4() : line " << __LINE__ << endl;
    if(_periodicity) {
        cout << "Cannot detect spanning if _periodicity if ON: line " << __LINE__ << endl;
        return false;
    }
    Index site = _last_placed_site;
    // first check if the site with a cluster id is already a spanning site
    for(const Index& ss: _spanning_sites){
        if(_lattice.getSite(ss).get_groupID() == _lattice.getSite(site).get_groupID()){
//            cout << "Already a spanning site : line " << __LINE__ << endl;
            return true;
        }
    }

    // only check for the newest site placed
    if(site.row_ == min_index){ // top index
        if(!check_if_id_matches(site, _top_edge)) {
            _top_edge.push_back(site);
        }
    }
    else if(site.row_ == max_index){
        if(!check_if_id_matches(site, _bottom_edge)){
            _bottom_edge.push_back(site);
        }
    }

    // checking column indices for Left-Right boundary
    if(site.column_ == min_index){ // left edge
        if(!check_if_id_matches(site, _left_edge)) {
            _left_edge.push_back(site);
        }
    }
    else if(site.column_ == max_index){
        if(!check_if_id_matches(site, _right_edge)) {
            _right_edge.push_back(site);
        }
    }

    if(_number_of_occupied_sites < length()){
//        cout << "Not enough site to span : line " << __LINE__ << endl;
        return false;
    }


    vector<Index>::iterator it_top = _top_edge.begin();
    vector<Index>::iterator it_bot = _bottom_edge.begin();
    bool found_spanning_site = false;
    int id = _lattice.getSite(site).get_groupID();

    if(_top_edge.size() < _bottom_edge.size()){
        // if matched found on the smaller edge look for match in the larger edge
        for(; it_top < _top_edge.end(); ++it_top){
            if(id == _lattice.getSite(*it_top).get_groupID()){
                for(; it_bot < _bottom_edge.end(); ++it_bot){
                    if(id == _lattice.getSite(*it_bot).get_groupID()){
                        // match found !
                        if(!check_if_id_matches(*it_top ,_spanning_sites)) {
//                            _spanning_occured = true;
                            _spanning_sites.push_back(*it_top);
                        }
                        found_spanning_site = true;
                        _bottom_edge.erase(it_bot);
                    }
                }

                if(found_spanning_site){
                    found_spanning_site = false;
                    _top_edge.erase(it_top);
                }

            }
        }
    }else{
        for (; it_bot < _bottom_edge.end(); ++it_bot) {
            if (id == _lattice.getSite(*it_bot).get_groupID()) {
                for (; it_top < _top_edge.end(); ++it_top) {
                    if (id == _lattice.getSite(*it_top).get_groupID()) {
                        // match found !
                        if (!check_if_id_matches(*it_top, _spanning_sites)) {
//                            _spanning_occured = true;
                            _spanning_sites.push_back(*it_top);
                        }
                        found_spanning_site = true;
                        _top_edge.erase(it_top);
                    }
                }
                if(found_spanning_site){
                    found_spanning_site = false;
                    _bottom_edge.erase(it_top);
                }
            }
        }

    }

    found_spanning_site = false;
    vector<Index>::iterator it_lft = _left_edge.begin();
    vector<Index>::iterator it_rht = _right_edge.begin();

    if(_left_edge.size() < _right_edge.size()){
        for(; it_lft < _left_edge.end(); ++it_lft) {
            if (id == _lattice.getSite(*it_lft).get_groupID()) {
                for (; it_rht < _right_edge.end(); ++it_rht) {
                    if (id == _lattice.getSite(*it_rht).get_groupID()) {
                        if (!check_if_id_matches(*it_lft, _spanning_sites)) {
                            _spanning_sites.push_back(*it_lft);
//                            _spanning_occured = true;
                        }
                        found_spanning_site = true;
                        _right_edge.erase(it_rht);
                    }
                }
                if (found_spanning_site) {
                    found_spanning_site = false;
                    _left_edge.erase(it_lft);
                }
            }
        }
    }else{
        for (; it_rht < _right_edge.end(); ++it_rht) {
            if (id == _lattice.getSite(*it_rht).get_groupID()) {
                for(; it_lft < _left_edge.end(); ++it_lft) {
                    if (id == _lattice.getSite(*it_lft).get_groupID()) {
                        if (!check_if_id_matches(*it_lft, _spanning_sites)) {
                            _spanning_sites.push_back(*it_lft);
//                            _spanning_occured = true;
                        }
                        found_spanning_site = true;
                        _left_edge.erase(it_lft);
                    }
                }
                if (found_spanning_site) {
                    found_spanning_site = false;
                    _right_edge.erase(it_rht);
                }
            }
        }
    }


    // now do the matching with left and right for horizontal spanning
    // meaning new site is added to _spanning_sites so remove them from top and bottom edges



    // filter spanning ids


    return !_spanning_sites.empty();

}


/**
 * todo scan the edges and save it in the edge indices
 */
void SitePercolation_ps_v10::scanEdges() {

}




/***********************************
 * Wrapping Detection
 **********************************/
/**
 * TODO : does not work probably because
 * Wrapping is detected here using the last placed site
 * @return bool. True if wrapping occured.
 */
bool SitePercolation_ps_v10::detectWrapping() {
    return detect_wrapping_v1();
}

bool SitePercolation_ps_v10::detect_wrapping_v1()  {
    Index site = lastPlacedSite();
    // only possible if the cluster containing 'site' has site_index_sequence >= length of the lattice
    if(_number_of_occupied_sites < length()){
        return false;
    }

    // check if it is already a wrapping site
    int id = _lattice.getSite(site).get_groupID();
    int tmp_id{};
    for (auto i: _wrapping_sites){
        tmp_id = _lattice.getSite(i).get_groupID();
        if(id == tmp_id ){
//            cout << "Already a wrappig cluster : line " << __LINE__ << endl;
            return true;
        }
    }

    // get four neighbors of site always. since wrapping is valid if periodicity is implied
    vector<Index> sites = _lattice.get_neighbor_sites(site);


    if(sites.size() < 2){ // at least two neighbor of  site is required
        cerr << "sites.size() < 2  : line "  << __LINE__ << endl;
        return false;
    }else{
        IndexRelative irel = _lattice.getSite(site).relativeIndex();
//        cout << "pivot's " << site << " relative " << irel << endl;
        IndexRelative b;
        for (auto a:sites){
            if(_lattice.getSite(a).get_groupID() != _lattice.getSite(site).get_groupID()){
                // different cluster
                continue;
            }
//            cout << "belongs to the same cluster : line " << __LINE__ << endl;

            b = _lattice.getSite(a).relativeIndex();
//            cout << "neibhbor " << a << " relative " << b << endl;
            if(abs(irel.x_ - b.x_) > 1 || abs(irel.y_ - b.y_) > 1){
//                cout << "Wrapping : line " << __LINE__ << endl;
                _wrapping_sites.push_back(site);
                return true;
            }
        }
    }

//    cout << "wrapping site_index_sequence " << _wrapping_indices << endl;
    // if %_wrapping_indices is not empty but wrapping is not detected for the current site (%site)
    // that means there is wrapping but not for the %site
    return !_wrapping_sites.empty();
}



/********************************************************************
 * Relabeling
 *
 *********************************************************************/


/**
 *
 * @param clstr
 * @param id
 */
void SitePercolation_ps_v10::relabel_sites(const Cluster_v3& clstr, int id) {
    const vector<Index> sites = clstr.getSiteIndices();
    for(auto a: sites){
        _lattice.getSite(a).set_groupID(id);
    }
}


/**
 * Relabels site and also reassign relative index to the relabeled sites
  *
  * @param site_a  : root index of the base cluster
  * @param clstr_b : 2nd cluster, which to be merged withe the root
  */
void SitePercolation_ps_v10::relabel_sites_v4(Index site_a, const Cluster_v3& clstr_b) {
    const vector<Index> sites = clstr_b.getSiteIndices();
    int id_a = _lattice.getSite(site_a).get_groupID();
    int id_b = clstr_b.get_ID();
    Index b = clstr_b.getRootSite();

    // get four site_b of site_a
    vector<Index> sites_neighbor_a = _lattice.get_neighbor_sites(site_a);
    Index site_b;
    IndexRelative relative_index_b_after;
    bool flag{false};
    // find which site_b has id_a of clstr_b
    for(auto n: sites_neighbor_a){
        if(id_b == _lattice.getSite(n).get_groupID()){
            // checking id_a equality is enough. since id_a is the id_a of the active site already.
            relative_index_b_after = getRelativeIndex(site_a, n);
            site_b = n;
//            cout << "neighbor  of" << site_a << " is " << site_b << endl;
            flag = true;
            break; // todo ?

        }
    }

    if(!flag){
        cout << "No neibhgor found! : line " << __LINE__ << endl;
    }


    IndexRelative relative_site_a = _lattice.getSite(site_a).relativeIndex();

    // with this delta_a and delta_y find the relative index of site_b while relative index of site_a is known
    IndexRelative relative_site_b_before = _lattice.getSite(site_b).relativeIndex();

    int delta_x_ab = relative_index_b_after.x_ - relative_site_b_before.x_;
    int delta_y_ab = relative_index_b_after.y_ - relative_site_b_before.y_;

//    cout << relative_index_b_after << " - " << relative_site_b_before << " = delta_x, delta_y = " << delta_x_ab << ", " << delta_y_ab << endl;
    int x, y;

    for (auto a : sites) {
        _lattice.getSite(a).set_groupID(id_a);
        relative_site_a = _lattice.getSite(a).relativeIndex();

//        cout << "Relative " << relative_site_a << endl;
        x = relative_site_a.x_ + delta_x_ab;
        y = relative_site_a.y_ + delta_y_ab;
//        cout << "(x,y) = (" << x << "," << y << ") : line " << __LINE__ << endl;
        _lattice.getSite(a).relativeIndex(x, y);
    }
}


/**
 * Relabels site and also reassign relative index to the relabeled sites
  *
  * @param site_a  : last added site index of the base cluster
  * @param clstr_b : 2nd cluster, which to be merged withe the root
  */
void SitePercolation_ps_v10::relabel_sites_v5(Index site_a, const Cluster_v3& clstr_b) {
    if (clstr_b.emptySite()) {
//        cout << "empty cluster cannot have root site : relabel_sites_v5 : line " << __LINE__ << endl;
        return;
    }
    const vector<Index> & sites = clstr_b.getSiteIndices();
    int id_a = _lattice.getGroupID(site_a);
    int id_b = clstr_b.get_ID();
    Index b = clstr_b.getRootSite();
    // get four site_b of site_a
    vector<Index> sites_neighbor_a = _lattice.get_neighbor_sites(site_a);
    Index site_b;
    IndexRelative relative_index_b_after;
    bool flag{false};
    // find which site_b has id_a of clstr_b
    for(auto n: sites_neighbor_a){
        if(id_b == _lattice.getSite(n).get_groupID()){
            // checking id_a equality is enough. since id_a is the id_a of the active site already.
            relative_index_b_after = getRelativeIndex(site_a, n);
            site_b = n;
//            cout << "neighbor  of" << site_a << " is " << site_b << endl;
            flag = true;
            break; // todo ?
        }
    }
    if(!flag){
//        cout << "No neibhgor found! : line " << __LINE__ << endl;
    }
    IndexRelative relative_site_a = _lattice.getSite(site_a).relativeIndex();

    // with this delta_a and delta_y find the relative index of site_b while relative index of site_a is known
    IndexRelative relative_site_b_before = _lattice.getSite(site_b).relativeIndex();

    int delta_x_ab = relative_index_b_after.x_ - relative_site_b_before.x_;
    int delta_y_ab = relative_index_b_after.y_ - relative_site_b_before.y_;
//    cout << relative_index_b_after << " - " << relative_site_b_before << " = delta_x, delta_y = " << delta_x_ab << ", " << delta_y_ab << endl;
    relabel_sites(sites, id_a, delta_x_ab, delta_y_ab);
}



/**
 * Relabels site and also reassign relative index to the relabeled sites
 *
 * @param site_a  : root index of the base cluster or last added site index of the base cluster
 * @param clstr_b : 2nd cluster, which to be merged withe the root
 * @param id : id to be used for relabeling sites
 */
void SitePercolation_ps_v10::relabel_sites_v6(Index site_a, const Cluster_v3& clstr_b, int id) {
    const vector<Index> sites = clstr_b.getSiteIndices();
    int id_a = _lattice.getSite(site_a).get_groupID();
    int id_b = clstr_b.get_ID();
    Index b = clstr_b.getRootSite();

    // get four site_b of site_a
    vector<Index> sites_neighbor_a = _lattice.get_neighbor_sites(site_a);
    Index site_b;
    IndexRelative relative_index_b_after;
    bool flag{false};
    // find which site_b has id_a of clstr_b
    for(auto n: sites_neighbor_a){
        if(id_b == _lattice.getSite(n).get_groupID()){
            // checking id_a equality is enough. since id_a is the id_a of the active site already.
            relative_index_b_after = getRelativeIndex(site_a, n);
            site_b = n;
//            cout << "neighbor  of" << site_a << " is " << site_b << endl;
            flag = true;
            break; // todo ?

        }
    }

    if(!flag){
        cout << "No neibhgor found! : line " << __LINE__ << endl;
    }


    IndexRelative relative_site_a = _lattice.getSite(site_a).relativeIndex();

    // with this delta_a and delta_y find the relative index of site_b while relative index of site_a is known
    IndexRelative relative_site_b_before = _lattice.getSite(site_b).relativeIndex();

    int delta_x_ab = relative_index_b_after.x_ - relative_site_b_before.x_;
    int delta_y_ab = relative_index_b_after.y_ - relative_site_b_before.y_;

//    cout << relative_index_b_after << " - " << relative_site_b_before << " = delta_x, delta_y = " << delta_x_ab << ", " << delta_y_ab << endl;

    relabel_sites(sites, id, delta_x_ab, delta_y_ab);

}



void SitePercolation_ps_v10::relabel_sites(const vector<Index> &sites, int id_a, int delta_x_ab, int delta_y_ab)  {
    int x, y;
    Index a;
    IndexRelative relative_site__a;
    for (value_type i = 0; i < sites.size(); ++i) {
        a = sites[i];
//        _lattice.getSite(a).set_groupID(id_a);
        _lattice.setGroupID(a, id_a);
        relative_site__a = _lattice.getSite(a).relativeIndex();
        x = relative_site__a.x_ + delta_x_ab;
        y = relative_site__a.y_ + delta_y_ab;
        _lattice.getSite(a).relativeIndex(x, y);
    }
}

/**
 *
 * @param clstr
 * @param id
 */
void SitePercolation_ps_v10::relabel_bonds(const Cluster_v3& clstr, int id) {
    vector<BondIndex> bonds = clstr.getBondIndices();
    cout << "relabel_bonds : base id " << id << " {";
    for(auto a: bonds){
        _lattice.setGroupID(a, id);
        cout << a << ",";
    }
    cout << "}" << endl;
}


/**********************************************
 * Information about current state of Class
 **********************************************/
/**
 * number of bonds in spanning cluster / total number of bonds
 * here, total number of bonds = (2*l*l - 2*l)
 *  since periodic boundary condition is turned off
 * @return
 */
double SitePercolation_ps_v10::spanningProbability() const {
    // TODO
    cout << "Not implemented : line " << __LINE__ << endl;
    return 0;
}


/**
 * Entropy calculation is performed here. The fastest method possible.
 * Cluster size is measured by bond.
 * @return current entropy of the lattice
 */
long double SitePercolation_ps_v10::entropy() {
//    return entropy_v1();
    return entropy_v2();
}



/**
 * number of bonds in the largest cluster / total number of bonds
 *
 * here cluster size is measured by number of bonds in the cluster
 * @return
 */
double SitePercolation_ps_v10::orderParameter() const{
    if(_number_of_bonds_in_the_largest_cluster != 0){
        // already calculated somewhere
        return double(_number_of_bonds_in_the_largest_cluster) / maxBonds();
    }
    value_type  len{}, nob{};
    for(const Cluster_v3& c: _clusters){
        nob = c.numberOfBonds();
        if (len < nob){
            len = nob;
        }
    }
    return double(len) / maxBonds();
}

/**
 * number of bonds in the largest cluster / total number of bonds
 *
 * here cluster size is measured by number of bonds in the cluster
 * @return
 */
double SitePercolation_ps_v10::orderParameter_v2() const{
//    double bond_num = _clusters[_index_largest_cluster].numberOfBonds() / double(_total_bonds);
    double tmp = _number_of_bonds_in_the_largest_cluster / double(maxBonds());
    return tmp;
}




/**
 * Only applicable if the number of bonds in the largest cluster is calculated when occupying the lattice.
 * Significantly efficient than the previous version numberOfBondsInTheLargestCluster()
 * @return
 */
value_type SitePercolation_ps_v10::numberOfBondsInTheLargestCluster() {
//    return _clusters[_index_largest_cluster].numberOfBonds();
    return _number_of_bonds_in_the_largest_cluster;
}



/**
 *
 * @return
 */
value_type SitePercolation_ps_v10::numberOfSitesInTheLargestCluster() {
    value_type  len{}, nob{};
    for(auto c: _clusters){
        nob = c.numberOfSites();
        if (len < nob){
            len = nob;
        }
    }
    _number_of_sites_in_the_largest_cluster = len;
    return len;
}


/**********************************
 * Spanning methods
 **********************************/


/**
 *
 * @return
 */
value_type SitePercolation_ps_v10::numberOfSitesInTheSpanningClusters() {

    if(! _spanning_sites.empty()){
        int id = _lattice.getSite(_spanning_sites.front()).get_groupID();
        if(id >= 0) {
            return _clusters[id].numberOfSites();
        }
    }
    return 0;
}


/**
 *
 * @return
 */
value_type SitePercolation_ps_v10::numberOfBondsInTheSpanningClusters() {
    if(!_spanning_sites.empty()){
//        cout << "number of spanning site_index_sequence " << _spanning_sites.size() << " : line " << __LINE__ << endl;
        int id = _lattice.getSite(_spanning_sites.front()).get_groupID();
        if(id >= 0) {
            return _clusters[id].numberOfBonds();
        }
    }
    return 0;
}

/**
 *
 * @return
 */
value_type SitePercolation_ps_v10::numberOfSitesInTheWrappingClusters(){
    value_type nos{};
    int id{};
    for(auto i: _wrapping_sites){
        id = _lattice.getSite(i).get_groupID();
        if(id >= 0) {
            nos += _clusters[id].numberOfSites();
        }
    }
    return nos;
}

/**
 *
 * @return
 */
value_type SitePercolation_ps_v10::numberOfBondsInTheWrappingClusters(){
    value_type nob{};
    int id{};
    for(auto i: _wrapping_sites){
        id = _lattice.getSite(i).get_groupID();
        if(id >= 0) {
            nob += _clusters[id].numberOfBonds();
        }
    }
    return nob;
}


/**
 * _first_spanning_cluster_id must be positiove
 * @return
 */
int SitePercolation_ps_v10::birthTimeOfSpanningCluster() const {
    if(!_spanning_sites.empty()) {
        cout << "number of spanning site_index_sequence " << _spanning_sites.size() << " : line " << __LINE__ << endl;
//        Index site = _spanning_sites.front();
//        int id = _lattice.getSite(site).set_groupID();
//        value_type clster_index = _cluster_index_from_id[id];
//        return _clusters[clster_index].birthTime();
        cout << "TODO : line " << __LINE__ << endl;
        return -1;
    }

    cout << "No spanning cluster : " << __LINE__ << endl;
    return -1;
}

/**
 * _first_spanning_cluster_id must be positiove
 * @return
 */
int SitePercolation_ps_v10::birthTimeOfACluster(int id) const {
    if (id < 0){
        cerr << "_first_spanning_cluster_id < 0 : line " << __LINE__ << endl;
        return -1;
    }
    for (const Cluster_v3& cls: _clusters){
        if(cls.get_ID() == id){
            return cls.birthTime();
        }
    }
    return -1;
}



/**
 * Box counting method for fractal dimension calculation
 * @param delta
 * @return
 */
value_type SitePercolation_ps_v10::box_counting(value_type delta) {
    value_type counter{};
    if(length() % delta == 0){
        for(value_type r{} ; r < length() ; r += delta){
            for(value_type c{}; c < length() ; c += delta){
                if(anyActiveSite(r, c, delta)){
                    ++ counter;
                }
            }
        }
    }
    else{
        cout << "Delta size is not corrent" << endl;
    }

    return counter;
}


/**
 * Box counting method for fractal dimension calculation
 * @param delta
 * @return
 */
array<value_type, 2> SitePercolation_ps_v10::box_counting_v2(value_type delta) {
    int spanning_cluster_id = _lattice.getSite(_spanning_sites.front()).get_groupID();
    value_type counter{}, spanning_counter{};
    bool foud_site{false}, found_spanning_site{false};
    int current_id{-1};
    if(length() % delta == 0){
        for(value_type r{} ; r < length() ; r += delta){
            for(value_type c{}; c < length() ; c += delta){
                for(value_type i{} ; i < delta ; ++i) {
                    for (value_type j{}; j < delta; ++j) {
                        current_id = _lattice.getSite({i + r, j + c}).get_groupID();
                        if(current_id >= 0){
                            foud_site = true;
                            if(current_id == spanning_cluster_id){
                                found_spanning_site = true;
                            }
                        }
                        if(foud_site && found_spanning_site){
                            break; // both are true no need to investigate further
                        }
                    }
                }

                if(foud_site){
                    ++ counter;
                    foud_site = false;
                }
                if(found_spanning_site){
                    ++spanning_counter;
                    found_spanning_site = false;
                }
            }
        }
    }
    else{
        cout << "Delta size is not corrent" << endl;
    }

    return array<value_type , 2>{counter, spanning_counter};
}




/**
 * Box counting method for fractal dimension calculation
 * @param delta
 * @return
 */
value_type SitePercolation_ps_v10::box_counting_spanning(value_type delta) {
    value_type counter{};
    if(length() % delta == 0){
        for(value_type r{} ; r < length() ; r += delta){
            for(value_type c{}; c < length() ; c += delta){
                if(anyActiveSpanningSite(r, c, delta)){
                    ++counter;
                }
            }
        }
    }
    else{
        cout << "Delta size is not corrent" << endl;
    }

    return counter;
}


/**
 * Return true if there is at least one active site in the region r, r + delta and c, c + delta
 * @param r
 * @param c
 * @param delta
 * @return
 */
bool SitePercolation_ps_v10::anyActiveSite(value_type row, value_type col, value_type delta) {
    for(value_type r{} ; r < delta ; ++r) {
        for (value_type c{}; c < delta; ++c) {
            if(_lattice.getSite({r+row, c + col}).isActive()){
                return true;
            }
        }
    }
    return false;
}


/**
 * Return true if there is at least one active site in the region r, r + delta and c, c + delta
 * Also this function required that, the site is in the spanning cluster
 * @param row
 * @param col
 * @param delta
 * @return
 */
bool SitePercolation_ps_v10::anyActiveSpanningSite(value_type row, value_type col, value_type delta) {
    int spanning_cluster_id = _lattice.getSite(_spanning_sites.front()).get_groupID();
    for(value_type r{} ; r < delta ; ++r) {
        for (value_type c{}; c < delta; ++c) {
            if(_lattice.getSite({r + row, c + col}).get_groupID() == spanning_cluster_id){
                return true;
            }
        }
    }
    return false;
}



/********************************
 *
 ********************************/

/**
 * Counts the number of active sites in the lattice
 */
value_type SitePercolation_ps_v10::count_number_of_active_site() {
    value_type counter{};
    for (value_type i{}; i != length(); ++i) {
        for (value_type j{}; j != length(); ++j) {
            if (_lattice.getSite({i, j}).isActive())
                ++counter;
        }
    }
    cout << "Number of active site : " << counter << endl;
    return counter;
}

std::string SitePercolation_ps_v10::getSignature() {
    string s = "sq_lattice_site_percolation";
    if(_periodicity)
        s += "_periodic";
    else
        s += "_non_periodic";
    return s;
}

/**
 *
 * @param filename
 * @param only_spanning
 */
void SitePercolation_ps_v10::writeVisualLatticeData(const string &filename, bool only_spanning) {
    std::ofstream fout(filename);
    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length()
                << R"(,"signature":")" << getSignature() << "\""
                << R"(,"x":")" << lastPlacedSite().column_ << "\""
                << R"(,"y":")" << lastPlacedSite().row_ << "\""
                << "}" ;

    fout << "#" << header_info.str() << endl;
//    fout << "#" << "spannig cluster"  << endl;
    fout << "#<x>,<y>,<color>" << endl;
    fout << "# color=0 -means-> spanning site" << endl;
    int id{-1}, spanning_id{-1};
    if(!_spanning_sites.empty()){
        spanning_id = _lattice.getSite(_spanning_sites.front()).get_groupID();
    }
    else if(!_wrapping_sites.empty()){
        spanning_id = _lattice.getSite(_wrapping_sites.front()).get_groupID();
    }

    if(only_spanning){
        if(spanning_id < 0){
            cerr << "id < 0 : line " << __LINE__ << endl;
        }
        vector<Index> sites = _clusters[id].getSiteIndices();
        for(auto s: sites){
            fout << s.column_ << ',' << s.row_ << ',' << 0 << endl;
        }
    }
    else {
        for (value_type y{}; y != length(); ++y) {
            for (value_type x{}; x != length(); ++x) {
                id = _lattice.getGroupID({y, x});
                if(id < 0) continue; // unoccupied
                if(id == spanning_id) {
                    fout << x << ',' << y << ',' << 0 << endl;
                }else{
                    fout << x << ',' << y << ',' << id + 1 << endl;
                }
            }
        }
    }
    fout.close();
}


/**
 * Not normalized
 * @return an array : index of the array is the cluster size.
 *                    value of the array is the number of times that cluster appears
 */
const std::vector<double> SitePercolation_ps_v10::clusterSizeDistribution() const{
    vector<double> cluster_counts;
    size_t n, sz;
//    cout << "number of max clusters " << _clusters.size() << endl;
    size_t mx_n{}, total_bonds{};
    for(size_t i{}; i < _clusters.size(); ++i){
        if(_clusters[i].empty())   continue;
//        cout << i << " th cluster ";
        n = _clusters[i].numberOfBonds();
        total_bonds += n;
        if(n > mx_n){
            mx_n = n;
        }
        sz = cluster_counts.size();
        if (n >= sz){
            cluster_counts.resize(n+1);
        }
//        cout << cluster_counts[n];
        ++cluster_counts[n]; // increase cluster count

//        cout << endl;
    }
    // since isolated bonds forms cluster of unit size and it might not be in the cluster
    cluster_counts[1] += maxBonds() - total_bonds;
//    cout << "before returning " << mx_n << endl;
    return cluster_counts;
}

void SitePercolation_ps_v10::init(bool gen_random) {
    if(gen_random) setRandomState(0, gen_random);

    // to improve performence
    number_of_sites_to_span.reserve(maxSites());
    number_of_bonds_to_span.reserve(maxSites());

    _top_edge.reserve(length());
    _bottom_edge.reserve(length());
    _left_edge.reserve(length());
    _right_edge.reserve(length());

//    initialize_index_sequence();
    initialize_cluster();
    randomize_v2();  // randomize
    number_of_clusters = maxBonds();
}

value_type SitePercolation_ps_v10::manageClusters(
        const std::vector<Index>& sites, const std::vector<BondIndex> &bonds) {
//    return merge_cluster_v2(sites, bonds);
    return merge_cluster_v3(sites, bonds);
}

/**
 *
 * @param site_a : root site or the currect site
 * @param clstr_b : the cluster whose elements should be relabeled.
 *      this cluster should not contain sites or bonds that comes from more than once cluster at a time
 * @param bond_pos : from which position on the clstr_b bonds should be relabeled
 * @param site_pos : from which position on the clstr_b sites should be relabeled
 */
void
SitePercolation_ps_v10::relabel_cluster(
        Index site_a, const Cluster_v3 &clstr_b, size_t bond_pos, size_t site_pos)
{

    int id_a = _lattice.getGroupID(site_a);
    int id_b = clstr_b.get_ID();

    const vector<Index>     & sites = clstr_b.getSiteIndices();
    const vector<BondIndex> & bonds = clstr_b.getBondIndices();
//    cout << "relabel_cluster : base id " << id_a << " {";
    for(size_t i{bond_pos}; i < bonds.size(); ++i){
//        cout << bonds[i] << ",";
        _lattice.setGroupID(bonds[i], id_a);
    }
//    cout << "}" << endl;

    // get four site_b of site_a
    vector<Index> sites_neighbor_a = _lattice.get_neighbor_sites(site_a);
    Index site_b;
    IndexRelative relative_index_b_after;
//    bool flag{false};
    // find which site_b has id_a of clstr_b
    for(auto n: sites_neighbor_a){
        if(id_b == _lattice.getSite(n).get_groupID()){
            // checking id_a equality is enough. since id_a is the id_a of the active site already.
            relative_index_b_after = getRelativeIndex(site_a, n);
            site_b = n;
//            cout << "neighbor  of" << site_a << " is " << site_b << endl;
//            flag = true;
            break;
        }
    }
//    if(!flag){
//        cout << "No neibhgor found! : line " << __LINE__ << endl;
//    }
    IndexRelative relative_site_a = _lattice.getRelativeIndex(site_a);

    // with this delta_a and delta_y find the relative index of site_b while relative index of site_a is known
    IndexRelative relative_site_b_before = _lattice.getRelativeIndex(site_b);

    int delta_x_ab = relative_index_b_after.x_ - relative_site_b_before.x_;
    int delta_y_ab = relative_index_b_after.y_ - relative_site_b_before.y_;
//    cout << relative_index_b_after << " - " << relative_site_b_before << " = delta_x, delta_y = " << delta_x_ab << ", " << delta_y_ab << endl;
//    relabel_sites(sites, id_a, delta_x_ab, delta_y_ab);

    int x, y;
    Index a;
    for (value_type i = site_pos; i < sites.size(); ++i) {
        a = sites[i];
        _lattice.setGroupID(a, id_a);
        relative_site_a = _lattice.getRelativeIndex(a);
        x = relative_site_a.x_ + delta_x_ab;
        y = relative_site_a.y_ + delta_y_ab;
//        _lattice.getSite(a).relativeIndex(x, y);
        _lattice.setRelativeIndex(a, {x,y});
    }
}

long double SitePercolation_ps_v10::entropy_v2() {
    long double H{};
    double number_of_cluster_with_size_one = maxBonds() - _bonds_in_cluster_with_size_two_or_more;
//    cout << " _bonds_in_cluster_with_size_two_or_more " << _bonds_in_cluster_with_size_two_or_more << " : line " << __LINE__ << endl;
    double mu = 1.0/double(maxBonds());
    H += number_of_cluster_with_size_one * logl(mu) * mu;
    H *= -1;
    _entropy_current =  _entropy + H;
    return _entropy_current;
}

long double SitePercolation_ps_v10::entropy_v1() {
    double nob, mu;
    long double H{};
    for(size_t i{}; i < _clusters.size(); ++i){
        nob = _clusters[i].numberOfBonds();
        if(nob == 0) continue;
        mu = nob/maxBonds();
        H += mu * logl(mu);
    }
    _entropy_current = -H;
    return _entropy_current;
}

void SitePercolation_ps_v10::viewRandomizedIndices() {
    cout << "randomized indices {" << endl;
    for(size_t i{}; i < randomized_index.size(); ++i){
        cout << i << " [" << randomized_index[i] << "] -> " << site_indices[randomized_index[i]] << endl;
    }
    cout << "}" << endl;
}

//void SitePercolation_ps_v10::track_largestCluster(value_type base) {
//    // calculating number of bonds in the largest cluster // by cluster index
//    // checking number of bonds
//    if(_clusters[base].numberOfBonds() > _number_of_bonds_in_the_largest_cluster){
//        _number_of_bonds_in_the_largest_cluster = _clusters[base].numberOfBonds();
//    }
//    if(_clusters[base].numberOfSites() > _number_of_sites_in_the_largest_cluster){
//        _number_of_sites_in_the_largest_cluster = _clusters[base].numberOfSites();
//    }
//}



