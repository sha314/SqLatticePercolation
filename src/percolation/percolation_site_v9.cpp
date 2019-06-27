//
// Created by shahnoor on 9/4/18.
//


#include <cstdlib>
#include <climits>
#include <unordered_set>
#include <mutex>

#include "percolation.h"

#include "../util/printer.h"
#include "site_position.h"
#include "../lattice/connector.h"
#include "../index/delta.h"
#include <omp.h>
#include <thread>
#include <algorithm>

#include "../util/checking.h"
#include "../util/time_tracking.h"

using namespace std;



/**
 *
 * @param length       : length of the lattice
 * @param impure_sites : number of impure sites. cannot be greater than length*length
 */
SitePercolation_ps_v9::SitePercolation_ps_v9(value_type length, bool periodicity)
        :SqLatticePercolation(length)
{
    std::cout << "Constructing SitePercolation_ps_v9 object : line " << __LINE__ << endl;
    SqLatticePercolation::set_type('s');

    _periodicity = periodicity;
    _index_sequence_position = 0;


    min_index = 0;
    max_index = length - 1;

    index_sequence.resize(maxSites());
    randomized_index.resize(maxSites());
    _max_iteration_limit = maxSites();

    initialize();

}


/**
 *
 */
void SitePercolation_ps_v9::initialize() {
    _lattice = SqLattice(length(), true, false, false, true);   // since it is a site percolation all bonds will be activated by default
    initialize_index_sequence();

    // to improve performence
    number_of_sites_to_span.reserve(maxSites());
    number_of_bonds_to_span.reserve(maxSites());

    _top_edge.reserve(length());
    _bottom_edge.reserve(length());
    _left_edge.reserve(length());
    _right_edge.reserve(length());

//    randomized_index_sequence = bond_index_sequence;
}


/**
 * Called only once when the object is constructed for the first time
 */
void SitePercolation_ps_v9::initialize_index_sequence() {
    value_type m{}, n{};
    for (value_type i{}; i != index_sequence.size(); ++i) {
        randomized_index[i] = i;
        index_sequence[i] = Index(m, n);
        ++n;
        if (n == length()) {
            n = 0;
            ++m;
        }
    }
    //for (value_type i{}; i != bond_index_sequence.size(); ++i) {cout << bond_index_sequence[i] << endl;}
}


/**
 * Reset all calculated values and then call initiate()
 * to initiallize for reuse
 *
 * caution -> it does not erase _calculation_flags, for it will be used for calculation purposes
 */
void SitePercolation_ps_v9::reset() {
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

    initialize();
    randomize_v2();
    time_relabel = 0;
    _total_relabeling = 0;
}


/**
 * Randomize the indices
 */
void SitePercolation_ps_v9::randomize_v2(){

    std::shuffle(randomized_index.begin(), randomized_index.end(), _random);
//    cout << "Index sequence : " << randomized_index_sequence << endl;
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
 * @param found_index_set : does not contain the pase index
 * @param base : the base index of cluster
 */
void SitePercolation_ps_v9::subtract_entropy_for_bond(const set<value_type> &found_index, int base){
    double nob, mu_bond, H{};
    if(base >= 0){
        nob = _clusters[base].numberOfBonds();
        mu_bond = nob / maxBonds();
        H += log(mu_bond) * mu_bond;
    }
    for(auto x : found_index){
        nob = _clusters[x].numberOfBonds();
        mu_bond = nob / maxBonds();
        H += log(mu_bond) * mu_bond;
    }
    _entropy -= -H;
}



/**
 * Must be called after merging the clusters
 * Cluster length is measured by bonds
 * @param index
 */
void SitePercolation_ps_v9::add_entropy_for_bond(value_type index){
    double nob = _clusters[index].numberOfBonds();
    double mu_bond = nob / maxBonds();
    double H = log(mu_bond) * mu_bond;
    _entropy += -H;
}



/**
 * Condition: must be called each time a site is placed
 */
void SitePercolation_ps_v9::track_numberOfBondsInLargestCluster() {

    // calculating number of bonds in the largest cluster // by cluster index
    // checking number of bonds
    if(_clusters[_index_last_modified_cluster].numberOfBonds() > _number_of_bonds_in_the_largest_cluster){
        _number_of_bonds_in_the_largest_cluster = _clusters[_index_last_modified_cluster].numberOfBonds();
    }

}

/**
 *
 */
void SitePercolation_ps_v9::track_numberOfSitesInLargestCluster(){

    // calculating number of bonds in the largest cluster // by cluster index
    // checking number of bonds
    if(_clusters[_index_last_modified_cluster].numberOfSites() > _number_of_sites_in_the_largest_cluster){
        _number_of_sites_in_the_largest_cluster = _clusters[_index_last_modified_cluster].numberOfSites();
    }
}


/**
 *
 * Find one row from _cluster to place 4 or less new bonds
 * Also remove the matched index values, because they will be inserted later.
 * This gives an advantage, i.e., you don't need to perform a checking.
 * todo takes so much time
 * @param hv_bonds
 * @return a set
 */
set<value_type>
SitePercolation_ps_v9::find_index_for_placing_new_bonds(const vector<Index> &neighbors) {
    set<value_type> found_index_set;    // use set to prevent repeated index
    for (auto n: neighbors) {
        int id = _lattice.getSite(n).get_groupID();
        if(id >=0) {
            found_index_set.insert(value_type(id));
        }
    }

    return found_index_set;
}


/**
 *
 * @param neighbors         :
 * @param found_index_set   : index of the clusters that will be merged together.
 *                            Does not contain the base cluster index or id.
 * @return                  : id of the base cluster
 */
int
SitePercolation_ps_v9::find_cluster_index_for_placing_new_bonds(
        const vector<Index> &neighbors, std::set<value_type> &found_index_set
){
    found_index_set.clear();
    value_type size{}, tmp{}, index, base{ULONG_MAX};
    int base_id{-1};
    int id;
    for (auto n: neighbors) {
        id = _lattice.getSite(n).get_groupID();
        if(id >=0) {
            index = value_type(id);
            tmp = _clusters[index].numberOfSites();
            if(tmp > size){
                size = tmp;
                base_id = id;
                base = index;
            }

            found_index_set.insert(index);

        }
    }
    found_index_set.erase(base);
    return base_id;
}



/**
 * Last placed site is added to a cluster. If this connects other clusters then merge all
 * cluster together to get one big cluster. All sites that are part of the other clusters
 * are relabled according to the id of the base cluster.
 * @param found_index_set : index of the clusters that are neighbors of the last placed site
 * @param hv_bonds        : bonds that connects the last placed site and its neighbors
 *                          and which are not part of any cluster of size larger than one
 * @param site            : last placed site
 * @return
 */
value_type SitePercolation_ps_v9::manage_clusters(
        const set<value_type> &found_index_set,
        vector<BondIndex> &hv_bonds,
        Index &site
)
{
    clock_t t = clock();
    vector<value_type> found_index(found_index_set.begin(), found_index_set.end());


    if (!found_index_set.empty()) {
        value_type base = found_index[0];
        _clusters[base].addSiteIndex(site);
        int id_base = _clusters[base].get_ID();
        vector<Index> neibhgors = _lattice.get_neighbor_sites(site);
        // find which of the neighbors are of id_base as the base cluster
        IndexRelative r;
        for(auto n: neibhgors){
            if(_lattice.getSite(n).get_groupID() == id_base){
                // find relative index with respect to this site
                r = getRelativeIndex(n, site);
                break; // since first time r is set running loop is doing no good
            }
        }

        // put_values_to_the_cluster new values in the 0-th found index
        _clusters[base].insert(hv_bonds);
        _lattice.getSite(site).relativeIndex(r);
        _lattice.getSite(site).set_groupID(id_base); // relabeling for 1 site

        // merge clusters with common values from all other cluster        // merge clusters with common values from all other cluster

        value_type ers{};
        for (value_type k{1}; k != found_index.size(); ++k) {
            ers = found_index[k];
            _total_relabeling += _clusters[ers].numberOfSites(); // only for debugging purposes
            // perform relabeling on the site_index_sequence
            relabel_sites_v5(site, _clusters[ers]);

            // store values of other found indices to the cluster
            _clusters[base].insert(_clusters[ers]);
            _cluster_count--; // reducing number of clusters
            _clusters[ers].clear(); // emptying the cluster
            
        }
        _index_last_modified_cluster = base;


    } else {
        // create new element for the cluster
        _clusters.push_back(Cluster_v3(_cluster_id));
        value_type _this_cluster_index = _clusters.size() -1;
        _lattice.getSite(site).set_groupID(_cluster_id); // relabeling for 1 site
        _cluster_count++; // increasing number of clusters
        _cluster_id++;
        _clusters.back().insert(hv_bonds);
        _clusters[_this_cluster_index].addSiteIndex(site);
        _index_last_modified_cluster = _this_cluster_index;   // last cluster is the place where new bonds are placed

    }

    return _index_last_modified_cluster;
}


/**
 * Last placed site is added to a cluster. If this connects other clusters then merge all
 * cluster together to get one big cluster. All sites that are part of the other clusters
 * are relabled according to the id of the base cluster.
 * @param found_index_set : index of the clusters that are neighbors of the last placed site
 * @param hv_bonds        : bonds that connects the last placed site and its neighbors
 *                          and which are not part of any cluster of size larger than one
 * @param site            : last placed site
 * @param base_id         : id of the base cluster
 * @return
 */
value_type SitePercolation_ps_v9::manage_clusters(
        const set<value_type> &found_index_set,
        vector<BondIndex> &hv_bonds,
        Index &site,
        int base_id
)
{


    if (base_id != -1) {
        value_type base = value_type(base_id); // converting here
        _clusters[base].addSiteIndex(site);
        int id_base = _clusters[base].get_ID();
        vector<Index> neibhgors = _lattice.get_neighbor_sites(site);
        // find which of the neighbors are of id_base as the base cluster
        IndexRelative r;
        for(auto n: neibhgors){
            if(_lattice.getSite(n).get_groupID() == id_base){
                // find relative index with respect to this site
                r = getRelativeIndex(n, site);
                break; // since first time r is set running loop is doing no good
            }
        }

        // put_values_to_the_cluster new values in the 0-th found index
        _clusters[base].insert(hv_bonds);
        _lattice.getSite(site).relativeIndex(r);
        _lattice.getSite(site).set_groupID(id_base); // relabeling for 1 site

        // merge clusters with common values from all other cluster        // merge clusters with common values from all other cluster


        for(value_type ers: found_index_set){

            _total_relabeling += _clusters[ers].numberOfSites(); // only for debugging purposes
            // perform relabeling on the site_index_sequence
            relabel_sites_v5(site, _clusters[ers]);

            // store values of other found indices to the cluster
            _clusters[base].insert_v2(_clusters[ers]);
            _cluster_count--; // reducing number of clusters
            _clusters[ers].clear(); // emptying the cluster

        }
        _index_last_modified_cluster = base;


    } else {
        // create new element for the cluster
        _clusters.push_back(Cluster_v3(_cluster_id));
        value_type _this_cluster_index = _clusters.size() -1;
        _lattice.getSite(site).set_groupID(_cluster_id); // relabeling for 1 site
        _cluster_count++; // increasing number of clusters
        _cluster_id++;
        _clusters.back().insert(hv_bonds);
        _clusters[_this_cluster_index].addSiteIndex(site);
        _index_last_modified_cluster = _this_cluster_index;   // last cluster is the place where new bonds are placed

    }
    return _index_last_modified_cluster;
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
value_type SitePercolation_ps_v9::relabel(value_type index_1, value_type index_2){
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
//IndexRelative SitePercolation_ps_v9::getRelativeIndex(Index root, Index site_new){
////    cout << "Entry \"SitePercolation_ps_v9::getRelativeIndex\" : line " << __LINE__ << endl;
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
void SitePercolation_ps_v9::connection_v1(Index site, vector<Index> &neighbors, vector<BondIndex> &bonds)
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
void SitePercolation_ps_v9::connection_v2(Index site, vector<Index> &site_neighbor, vector<BondIndex> &bond_neighbor)
{

    site_neighbor.clear();
    bond_neighbor.clear();
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
 * This function must be called after a site is placed in the
 * newest cluster
 * and an id is given from the lattice according to the cluster they are on.
 * @param site
 */
void SitePercolation_ps_v9::save_index_if_in_boundary_v2(const Index& site){
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
bool SitePercolation_ps_v9::check_if_id_matches(Index site, const vector<Index> &edge){
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
bool SitePercolation_ps_v9::check_if_id_matches_and_erase(Index site, vector<Index> &edge){
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
bool SitePercolation_ps_v9::occupy() {
    if(_index_sequence_position >= maxSites()){
        return false;
    }


    Index site = selectSite();

//    placeSite(site);
    placeSite_weighted(site);


    _occuption_probability = occupationProbability(); // for super class
    return true;
}

/***
 * Index of the selected site must be provided with the argument
 *
 * Wrapping and spanning index arrangement is enabled.
 * Entropy is calculated smoothly.
 * Entropy is measured by site and bond both.
 * @param current_site
 * @return
 */
value_type SitePercolation_ps_v9::placeSite(Index current_site) {
    // randomly choose a site
    if (_number_of_occupied_sites == maxSites()) {
        return ULONG_MAX;// unsigned long int maximum value
    }

    _last_placed_site = current_site;
//    cout << "placing site " << current_site << endl;

    _lattice.activate_site(current_site);

    ++_number_of_occupied_sites;


    // find the bonds for this site
    vector<BondIndex> bonds;
    vector<Index>     sites;

//    connection_v1(current_site, site_index_sequence, bonds);
    connection_v2(current_site, sites, bonds);
    _bonds_in_cluster_with_size_two_or_more += bonds.size();
    // find one of hv_bonds in _clusters and add ever other value to that place. then erase other position
    set<value_type> found_index_set = find_index_for_placing_new_bonds(sites);

//    cout << "Found indices " << found_index_set << endl;

    subtract_entropy_for_bond(found_index_set);  // tracking entropy change

    auto t0 = chrono::system_clock::now();
    value_type merged_cluster_index = manage_clusters(
            found_index_set, bonds, current_site
    );
    auto t1 = chrono::system_clock::now();
    time_relabel += chrono::duration<double>(t1 - t0).count();


    add_entropy_for_bond(merged_cluster_index); // tracking entropy change

    // running tracker
    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster
    track_numberOfSitesInLargestCluster();
    return merged_cluster_index;
}

/***
 * Index of the selected site must be provided with the argument
 *
 * Wrapping and spanning index arrangement is enabled.
 * Entropy is calculated smoothly.
 * Entropy is measured by site and bond both.
 * @param current_site
 * @return
 */
value_type SitePercolation_ps_v9::placeSite(
        Index current_site,
        vector<Index>& neighbor_sites,
        vector<BondIndex>& neighbor_bonds
) {

    if (_number_of_occupied_sites == maxSites()) {
        return ULONG_MAX;// unsigned long int maximum value
    }
    _bonds_in_cluster_with_size_two_or_more += neighbor_bonds.size();
    _last_placed_site = current_site;
    _lattice.activate_site(current_site);
    ++_number_of_occupied_sites;

    // find one of hv_bonds in _clusters and add ever other value to that place. then erase other position
    set<value_type> found_index_set = find_index_for_placing_new_bonds(neighbor_sites);

    subtract_entropy_for_bond(found_index_set);  // tracking entropy change
    value_type merged_cluster_index = manage_clusters(
            found_index_set, neighbor_bonds, current_site
    );
    add_entropy_for_bond(merged_cluster_index); // tracking entropy change
    // running tracker
    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster
    track_numberOfSitesInLargestCluster();
    return merged_cluster_index;
}

/***
 * Index of the selected site must be provided with the argument
 *
 * Wrapping and spanning index arrangement is enabled.
 * Entropy is calculated smoothly.
 * Entropy is measured by site and bond both.
 * @param current_site
 * @return
 */
value_type SitePercolation_ps_v9::placeSite_weighted(Index current_site) {
    // randomly choose a site
    if (_number_of_occupied_sites == maxSites()) {
        return ULONG_MAX;// unsigned long int maximum value
    }

    _last_placed_site = current_site;
//    cout << "placing site " << current_site << endl;

    _lattice.activate_site(current_site);

    ++_number_of_occupied_sites;


    // find the bonds for this site
    vector<BondIndex> bonds;
    vector<Index>     sites;

//    auto t0 = chrono::system_clock::now();

    connection_v2(current_site, sites, bonds);
    _bonds_in_cluster_with_size_two_or_more += bonds.size();
//    auto t1 = chrono::system_clock::now();
//    time_relabel += chrono::duration<double>(t1 - t0).count();

    // find one of hv_bonds in _clusters and add ever other value to that place. then erase other position
    set<value_type> found_index_set;
    int  base_id = find_cluster_index_for_placing_new_bonds(sites, found_index_set);
//    if (found_index_set.find(value_type(base_id)) != found_index_set.end()) {
//        cout << "base " << base_id << " and Found indices " << found_index_set << endl;
//        exit(1);
//    }

    subtract_entropy_for_bond(found_index_set, base_id);  // tracking entropy change
    value_type merged_cluster_index = manage_clusters(
            found_index_set, bonds, current_site, base_id
    );
    add_entropy_for_bond(merged_cluster_index); // tracking entropy change
    // running tracker
    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster
    track_numberOfSitesInLargestCluster();

    return merged_cluster_index;
}

/***
 * Index of the selected site must be provided with the argument
 *
 * Wrapping and spanning index arrangement is enabled.
 * Entropy is calculated smoothly.
 * Entropy is measured by site and bond both.
 * @param current_site
 * @return
 */
value_type SitePercolation_ps_v9::placeSite_weighted(
        Index current_site,
        vector<Index>& neighbor_sites,
        vector<BondIndex>& neighbor_bonds
) {
    // randomly choose a site
    if (_number_of_occupied_sites == maxSites()) {
        return ULONG_MAX;// unsigned long int maximum value
    }

    _bonds_in_cluster_with_size_two_or_more += neighbor_bonds.size();
    _last_placed_site = current_site;
//    cout << "placing site " << current_site << endl;

    _lattice.activate_site(current_site);

    ++_number_of_occupied_sites;


    // find one of hv_bonds in _clusters and add ever other value to that place. then erase other position
    set<value_type> found_index_set;
    int  base_id = find_cluster_index_for_placing_new_bonds(neighbor_sites, found_index_set);
//    if (found_index_set.find(value_type(base_id)) != found_index_set.end()) {
//        cout << "base " << base_id << " and Found indices " << found_index_set << endl;
//        exit(1);
//    }

    subtract_entropy_for_bond(found_index_set, base_id);  // tracking entropy change
    value_type merged_cluster_index = manage_clusters(
            found_index_set, neighbor_bonds, current_site, base_id
    );
    add_entropy_for_bond(merged_cluster_index); // tracking entropy change
    // running tracker
    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster
    track_numberOfSitesInLargestCluster();
    return merged_cluster_index;
}



/**
 *
 * @return
 */
Index SitePercolation_ps_v9::selectSite(){
//    Index current_site = randomized_index_sequence[_index_sequence_position]; // old
    value_type index = randomized_index[_index_sequence_position];
    Index current_site = index_sequence[index]; // new process
    ++_index_sequence_position;

    return current_site;
}


/***************************************************
 * View methods
 ****************************************/


/**
 *
 */
void SitePercolation_ps_v9::spanningIndices() const {
    cout << "Spanning Index : id" << endl;
    for(auto i: _spanning_sites){
        cout << "Index " << i << " : id " << _lattice.getSite(i).get_groupID() << endl;
    }
}

void SitePercolation_ps_v9::wrappingIndices() const {
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
bool SitePercolation_ps_v9::detectSpanning_v5(const Index& site) {
//    cout << "Entry -> detectSpanning_v4() : line " << __LINE__ << endl;
    if(_periodicity) {
        cout << "Cannot detect spanning if _periodicity if ON: line " << __LINE__ << endl;
        return false;
    }


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
bool SitePercolation_ps_v9::detectSpanning_v6(const Index& site) {
//    cout << "Entry -> detectSpanning_v4() : line " << __LINE__ << endl;
    if(_periodicity) {
        cout << "Cannot detect spanning if _periodicity if ON: line " << __LINE__ << endl;
        return false;
    }

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
void SitePercolation_ps_v9::scanEdges() {

}




/***********************************
 * Wrapping Detection
 **********************************/
/**
 * Wrapping is detected here using the last placed site
 * @return bool. True if wrapping occured.
 */
bool SitePercolation_ps_v9::detectWrapping() {
    return detect_wrapping_v1();
}

bool SitePercolation_ps_v9::detect_wrapping_v1() {
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
void SitePercolation_ps_v9::relabel_sites(const Cluster_v3& clstr, int id) {
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
void SitePercolation_ps_v9::relabel_sites_v4(Index site_a, const Cluster_v3& clstr_b) {
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
void SitePercolation_ps_v9::relabel_sites_v5(Index site_a, const Cluster_v3& clstr_b) {
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

    relabel_sites(sites, id_a, delta_x_ab, delta_y_ab);

}



/**
 * Relabels site and also reassign relative index to the relabeled sites
 *
 * @param site_a  : root index of the base cluster or last added site index of the base cluster
 * @param clstr_b : 2nd cluster, which to be merged withe the root
 * @param id : id to be used for relabeling sites
 */
void SitePercolation_ps_v9::relabel_sites_v6(Index site_a, const Cluster_v3& clstr_b, int id) {
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



void SitePercolation_ps_v9::relabel_sites(const vector<Index> &sites, int id_a, int delta_x_ab, int delta_y_ab)  {
    int x, y;
    Index a;
    IndexRelative relative_site__a;
    for (value_type i = 0; i < sites.size(); ++i) {
        a = sites[i];
        _lattice.getSite(a).set_groupID(id_a);
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
void SitePercolation_ps_v9::relabel_bonds(const Cluster_v3& clstr, int id) {
    vector<BondIndex> bonds = clstr.getBondIndices();
    for(auto a: bonds){
        _lattice.getBond(a).set_groupID(id);
    }
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
double SitePercolation_ps_v9::spanningProbability() const {
    // TODO
    cout << "Not implemented : line " << __LINE__ << endl;
    return 0;
}


/**
 * Entropy calculation is performed here. The fastest method possible.
 * Cluster size is measured by bond.
 * @return current entropy of the lattice
 */
double SitePercolation_ps_v9::entropy() {
    double H{};
    double number_of_cluster_with_size_one = maxBonds() - _bonds_in_cluster_with_size_two_or_more;
//    cout << " _bonds_in_cluster_with_size_two_or_more " << _bonds_in_cluster_with_size_two_or_more << " : line " << __LINE__ << endl;
    double mu = 1.0/double(maxBonds());
    H += number_of_cluster_with_size_one * log(mu) * mu;
    H *= -1;
    _entropy_current =  _entropy + H;
    return _entropy_current;
}



/**
 * number of bonds in the largest cluster / total number of bonds
 *
 * here cluster size is measured by number of bonds in the cluster
 * @return
 */
double SitePercolation_ps_v9::orderParameter() const{
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
double SitePercolation_ps_v9::orderParameter_v2() const{
//    double bond_num = _clusters[_index_largest_cluster].numberOfBonds() / double(_total_bonds);
    double tmp = _number_of_bonds_in_the_largest_cluster / double(maxBonds());
    return tmp;
}



/**
 * Finds out the number of bonds in the largest cluster by scanning all the clusters.
 * Best if no calculation is done for this previously
 * @return
 */
value_type SitePercolation_ps_v9::numberOfBondsInTheLargestCluster() {
    cout << "Inefficient. user *_v2() : line " << __LINE__ << endl;
    value_type  len{}, nob{};
    for(auto c: _clusters){
        nob = c.numberOfBonds();
        if (len < nob){
            len = nob;
        }
    }
    _number_of_bonds_in_the_largest_cluster = len;
    return len;
}


/**
 * Only applicable if the number of bonds in the largest cluster is calculated when occupying the lattice.
 * Significantly efficient than the previous version numberOfBondsInTheLargestCluster()
 * @return
 */
value_type SitePercolation_ps_v9::numberOfBondsInTheLargestCluster_v2() {
//    return _clusters[_index_largest_cluster].numberOfBonds();
    return _number_of_bonds_in_the_largest_cluster;
}


/**
 * Only applicable if the number of sites in the largest cluster is calculated when occupying the lattice.
 * Significantly efficient than the previous version numberOfBondsInTheLargestCluster()
 * @return
 */
value_type SitePercolation_ps_v9::numberOfBondsInTheSpanningCluster() {
    cout << "Not Implemented : line " << __LINE__ << endl;
    return _number_of_bonds_in_the_largest_cluster;
}


/**
 *
 * @return
 */
value_type SitePercolation_ps_v9::numberOfSitesInTheLargestCluster() {
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
double SitePercolation_ps_v9::numberOfSitesInTheSpanningClusters() {
    double nos{};
//    for(auto i: spanning_cluster_ids){
//        nos += _clusters[_cluster_index_from_id[i]].numberOfSites();
//    }

    int id{};
    for(auto i: _spanning_sites){
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
double SitePercolation_ps_v9::numberOfBondsInTheSpanningClusters() {
    double nos{};
//    for(auto i: spanning_cluster_ids){
//        nos += _clusters[_cluster_index_from_id[i]].numberOfBonds();
//    }
    int id{};
    if(_spanning_sites.size() > 1){
        cout << endl;
        for(auto i: _spanning_sites){

            id = _lattice.getSite(i).get_groupID();
            if(id >=0) {
//            cout << id << ", ";
                nos += _clusters[id].numberOfBonds();
            }
        }
        cout << endl;
    }
    else{
        id = _lattice.getSite(_spanning_sites.front()).get_groupID();
        if(id >= 0) {
//        cout << endl << "only one spanning cluster " << id << endl;
            nos = _clusters[id].numberOfBonds();
        }
    }

    return nos;
}



/**
 *
 * @return
 */
value_type SitePercolation_ps_v9::numberOfSitesInTheSpanningClusters_v2() {

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
value_type SitePercolation_ps_v9::numberOfBondsInTheSpanningClusters_v2() {
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
value_type SitePercolation_ps_v9::numberOfSitesInTheWrappingClusters(){
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
value_type SitePercolation_ps_v9::numberOfBondsInTheWrappingClusters(){
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
int SitePercolation_ps_v9::birthTimeOfSpanningCluster() const {
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
int SitePercolation_ps_v9::birthTimeOfACluster(int id) const {
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
value_type SitePercolation_ps_v9::box_counting(value_type delta) {
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
array<value_type, 2> SitePercolation_ps_v9::box_counting_v2(value_type delta) {
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
value_type SitePercolation_ps_v9::box_counting_spanning(value_type delta) {
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
bool SitePercolation_ps_v9::anyActiveSite(value_type row, value_type col, value_type delta) {
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
bool SitePercolation_ps_v9::anyActiveSpanningSite(value_type row, value_type col, value_type delta) {
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
value_type SitePercolation_ps_v9::count_number_of_active_site() {
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

std::string SitePercolation_ps_v9::getSignature() {
    string s = "sq_lattice_site_percolation";
    if(_periodicity)
        s += "_periodic_";
    else
        s += "_non_periodic_";
    return s;
}

/**
 *
 * @param filename
 * @param only_spanning
 */
void SitePercolation_ps_v9::writeVisualLatticeData(const string &filename, bool only_spanning) {
    std::ofstream fout(filename);
    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length()
                << ",\"signature\":\"" << getSignature() << "\""
                << ",\"x\":\"" << lastPlacedSite().column_ << "\""
                << ",\"y\":\"" << lastPlacedSite().row_ << "\""
                << "}" ;

    fout << "#" << header_info.str() << endl;
    fout << "#<x>,<y>,<color>" << endl;
    fout << "# color=0 -means-> unoccupied site" << endl;
    int id{-1};
    if(!_spanning_sites.empty()){
        id = _lattice.getSite(_spanning_sites.front()).get_groupID();
    }
    else if(!_wrapping_sites.empty()){
        id = _lattice.getSite(_wrapping_sites.front()).get_groupID();
    }

    if(only_spanning){
        if(id < 0){
            cerr << "id < 0 : line " << __LINE__ << endl;
        }
        vector<Index> sites = _clusters[id].getSiteIndices();
        for(auto s: sites){
            fout << s.column_ << ',' << s.row_ << ',' << id << endl;
        }
    }
    else {
        for (value_type y{}; y != length(); ++y) {
            for (value_type x{}; x != length(); ++x) {
                id = _lattice.getSite({y, x}).get_groupID();
                if(id != -1) {
                    fout << x << ',' << y << ',' << id << endl;
                }
            }
        }
    }
    fout.close();
}

/**
 *
 * @return an array : index of the array is the cluster size.
 *                    value of the array is the number of times that cluster appears
 */
const std::vector<double> SitePercolation_ps_v9::clusterSizeDistribution() const{
    vector<double> cluster_counts;
    size_t n, sz;
    cout << "number of max clusters " << _clusters.size() << endl;
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
    cout << "before returning " << mx_n << endl;
    return cluster_counts;
}




