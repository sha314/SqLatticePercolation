//
// Created by shahnoor on 10/31/2017.
//

#include "percolation.h"


#include <cstdlib>
#include <climits>
#include <unordered_set>

#include "percolation.h"

#include "../util/printer.h"
#include "site_position.h"


using namespace std;



/**
 *
 * @param length
 */
BondPercolation_pb_v1::BondPercolation_pb_v1(value_type length, bool periodicity)
        : SqLatticePercolation(length) {
    std::cout << "Constructing BondPercolation_pb_v1 object : line " << __LINE__ << endl;
    SqLatticePercolation::set_type('b');
    _periodicity = periodicity;
    _index_sequence_position = 0;
    _number_of_occupied_bonds = 0;
    _lattice = SqLattice(length, false, true, true, true);

    index_sequence.reserve(maxBonds());

    // there are 2*L*L cluster initially but only clusters with size larger than 1 should be counted
    _clusters = vector<Cluster>();
    _max_iteration_limit = maxBonds();
    initialize_index_sequence();
    randomized_index.resize(maxIterationLimit());
    for(value_type i{}; i < maxIterationLimit(); ++i){randomized_index[i] = i;}
    initialize();

//    initialize_indices();

    randomize_v2();  // randomize the untouched_site_indices
}

/**
 * Initialize the Class Object
 */
void BondPercolation_pb_v1::initialize() {

    // to improve performence
    number_of_sites_to_span.reserve(maxSites());
    number_of_bonds_to_span.reserve(maxSites());


//    _top_edge.reserve(length());
//    _bottom_edge.reserve(length());
//    _left_edge.reserve(length());
//    _right_edge.reserve(length());

}


/**
 * Called only once when the object is constructed for the first time
 */
void BondPercolation_pb_v1::initialize_index_sequence() {
    value_type row{}, col{};

    if (!_periodicity) {
        value_type limit = maxSites() - 1; // since if row==4 and col ==4 we don't want that bond
        for (value_type i{}; i < limit; ++i) {
            /*if (row == max_index && col == max_index){
                continue; // it is the last step. so decrease the loop limit by 1 will do the job
            }else */
            if (row == max_index) {
                index_sequence.push_back({BondType::Horizontal, row, col});
            } else if (col == max_index) {
                index_sequence.push_back({BondType::Vertical, row, col});
            } else {
                index_sequence.push_back({BondType::Horizontal, row, col});
                index_sequence.push_back({BondType::Vertical, row, col});
            }

            ++col;
            if (col == length()) {
                col = 0;
                ++row;
            }
        }
    } else {
        for (value_type i{}; i < maxSites(); ++i) {

            index_sequence.push_back({BondType::Horizontal, row, col});
            index_sequence.push_back({BondType::Vertical, row, col});

            ++col;
            if (col == length()) {
                col = 0;
                ++row;
            }
        }
    }
}

/**
 *
 */
void BondPercolation_pb_v1::initialize_indices() {
    indices = vector<value_type>(index_sequence.size());
    for (value_type i{}; i != indices.size(); ++i) {
        indices[i] = i; // assign index first
    }
}


/**
 * Randomize the site placing order
 * Takes 3.031000 sec for 1000 x 1000 sites
 */
void BondPercolation_pb_v1::randomize() {
    randomized_index_sequence = index_sequence;
    value_type len = randomized_index_sequence.size();
    value_type j{};
    BondIndex tmp;
    for (value_type i{}; i != len; ++i) {
        // select a j from the array. which must come from the ordered region
        j = i + std::rand() % (len - i);
        // perform the swapping with i-th and j-th value
        tmp = randomized_index_sequence[i];
        randomized_index_sequence[i] = randomized_index_sequence[j]; // todo problematic
        randomized_index_sequence[j] = tmp;
    }

//    cout << "Index sequence : " << randomized_index_sequence.size()
//         << " : " << randomized_index_sequence << endl;

}

/**
 * Randomize the site placing order
 * Takes 3.031000 sec for 1000 x 1000 sites
 */
void BondPercolation_pb_v1::randomize_v2() {

    std::shuffle(randomized_index.begin(), randomized_index.end(), _random_generator);
//    cout << "Index sequence : " << randomized_index_sequence << endl;
}


/**
 * Reset all calculated values and then call initiate()
 * to initiallize for reuse
 *
 * caution -> it does not erase _calculation_flags, for it will be used for calculation purposes
 */
void BondPercolation_pb_v1::reset() {
    SqLatticePercolation::reset();
    _total_number_of_active_bonds = 0;
    _index_sequence_position = 0;
    _number_of_occupied_bonds = 0;
    number_of_bonds_to_span.clear();
    _sites_required_to_min_span = 0;
    sites_in_cluster_with_size_greater_than_one = 0;
    _wrapping_indices.clear();

    initialize();

//    cout << "use suffle : line " << __LINE__ << endl;
    randomize_v2();
}


/**
 * Create a custome configuration for the lattice
 */
void BondPercolation_pb_v1::configure(std::vector<BondIndex> bond_indices) {
    cout << "Entry -> configure() : line " << endl;

}


///**
// *
// */
//void BondPercolation_pb_v1::calculate_occupation_probability() {
//    // occupation probability = number of site present / total number of site
//    _number_of_occupied_bonds.push_back(_total_number_of_active_bonds);
//}

/**
 *
 */
double BondPercolation_pb_v1::entropy_slow() {
    double mu{}; // probability for the cluster
    double H{}; // entropyDistribution
    double x;
    double count{};
    for (value_type i{}; i < _clusters.size(); ++i) {
        x = _clusters[i].numberOfSites();
        count += x;
        mu = x / double(maxSites());
//        cout << "mu " << mu << endl;
        H += mu * log(mu);
    }
    // for cluster with numberOfBonds 1
//    cout << sites_in_cluster_with_size_greater_than_one  << " == ? " << count << endl;
    double number_of_cluster_with_size_one = maxSites() - sites_in_cluster_with_size_greater_than_one;
//    cout << number_of_cluster_with_size_one << endl;
    mu = 1 / double(maxSites());
    H += number_of_cluster_with_size_one * mu * log(mu);
    H *= -1;    // since S = - mu * log(mu)
//    cout << "H = " << H << endl;
    _entropy_current = H;
    return H;
}



/**
 *
 * @param clstr
 * @param id
 */
void BondPercolation_pb_v1::relabel_sites(const Cluster &clstr, int id) {

    const vector<Index> sites = clstr.getSiteIndices();
    for (auto a: sites) {
        _lattice.getSite(a).set_groupID(id);
    }

}

/**
 *
 * @param clstr
 * @param id
 */
void BondPercolation_pb_v1::relabel_bonds(const Cluster &clstr, int id) {
    vector<BondIndex> bonds = clstr.getBondIndices();
    for (auto a: bonds) {
        _lattice.getBond(a).set_groupID(id);
    }
}

/**
 * Relabels site and also reassign relative index to the relabeled sites
  *
  * @param bond_a  : : last added bond index of the base cluster
  * @param clstr_b : 2nd cluster, which to be merged withe the root
  */
void BondPercolation_pb_v1::relabel_bonds_v1(BondIndex bond_a, const Cluster &clstr_b) {
    const vector<BondIndex> bonds = clstr_b.getBondIndices();
    int id_a = _lattice.getBond(bond_a).get_groupID();
    int id_b = clstr_b.get_ID();
    BondIndex b = clstr_b.getRootBond();

    // get four site_b of bond_a
    vector<BondIndex> sites_neighbor_a = _lattice.get_neighbor_bond_indices(bond_a);
    BondIndex bond_b;
    IndexRelative relative_index_b_after;
    bool flag{false};
    // find which site_b has id_a of clstr_b
    for(auto n: sites_neighbor_a){
        if(id_b == _lattice.getBond(n).get_groupID()){
            // checking id_a equality is enough. since id_a is the id_a of the active site already.
            relative_index_b_after = getRelativeIndex_v2(bond_a, n);
            bond_b = n;
//            cout << "neighbor  of" << bond_a << " is " << site_b << endl;
            flag = true;
            break; // todo ?

        }
    }

    if(!flag){
        cout << "No neibhgor found! : line " << __LINE__ << endl;
    }

    IndexRelative relative_site_a = _lattice.getBond(bond_a).relativeIndex();

    // with this delta_a and delta_y find the relative index of site_b while relative index of bond_a is known
    IndexRelative relative_site_b_before = _lattice.getBond(bond_b).relativeIndex();

    int delta_x_ab = relative_index_b_after.x_ - relative_site_b_before.x_;
    int delta_y_ab = relative_index_b_after.y_ - relative_site_b_before.y_;

//    cout << relative_index_b_after << " - " << relative_site_b_before << " = delta_x, delta_y = " << delta_x_ab << ", " << delta_y_ab << endl;

    relabel_bonds(bonds, id_a, delta_x_ab, delta_y_ab);

}



void BondPercolation_pb_v1::relabel_bonds(const vector<BondIndex> &sites, int id_a, int delta_x_ab, int delta_y_ab)  {
    int x, y;
    BondIndex a;
    IndexRelative relative_site__a;
    for (value_type i = 0; i < sites.size(); ++i) {
        a = sites[i];
        _lattice.getBond(a).set_groupID(id_a);
        relative_site__a = _lattice.getBond(a).relativeIndex();
        x = relative_site__a.x_ + delta_x_ab;
        y = relative_site__a.y_ + delta_y_ab;
        _lattice.getBond(a).relativeIndex(x, y);
    }
}

/**
 *
 * @param bond     : last placed bond
 */
void BondPercolation_pb_v1::relabel_cluster(BondIndex bond,  const vector<Index>& sites){

    Index site_a = {bond.row_, bond.column_};
    Index site_b;
    if(bond.horizontal()){
        site_b = {bond.row_, (bond.column_ + 1 ) % length()};
    }else{
        site_b = {(bond.row_ + 1 ) % length(), bond.column_};
    }
    // basically site_a and site_b must be added to cluster_b

    int id = _lattice.getBond(bond).get_groupID();
    int id_a = _lattice.getSite(site_a).get_groupID();
    int id_b = _lattice.getSite(site_b).get_groupID();

    Index site;
    if(sites.size() == 2){ // brand new cluster
        cout << "brand new cluster" << endl;
        site = sites[0];
        site_b = sites[1];
        _lattice.getSite(site).set_groupID(id);
    }
    if(sites.size() == 1){
        if(_lattice.getSite(sites[0]).get_groupID() == id_a){
            site = site_b;
            site_b = site_a;
        }
    }

    // relabeling sites
    int delta_x_ab {};
    int delta_y_ab {};

    // what it would be the relative index of the site
    IndexRelative relative_index_b_after = getRelativeIndex(site, site_b);
    // with this delta_a and delta_y find the relative index of site_b while relative index of site_a is known
    IndexRelative relative_site_b_before = _lattice.getSite(site_b).relativeIndex();
    delta_x_ab = relative_index_b_after.x_ - relative_site_b_before.x_;
    delta_y_ab = relative_index_b_after.y_ - relative_site_b_before.y_;


    IndexRelative relative_site_a = _lattice.getSite(site).relativeIndex();

//    cout << relative_index_b_after << " - " << relative_site_b_before << " = delta_x, delta_y = " << delta_x_ab << ", " << delta_y_ab << endl;

    int x, y;
    Index a;
    IndexRelative relative_site__a;
    for (value_type i{}; i < sites.size(); ++i) {
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
 * @param bond     : last placed bond
 * @param clstr_b  : given cluster
 * @param bond_pos : position of the bonds from where bonds will be relabeled
 *                   according to bond_pos-1 bond ?
 * @param site_pos : position of the sites from where sites will be relabeled
 *                   according to site_pos-1 site ?
 */
void BondPercolation_pb_v1::relabel_cluster(BondIndex bond, const Cluster& clstr_b, size_t bond_pos, size_t site_pos){
    const vector<Index> & sites = clstr_b.getSiteIndices();
    const vector<BondIndex> & bonds = clstr_b.getBondIndices();

    int id = _lattice.getBond(bond).get_groupID();
    Index site_a = {bond.row_, bond.column_};
    Index site_b;
    if(bond.horizontal()){
        site_b = {bond.row_, (bond.column_ + 1 ) % length()};
    }else{
        site_b = {(bond.row_ + 1 ) % length(), bond.column_};
    }

    int id_a = _lattice.getSite(site_a).get_groupID();
    int id_b = _lattice.getSite(site_b).get_groupID();

    Index base_site;
    if(_lattice.getSite(site_a).get_groupID() == id){
        base_site = site_a;
    }else{
        base_site = site_b;
        site_b = site_a;
    }

    /*// BEGIN debugging
    cout << "base site " << base_site << " id " << _lattice.getSite(base_site).get_groupID()
         << " relative index " << _lattice.getSite(base_site).relativeIndex() << endl;
    cout << "sites to be relabeled :" ;
    copy(sites.begin()+site_pos, sites.end(), std::ostream_iterator<Index>(cout));
    cout << endl;
    cout << "site_b " << site_b << endl;
    // END debugging*/

    // relabeling bonds
    for(size_t i{bond_pos}; i < bonds.size(); ++i){
        _lattice.getBond(bonds[i]).set_groupID(id);
    }

    // relabeling sites
    int delta_x_ab {};
    int delta_y_ab {};

    // what it would be the relative index of the site
    IndexRelative relative_site_b_before = _lattice.getSite(site_b).relativeIndex();
//    cout << "before " << relative_site_b_before << endl;
    IndexRelative relative_index_b_after = getRelativeIndex(base_site, site_b);
//    cout << "after " << relative_index_b_after << endl;
    // with this delta_a and delta_y find the relative index of site_b while relative index of site_a is known

    delta_x_ab = relative_index_b_after.x_ - relative_site_b_before.x_;
    delta_y_ab = relative_index_b_after.y_ - relative_site_b_before.y_;

    IndexRelative relative_site_a = _lattice.getSite(base_site).relativeIndex();
//    cout << "site_b " << site_b << endl;
//    cout << relative_index_b_after << " - " << relative_site_b_before
//         << " = (delta_x, delta_y) = (" << delta_x_ab << ", " << delta_y_ab << ")" << endl;

    int x, y;
    Index a;
    IndexRelative relative_site__a;
    for (value_type i{site_pos}; i < sites.size(); ++i) {
        a = sites[i];
        _lattice.getSite(a).set_groupID(id);
        relative_site__a = _lattice.getSite(a).relativeIndex();
//        cout << "index " << a << " relative index " << relative_site__a << " to ";
        x = relative_site__a.x_ + delta_x_ab;
        y = relative_site__a.y_ + delta_y_ab;
        _lattice.getSite(a).relativeIndex(x, y);
//        cout << _lattice.getSite(a).relativeIndex() << endl;
    }
}

/**
 *
 * @param last_bond
 * @param clstr_b
 */
void BondPercolation_pb_v1::relabel_v1(BondIndex last_bond, const Cluster &clstr_b) {
    Index site_a = {last_bond.row_, last_bond.column_};
    const vector<Index> &sites = clstr_b.getSiteIndices();
    int id_a = _lattice.getBond(last_bond).get_groupID();
    int id_b = clstr_b.get_ID();
    Index b = clstr_b.getRootSite();

    // get four site_b of site_a
    vector<Index> sites_neighbor_a = _lattice.get_neighbor_site_indices(site_a);
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
    relabel_bonds(clstr_b, id_a);
}

void BondPercolation_pb_v1::relabel_sites(const vector<Index> &sites, int id_a, int delta_x_ab, int delta_y_ab)  {
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
 * All site placing method in one place
 *
 * @return true if operation is successfull
 */
bool BondPercolation_pb_v1::occupy() {
    if (_index_sequence_position >= maxBonds()) {
        return false;
    }

    value_type v = placeBond_v1();
    _occuption_probability = occupationProbability(); // for super class
    auto s = SqLattice::get_neighbor_indices(length(), _last_placed_bond);
    if(_lattice.getSite(s[0]).get_groupID() != _lattice.getSite(s[1]).get_groupID()){
        cout << "relabeling is not perfect for site " << endl;
        for(auto a : s) {
            cout << a << " id " << _lattice.getSite(a).get_groupID() << endl;
        }
        return false; // to reminate iteration
    }
    return v != ULONG_MAX;
}


/**
 * Place one site at a time. Actually activates one site and 4 bonds at a time.
 * Also label and relabel each bond and site as per their corresponding cluster
 *
 * @return merged_cluster_index -> the index of _cluster
 *                                      where new site and bond indicis are placed
 *                                      or the merged cluster index
 *
 */
value_type BondPercolation_pb_v1::placeBond_v0() {

//    if (_index_sequence_position == randomized_index_sequence.size()) {
    if (_number_of_occupied_bonds == maxBonds()){
        return ULONG_MAX;// unsigned long int maximum value
    }
    value_type index = randomized_index[_index_sequence_position];
    BondIndex current_bond = index_sequence[index];
    _last_placed_bond = current_bond;
//    cout << "placing bond " << current_bond << " : " << __LINE__ << endl;

    _lattice.activateBond(current_bond);

    ++_number_of_occupied_bonds;
    ++_index_sequence_position;

    // find the bonds for this site
    vector<BondIndex> bonds;
    vector<Index> sites;

//    connection_v1(current_bond, sites, bonds); // problem
    connection_v2(current_bond, sites, bonds);
//    cout << "Found bonds : " << bonds << endl;
//    cout << "Found sites : " << sites << endl;

    sites_in_cluster_with_size_greater_than_one += sites.size();
    // find one of hv_bonds in _clusters and add ever other value to that place. then erase other position
    set<value_type> found_index_set = find_index_for_placing_new_bonds(bonds);


//    cout << "Found indices " << found_index_set << endl;

    subtract_entropy_for_site(found_index_set);  // tracking entropy change
    value_type merged_cluster_index = manage_clusters(
            found_index_set, sites, current_bond
    );
    add_entropy_for_site(merged_cluster_index); // tracking entropy change

    // running tracker
//    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster
    return merged_cluster_index;
}


/**
 * Place one site at a time. Actually activates one site and 4 bonds at a time.
 * Also label and relabel each bond and site as per their corresponding cluster
 *
 * @return merged_cluster_index -> the index of _cluster
 *                                      where new site and bond indicis are placed
 *                                      or the merged cluster index
 *
 */
value_type BondPercolation_pb_v1::placeBond_v1() {

//    if (_index_sequence_position == randomized_index_sequence.size()) {
    if (_number_of_occupied_bonds == maxBonds()){
        return ULONG_MAX;// unsigned long int maximum value
    }
    value_type index = randomized_index[_index_sequence_position];
    ++_index_sequence_position;
    BondIndex current_bond = index_sequence[index];
    _last_placed_bond = current_bond;
//    cout << "placing bond " << current_bond << " : " << __LINE__ << endl;

    _lattice.activateBond(current_bond);

    ++_number_of_occupied_bonds;

    // find the bonds for this site
    vector<BondIndex> bonds;
    vector<Index> sites;

//    connection_v1(current_bond, sites, bonds); // problem
    connection_v2(current_bond, sites, bonds);
//    cout << "Found bonds : " << bonds << endl;
//    cout << "Found sites : " << sites << endl;

    sites_in_cluster_with_size_greater_than_one += sites.size();
    // find one of hv_bonds in _clusters and add ever other value to that place. then erase other position
    set<value_type> found_index_set;
    int base_id = find_cluster_index_for_placing_new_bonds_v2(bonds, found_index_set);

//    if (found_index_set.find(value_type(base_id)) != found_index_set.end()) {
//        cout << "base " << base_id << " and Found indices " << found_index_set << endl;
//        exit(1);
//    }


    subtract_entropy_for_site(found_index_set);  // tracking entropy change
    value_type merged_cluster_index = manage_clusters(
            found_index_set, sites, current_bond, base_id
    );
    add_entropy_for_site(merged_cluster_index); // tracking entropy change

    // running tracker
//    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster
    return merged_cluster_index;
}


/**
 * For marking sites that is a part of cluster with size > 1
 * @param sites
 */
void BondPercolation_pb_v1::mark_sites(vector<Index> &sites) {
    for(value_type i{}; i < sites.size(); ++i){
        if (_lattice.getSite(sites[i]).get_groupID() >= 0) {
            // if id >= 0 then it is counted in the cluster so erase it
//            sites.erase(it); // commenting for now.
        } else {
            _lattice.getSite(sites[i]).set_groupID(0); // else set id == 0
        }
    }
//    for (auto it{sites.begin()}; it < sites.end(); ++it) {
//        if (_lattice.getSite(*it).get_groupID() >= 0) {
//            // if id >= 0 then it is counted in the cluster so erase it
////            sites.erase(it); // commenting for now.
//        } else {
//            _lattice.getSite(*it).set_groupID(0); // else set id == 0
//        }
//    }
}


/***********************
 * new methods
 *
 ********************************************/


///**
// * Functions that will give correct value
// *      entropy()
// *      entropy_v3()
// * @param found_index_set
// * @param sites
// * @param site
// * @return
// */
//value_type BondPercolation_pb_v1::manage_clusters(
//        const set<value_type> &found_index_set,
//        vector<Index> &sites,
//        BondIndex &bond
//) {
//    vector<value_type> found_index(found_index_set.begin(), found_index_set.end());
//    value_type merged_cluster_index{};
//
//    if (!found_index_set.empty()) {
//        value_type base = found_index[0];
//        _clusters[base].addBondIndex(bond);
//        _lattice.getBond(bond).set_groupID(_clusters[base].get_ID()); // relabeling for 1 site
//
//        // put_values_to_the_cluster new values in the 0-th found index
//        _clusters[base].insert(sites); // todo
//        // merge clusters with common values from all other cluster
//        int tmp_id;
//        value_type ers;
//        for (value_type k{1}; k != found_index.size(); ++k) {
//            tmp_id = _clusters[base].get_ID();
//
//            // erase the redundant cluster ids
//            ers = found_index[k];
//                        
//            // perform relabeling on the sites
//            relabel_bonds(_clusters[ers], tmp_id);
//
//            // store values of other found indices to the cluster
//            _clusters[base].insert(_clusters[ers]);
//
//            _clusters[ers].clear();
//        }
//        merged_cluster_index = base;
//
//    } else {
//        // create new element for the cluster
//        _clusters.push_back(Cluster(_cluster_id));
//        merged_cluster_index = _clusters.size() - 1;  // this new cluster index
////        _cluster_index_from_id[_cluster_id] = _clusters.size() - 1; // keeps track of cluster id and cluster index
//        _cluster_index_from_id.insert(_cluster_id); // new version
//        _lattice.getBond(bond).set_groupID(_cluster_id); // relabeling for 1 site
//        _cluster_id++;  // increase the cluster id for next round
//        _clusters.back().insert(sites);
//        _clusters[merged_cluster_index].addBondIndex(bond);
//
//    }
//
//
//    // data for short cut calculation
//    _index_last_modified_cluster = merged_cluster_index;
//
//    return merged_cluster_index;
//}

/**
 * Relative index of site_new with respect to root
 * @param root
 * @param site_new
 * @return
 */
IndexRelative
BondPercolation_pb_v1::getRelativeIndex(BondIndex root_bond, BondIndex bond_new){
    Index root = {root_bond.row_, root_bond.column_};
    Index site_new = {bond_new.row_, bond_new.column_};

//    cout << "Entry \"SitePercolation_ps_v8::getRelativeIndex\" : line " << __LINE__ << endl;
    int delta_x = -int(root.column_) + int(site_new.column_); // if +1 then root is on the right ??
    int delta_y = int(root.row_) - int(site_new.row_); // if +1 then root is on the top ??


    // normalizing delta_x
    if(delta_x > 1){
        delta_x /= -delta_x;
    }
    else if(delta_x < -1){
        delta_x /= delta_x;
    }

    // normalizing delta_y
    if(delta_y > 1){
        delta_y /= -delta_y;
    }else if(delta_y < -1){
        delta_y /= delta_y;
    }

    IndexRelative indexRelative_root = _lattice.getSite(root).relativeIndex();
//    cout << "Relative index of root " << indexRelative_root << endl;
//    cout << "Delta x,y " << delta_x << ", " << delta_y << endl;
    IndexRelative r =  {indexRelative_root.x_ + delta_x, indexRelative_root.y_ + delta_y};
//    cout << "Relative index of site_new " << r << endl;
    return r;
}

/**
 * Relative index of site_new with respect to root
 * @param root
 * @param site_new
 * @return
 */
IndexRelative BondPercolation_pb_v1::getRelativeIndex(Index root, Index site_new){
//    cout << "Entry \"SitePercolation_ps_v8::getRelativeIndex\" : line " << __LINE__ << endl;
    int delta_x = -int(root.column_) + int(site_new.column_); // if +1 then root is on the right ??
    int delta_y = int(root.row_) - int(site_new.row_); // if +1 then root is on the top ??


    // normalizing delta_x
    if(delta_x > 1){
        delta_x /= -delta_x;
    }
    else if(delta_x < -1){
        delta_x /= delta_x;
    }

    // normalizing delta_y
    if(delta_y > 1){
        delta_y /= -delta_y;
    }else if(delta_y < -1){
        delta_y /= delta_y;
    }

    IndexRelative indexRelative_root = _lattice.getSite(root).relativeIndex();
//    cout << "root " << root << " Relative index " << indexRelative_root << endl;
//    cout << "Delta x,y " << delta_x << ", " << delta_y << endl;
    IndexRelative r =  {indexRelative_root.x_ + delta_x, indexRelative_root.y_ + delta_y};
//    cout << "site_new " << site_new << " Relative index " << r << endl;
//    cout << "Exit \"SitePercolation_ps_v8::getRelativeIndex\" : line " << __LINE__ << endl;
    return r;
}


/**
 * Relative index of new_bond with respect to root.
 * Considering sites as bases.
 * relative index of a bond is calculated with respect to
 * the index of the site it is connected to.
 * @param root
 * @param site_new
 * @return
 */
IndexRelative
BondPercolation_pb_v1::getRelativeIndex_v2(BondIndex root, BondIndex bond_new){
    cout << "Entry getRelativeIndex_v2 ; line " << __LINE__ << endl;
    vector<Index> sites_root = _lattice.get_neighbor_indices(root);
    cout << "sites of " << root << " are " << sites_root << endl;
    vector<Index> sites_new = _lattice.get_neighbor_indices(bond_new);
    cout << "sites of " << bond_new << " are " << sites_new << endl;
//    cout << "Entry \"SitePercolation_ps_v8::getRelativeIndex\" : line " << __LINE__ << endl;
    // finding common site
    Index site_common;
    for(size_t i{}; i < sites_root.size(); ++i){
        for(size_t j{}; j < sites_new.size(); ++j){
            if(sites_root[i] == sites_new[j]){
                site_common = sites_new[j];
                sites_new.erase(sites_new.begin() + j);
                sites_root.erase(sites_root.begin() + i);
                break;
            }
        }
    }

    cout << "common site " << site_common << endl;
    cout << "sites of root " << sites_root << endl;
    cout << "sites of new " << sites_new << endl;

    int delta_x {}; // if +1 then root is on the right ??
    int delta_y {}; // if +1 then root is on the top ??

// root is horizontal bond
    if(root.horizontal()){
        if(bond_new.horizontal()){
            // both root and new bond is horizontal
            if(bond_new.column_ == site_common.column_){
                // if column of new bond equals column of common site then new bond is on the right side of root bond
                delta_x = +1;
            }else{
                // left of root
                delta_x = -1;
            }
            delta_y = 0;
        }
        else{
            // root bond is horizontal but new bond is vertical
            if(site_common.row_ == bond_new.row_){
                // bottom of root bond
                delta_y = -1;
            }else{
                // top of root bond
                delta_y = +1;
            }
            if(site_common.column_ == root.column_){
                // left of the root bond
                delta_x = -1;
            }else{
                // right
                delta_x = +1;
            }
        }

    }// root is vertical bond
    else{
        if(bond_new.horizontal()){
            // root bond is vertical but new bond is horizontal
            if(bond_new.column_ == site_common.column_){
                // if column of new bond equals column of common site then new bond is on the right side of root bond
                delta_x = +1;
            }else{
                // left of root
                delta_x = -1;
            }
            if(root.row_ == site_common.row_){
                // top side of the root bond
                delta_y = +1;
            }else{
                // bottom
                delta_y = -1;
            }

        }
        else{
            // both bonds are vertical
            if(site_common.row_ == bond_new.row_){
                // bottom of root bond
                delta_y = -1;
            }else{
                // top of root bond
                delta_y = +1;
            }
            delta_x = 0;
        }

    }

    cout << "got delta_x " << delta_x << " and delta_y " << delta_y << endl;

    IndexRelative indexRelative_root = _lattice.getBond(root).relativeIndex();
//    cout << "Relative index of root " << indexRelative_root << endl;
//    cout << "Delta x,y " << delta_x << ", " << delta_y << endl;
    IndexRelative r =  {indexRelative_root.x_ + delta_x, indexRelative_root.y_ + delta_y};
//    cout << "Relative index of site_new " << r << endl;
    return r;
}



/**
 * Functions that will give correct value.
 * Relabel sites and bonds and wrapping is detected using sites
 * @param found_index_set
 * @param sites
 * @param site
 * @return
 */
value_type BondPercolation_pb_v1::manage_clusters(
        const set<value_type> &found_index_set,
        vector<Index> &sites,
        BondIndex &bond
) {
//    cout << "on test line " << __LINE__ << endl;
    vector<value_type> found_index(found_index_set.begin(), found_index_set.end());
    value_type merged_cluster_index{};

    if (!found_index_set.empty()) {
        value_type base = found_index[0];
        BondIndex root = _clusters[base].getRootBond(); // root index of the base cluster
        int id_base = _clusters[base].get_ID();

        _clusters[base].addBondIndex(bond);

        _lattice.getBond(bond).set_groupID(id_base); // relabeling for 1 site

        // put_values_to_the_cluster new values in the 0-th found index
        size_t site_pos = _clusters[base].numberOfSites();
        size_t bond_pos = _clusters[base].numberOfBonds();
        _clusters[base].insert(sites);
        relabel_cluster(bond, _clusters[base], bond_pos, site_pos);
        // merge clusters with common values from all other cluster
        int tmp_id;
        value_type ers{};
        for (value_type k{1}; k != found_index.size(); ++k) {
            ers = found_index[k];
            _total_relabeling += _clusters[ers].numberOfSites();
            // store values of other found indices to the cluster
            bond_pos = _clusters[base].numberOfBonds();
            site_pos = _clusters[base].numberOfSites();


            _clusters[base].insert(_clusters[ers]);

            relabel_cluster(bond, _clusters[base], bond_pos, site_pos);
            // delete the merged cluster
            _clusters[ers].clear();

        }

        merged_cluster_index = base;

    } else {
        // create new element for the cluster

        relabel_new_sites_relative(sites, _cluster_id);

        _clusters.push_back(Cluster(_cluster_id));
        merged_cluster_index = _clusters.size() - 1;  // this new cluster index

        _cluster_index_from_id.insert(_cluster_id); // new version
        _lattice.getBond(bond).set_groupID(_cluster_id); // relabeling for 1 site
        _cluster_id++;  // increase the cluster id for next round
        _clusters.back().insert(sites);
        _clusters[merged_cluster_index].addBondIndex(bond);

    }


    // data for short cut calculation
    _index_last_modified_cluster = merged_cluster_index;

    return merged_cluster_index;
}

/**
 * Functions that will give correct value.
 * Relabel sites and bonds and wrapping is detected using sites
 * @param found_index_set
 * @param sites
 * @param site
 * @return
 */
value_type BondPercolation_pb_v1::manage_clusters(
        const set<value_type> &found_index_set,
        vector<Index> &sites,
        BondIndex &bond,
        int base_id
) {
//    cout << "on test line " << __LINE__ << endl;

    value_type merged_cluster_index{};

    if (base_id != -1) {
        value_type base = value_type(base_id);
        BondIndex root = _clusters[base].getRootBond(); // root index of the base cluster

        _clusters[base].addBondIndex(bond);

        _lattice.getBond(bond).set_groupID(base_id); // relabeling for 1 site

        // put_values_to_the_cluster new values in the 0-th found index
        size_t site_pos = _clusters[base].numberOfSites();
        size_t bond_pos = _clusters[base].numberOfBonds();
        _clusters[base].insert(sites);
        relabel_cluster(bond, _clusters[base], bond_pos, site_pos);
        // merge clusters with common values from all other cluster


        for(value_type ers : found_index_set){
            // store values of other found indices to the cluster
            bond_pos = _clusters[base].numberOfBonds();
            site_pos = _clusters[base].numberOfSites();
            _total_relabeling += _clusters[ers].numberOfSites(); // only for debugging purposes
            _clusters[base].insert_v2(_clusters[ers]);

            relabel_cluster(bond, _clusters[base], bond_pos, site_pos);
            // delete the merged cluster
            _clusters[ers].clear();

        }

        merged_cluster_index = base;

    } else {
        // create new element for the cluster

        relabel_new_sites_relative(sites, _cluster_id);

        _clusters.push_back(Cluster(_cluster_id));
        merged_cluster_index = _clusters.size() - 1;  // this new cluster index

        _cluster_index_from_id.insert(_cluster_id); // new version
        _lattice.getBond(bond).set_groupID(_cluster_id); // relabeling for 1 site
        _cluster_id++;  // increase the cluster id for next round
        _clusters[merged_cluster_index].insert(sites);
        _clusters[merged_cluster_index].addBondIndex(bond);

    }


    // data for short cut calculation
    _index_last_modified_cluster = merged_cluster_index;

    return merged_cluster_index;
}

/**
 * Applicable only when a new cluster is created
 * @param sites
 * @param id
 */
void BondPercolation_pb_v1::relabel_new_sites_relative(const vector<Index> &sites, int id) {
    for(size_t i{0}; i < sites.size(); ++i){
            _lattice.getSite(sites[i]).set_groupID(id);
            if(i > 0) {
                IndexRelative r = this->getRelativeIndex(sites[0], sites[i]);
                _lattice.getSite(sites[i]).relativeIndex(r);
            }
        }
}


/**
 * Take a bond index only if the corresponding site is active
 *
 * @param site
 * @param site_neighbor
 * @param bond_neighbor
 */
void
BondPercolation_pb_v1::connection_v1(BondIndex bond, vector<Index> &site_neighbor, vector<BondIndex> &bond_neighbor) {
    // for a bond there is two site. each site connects to maximum of 3 bonds.
    // select the site if one of the three bonds is active
    value_type next_column = (bond.column_ + 1) % length();
    value_type prev_column = (bond.column_ + 1 + length()) % length();
    value_type prev_row = (bond.row_ + 1 + length()) % length();
    value_type next_row = (bond.row_ + 1) % length();
    if (!_periodicity) {
        cout << "Testing ... ... .. .. .. . . . : line " << __LINE__ << endl;
        if (bond.horizontal()) {
            site_neighbor.resize(2);
            site_neighbor[0] = {bond.row_, bond.column_}; // left site
            site_neighbor[1] = {bond.row_, next_column}; // right site

            // without periodicity
            if (bond.row_ == min_index) { // top edge including corners
                if (bond.column_ == min_index) {
                    // upper left corner

                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, bond.row_, bond.column_}); // left end of bond
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, bond.row_, next_column}); // right end bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, next_column}); // right end bond
                    }

                } else if (bond.column_ == max_index) {
                    // upper right corner
                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, bond.row_, bond.column_}); // left end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column}); // left end bond
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, bond.row_, next_column}); // right end bond
                    }

                } else {
                    // top edge excluding corners

                    bond_neighbor.reserve(4);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, bond.row_, bond.column_}); // left end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column}); // left end bond
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, bond.row_, next_column}); // right end bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, next_column}); // right end bond
                    }
                }
            } else if (bond.row_ == max_index) { // bottom edge including corners

                if (bond.column_ == min_index) {
                    // lower left corner
                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_}); // left end of bond
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, prev_row, next_column}); // right end bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, next_column}); // right end bond
                    }

                } else if (bond.column_ == max_index) {
                    // lower right corner
                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_}); // left end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column}); // left end of bond
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, prev_row, next_column}); // right end bond
                    }

                } else {
                    // bottom edge excluding corners
                    //  bottom edge
                    bond_neighbor.reserve(4);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_}); // left end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column}); // left end of bond
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, prev_row, next_column}); // right end bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, next_column}); // right end bond
                    }

                }
            }
                /* site.x_ > min_index && site.x_ < max_index &&  is not possible anymore*/
            else if (bond.column_ == min_index) { // left edge excluding corners
                bond_neighbor.reserve(5);
                if (_lattice.getSite(site_neighbor[0]).isActive()) {
                    bond_neighbor.push_back({BondType::Vertical, bond.row_, bond.column_}); // left end of bond
                    bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_}); // left end of bond
                }
                if (_lattice.getSite(site_neighbor[1]).isActive()) {
                    bond_neighbor.push_back({BondType::Vertical, prev_row, next_column}); // right end bond
                    bond_neighbor.push_back({BondType::Vertical, bond.row_, next_column}); // right end bond
                    bond_neighbor.push_back({BondType::Horizontal, bond.row_, next_column}); // right end bond
                }

            } else if (bond.column_ == max_index) {
                // right edge no corners
                bond_neighbor.reserve(5);
                if (_lattice.getSite(site_neighbor[0]).isActive()) {
                    bond_neighbor.push_back({BondType::Vertical, bond.row_, bond.column_}); // left end of bond
                    bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_}); // left end of bond
                    bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column}); // left end of bond
                }
                if (_lattice.getSite(site_neighbor[1]).isActive()) {
                    bond_neighbor.push_back({BondType::Vertical, prev_row, next_column}); // right end bond
                    bond_neighbor.push_back({BondType::Vertical, bond.row_, next_column}); // right end bond
                }

            }

        } else if (bond.vertical()) {
            // without periodicity
            site_neighbor.resize(2);
            site_neighbor[0] = {bond.row_, bond.column_}; // top end of bond
            site_neighbor[1] = {next_row, bond.column_}; // bottom end of bond

            if (bond.row_ == min_index) { // top edge including corners
                if (bond.column_ == min_index) {
                    // upper left corner
                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        // top end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, bond.column_});
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        // bottom end of bond
                        bond_neighbor.push_back({BondType::Horizontal, next_row, bond.column_});
                        bond_neighbor.push_back({BondType::Vertical, next_row, bond.column_});
                    }

                } else if (bond.column_ == max_index) {
                    // upper right corner

                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        // top end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column});
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        // bottom end of bond
                        bond_neighbor.push_back({BondType::Horizontal, next_row, prev_column});
                        bond_neighbor.push_back({BondType::Vertical, next_row, bond.column_});
                    }

                } else {
                    // top edge excluding corners

                    bond_neighbor.reserve(5);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        // top end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, bond.column_});
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column});
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        // bottom end of bond
                        bond_neighbor.push_back({BondType::Horizontal, next_row, prev_column});
                        bond_neighbor.push_back({BondType::Horizontal, next_row, bond.column_});
                        bond_neighbor.push_back({BondType::Vertical, next_row, bond.column_});
                    }

                }
            } else if (bond.row_ == max_index) { // bottom edge including corners
                if (bond.column_ == min_index) {
                    // lower left corner
                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        // top end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, bond.column_});
                        bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_});
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        // bottom end of bond
                        bond_neighbor.push_back({BondType::Horizontal, next_row, bond.column_});

                    }

                } else if (bond.column_ == max_index) {
                    // lower right corner
                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        // top end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column});
                        bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_});
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        // bottom end of bond
                        bond_neighbor.push_back({BondType::Horizontal, next_row, prev_column});
                    }

                } else {
                    // bottom edge excluding corners
                    //  bottom edge
                    bond_neighbor.reserve(5);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        // top end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, bond.column_});
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column});
                        bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_});
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        // bottom end of bond
                        bond_neighbor.push_back({BondType::Horizontal, next_row, bond.column_});
                        bond_neighbor.push_back({BondType::Horizontal, next_row, prev_column});
                    }
                }
            }
                /* site.x_ > min_index && site.x_ < max_index &&  is not possible anymore*/
            else if (bond.column_ == min_index) {
                // left edge not in the corners
                bond_neighbor.reserve(4);
                if (_lattice.getSite(site_neighbor[0]).isActive()) {
                    // top end of bond
                    bond_neighbor.push_back({BondType::Horizontal, bond.row_, bond.column_});
                    bond_neighbor.push_back({BondType::Vertical, bond.row_, bond.column_});
                }
                if (_lattice.getSite(site_neighbor[1]).isActive()) {
                    // bottom end of bond
                    bond_neighbor.push_back({BondType::Horizontal, next_row, bond.column_});
                    bond_neighbor.push_back({BondType::Vertical, next_row, bond.column_});
                }

            } else if (bond.column_ == max_index) {
                // right edge no corners
                bond_neighbor.reserve(4);
                if (_lattice.getSite(site_neighbor[0]).isActive()) {
                    // top end of bond
                    bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column});
                    bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_});
                }
                if (_lattice.getSite(site_neighbor[1]).isActive()) {
                    // bottom end of bond
                    bond_neighbor.push_back({BondType::Vertical, next_row, bond.column_});
                    bond_neighbor.push_back({BondType::Horizontal, next_row, prev_column});
                }
            }
        }
        return;
    }
    // A bond can have 2 site neighbor and 6 bond neighbor
    site_neighbor = vector<Index>(2);
    site_neighbor[0] = {bond.row_, bond.column_};
    // as our BondIndex suggests neighbor[1] should be downward or rightward, we set it below.


    // now we have all the neighboring site
    // it's time to find the 6 neighboring bond_neighbor
    // without periodicity
    if (bond.horizontal()) {
        // increase column index for the right neighbor
        site_neighbor[1] = {bond.row_, next_column};

        bond_neighbor.reserve(6);
        if (_lattice.getSite(site_neighbor[0]).isActive()) {
            // left end of bond
            bond_neighbor.push_back({BondType::Vertical, bond.row_, bond.column_});
            bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_});
            bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column});
        }
        if (_lattice.getSite(site_neighbor[1]).isActive()) {
            // right end bond
            bond_neighbor.push_back({BondType::Vertical, prev_row, next_column});
            bond_neighbor.push_back({BondType::Vertical, bond.row_, next_column});
            bond_neighbor.push_back({BondType::Horizontal, bond.row_, next_column});
        }
    } else if (bond.vertical()) {
        // increase row index
        site_neighbor[1] = {next_row, bond.column_};

        bond_neighbor.reserve(6);
        if (_lattice.getSite(site_neighbor[0]).isActive()) {
            // top end of bond
            bond_neighbor.push_back({BondType::Horizontal, bond.row_, bond.column_});
            bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column});
            bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_});
        }
        if (_lattice.getSite(site_neighbor[1]).isActive()) {
            // bottom end of bond
            bond_neighbor.push_back({BondType::Horizontal, next_row, bond.column_});
            bond_neighbor.push_back({BondType::Horizontal, next_row, prev_column});
            bond_neighbor.push_back({BondType::Vertical, next_row, bond.column_});
        }

    }

}


/**
 * Take a bond index only if the corresponding site is active
 *
 * @param site
 * @param site_neighbor
 * @param bond_neighbor
 */
void
BondPercolation_pb_v1::connection_v2(BondIndex bond, vector<Index> &site_neighbor, vector<BondIndex> &bond_neighbor) {
    // for a bond there is two site. each site connects to maximum of 3 bonds.
    // select the site if one of the three bonds is active
    value_type next_column = (bond.column_ + 1) % length();
    value_type prev_column = (bond.column_ - 1 + length()) % length();
    value_type prev_row = (bond.row_ - 1 + length()) % length();
    value_type next_row = (bond.row_ + 1) % length();
    if (!_periodicity) {

        cout << "Testing ... ... .. .. .. . . . : line " << __LINE__ << endl;

        if (bond.horizontal()) {
            // without periodicity
            connection_2_horizontal_no_periodicity(bond, site_neighbor, bond_neighbor,
                                                   next_column, prev_column, prev_row);

        } else if (bond.vertical()) {
            // without periodicity
            connection_2_vertical_no_peridicity(bond, site_neighbor, bond_neighbor,
                                                prev_column, prev_row, next_row);
        }

    } else {
        // with periodicity
        connection_2_periodic(bond, site_neighbor, bond_neighbor,
                              prev_column, next_column, prev_row, next_row);

    }
}

/**
 * Can also be used for any bond that is not in the boundary.
 *
 * Take a site index only if the corresponding bond is active. // todo
 *
 * @param bond
 * @param site_neighbor
 * @param bond_neighbor
 * @param next_column
 * @param prev_column
 * @param prev_row
 * @param next_row
 */
void BondPercolation_pb_v1::connection_2_periodic(const BondIndex &bond, vector<Index> &site_neighbor,
                                                  vector<BondIndex> &bond_neighbor, value_type prev_column,
                                                  value_type next_column, value_type prev_row, value_type next_row) {
    // A bond can have 2 site neighbor and 6 bond neighbor
    site_neighbor.reserve(2);

//    vector<BondIndex> bond_tmp(6);
    // as our BondIndex suggests neighbor[1] should be downward or rightward, we set it below.


    // now we have all the neighboring site
    // it's time to find the 6 neighboring bond_neighbor
    // without periodicity
    if (bond.horizontal()) {
        // increase column index for the right neighbor

        bond_neighbor.resize(6);

        // left end of bond
        bond_neighbor[0] = {BondType::Vertical, bond.row_, bond.column_};
        bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_};
        bond_neighbor[2] = {BondType::Horizontal, bond.row_, prev_column};
        if(!(_lattice.getBond(bond_neighbor[0]).isActive() ||
           _lattice.getBond(bond_neighbor[1]).isActive() ||
           _lattice.getBond(bond_neighbor[2]).isActive())
                ){
            site_neighbor.push_back({bond.row_, bond.column_});
        }

        // right end bond
        bond_neighbor[3] = {BondType::Vertical, prev_row, next_column};
        bond_neighbor[4] = {BondType::Vertical, bond.row_, next_column};
        bond_neighbor[5] = {BondType::Horizontal, bond.row_, next_column};
        if(!(_lattice.getBond(bond_neighbor[3]).isActive() ||
           _lattice.getBond(bond_neighbor[4]).isActive() ||
           _lattice.getBond(bond_neighbor[5]).isActive())
                ){
            site_neighbor.push_back({bond.row_, next_column});
        }


    } else if (bond.vertical()) {
        // increase row index

        bond_neighbor.resize(6);

        // top end of bond
        bond_neighbor[0] = {BondType::Horizontal, bond.row_, bond.column_};
        bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column};
        bond_neighbor[2] = {BondType::Vertical, prev_row, bond.column_};
        if(!(_lattice.getBond(bond_neighbor[0]).isActive() ||
             _lattice.getBond(bond_neighbor[1]).isActive() ||
             _lattice.getBond(bond_neighbor[2]).isActive())
                ){
            site_neighbor.push_back({bond.row_, bond.column_});
        }

        // bottom end of bond
        bond_neighbor[3] = {BondType::Horizontal, next_row, bond.column_};
        bond_neighbor[4] = {BondType::Horizontal, next_row, prev_column};
        bond_neighbor[5] = {BondType::Vertical, next_row, bond.column_};
        if(!(_lattice.getBond(bond_neighbor[3]).isActive() ||
             _lattice.getBond(bond_neighbor[4]).isActive() ||
             _lattice.getBond(bond_neighbor[5]).isActive())
                ){
            site_neighbor.push_back({next_row, bond.column_});
        }

    }
}

/**
 * Non-Periodic Connection for vertical bond
 *
 * @param bond : vertical bond
 * @param site_neighbor : neighboring sites
 * @param bond_neighbor : neighboring bonds
 * @param prev_column
 * @param prev_row
 * @param next_row
 */
void BondPercolation_pb_v1::connection_2_vertical_no_peridicity(const BondIndex &bond, vector<Index> &site_neighbor,
                                                                vector<BondIndex> &bond_neighbor,
                                                                value_type prev_column, value_type prev_row,
                                                                value_type next_row) {
    bool top_site{false}, bottom_site{false};

    if (bond.row_ == min_index) { // top edge including corners
        if (bond.column_ == min_index) {
            // upper left corner
            bond_neighbor.resize(3);
            // top end of bond
            bond_neighbor[0] = {BondType::Horizontal, bond.row_, bond.column_};

            top_site = !_lattice.getBond(bond_neighbor[0]).isActive();

            // bottom end of bond
            bond_neighbor[1] = {BondType::Horizontal, next_row, bond.column_};
            bond_neighbor[2] = {BondType::Vertical, next_row, bond.column_};

            bottom_site = !(_lattice.getBond(bond_neighbor[1]).isActive() ||
                          _lattice.getBond(bond_neighbor[2]).isActive());

        } else if (bond.column_ == max_index) {
            // upper right corner

            bond_neighbor.resize(3);
            // top end of bond
            bond_neighbor[0] = {BondType::Horizontal, bond.row_, prev_column};
            top_site = !_lattice.getBond(bond_neighbor[0]).isActive();

            // bottom end of bond
            bond_neighbor[1] = {BondType::Horizontal, next_row, prev_column};
            bond_neighbor[2] = {BondType::Vertical, next_row, bond.column_};

            bottom_site = !(_lattice.getBond(bond_neighbor[1]).isActive() ||
                          _lattice.getBond(bond_neighbor[2]).isActive());
        } else {
            // top edge excluding corners

            bond_neighbor.resize(5);

            // top end of bond
            bond_neighbor[0] = {BondType::Horizontal, bond.row_, bond.column_};
            bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column};

            top_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                       _lattice.getBond(bond_neighbor[1]).isActive());

            // bottom end of bond
            bond_neighbor[2] = {BondType::Horizontal, next_row, prev_column};
            bond_neighbor[3] = {BondType::Horizontal, next_row, bond.column_};
            bond_neighbor[4] = {BondType::Vertical, next_row, bond.column_};

            bottom_site = !(_lattice.getBond(bond_neighbor[2]).isActive() ||
                          _lattice.getBond(bond_neighbor[3]).isActive() ||
                          _lattice.getBond(bond_neighbor[4]).isActive());

        }
    } else if (bond.row_ == max_index) { // bottom edge including corners
        if (bond.column_ == min_index) {
            // lower left corner
            bond_neighbor.resize(3);

            // top end of bond
            bond_neighbor[0] = {BondType::Horizontal, bond.row_, bond.column_};
            bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_};

            top_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                       _lattice.getBond(bond_neighbor[1]).isActive());

            // bottom end of bond
            bond_neighbor[2] = {BondType::Horizontal, next_row, bond.column_};

            bottom_site = !_lattice.getBond(bond_neighbor[2]).isActive();


        } else if (bond.column_ == max_index) {
            // lower right corner
            bond_neighbor.resize(3);

            // top end of bond
            bond_neighbor[0] = {BondType::Horizontal, bond.row_, prev_column};
            bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_};

            top_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                       _lattice.getBond(bond_neighbor[1]).isActive());

            // bottom end of bond
            bond_neighbor[2] = {BondType::Horizontal, next_row, prev_column};

            bottom_site = !_lattice.getBond(bond_neighbor[2]).isActive();

        } else {
            // bottom edge excluding corners
            //  bottom edge
            bond_neighbor.resize(5);

            // top end of bond
            bond_neighbor[0] = {BondType::Horizontal, bond.row_, bond.column_};
            bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column};
            bond_neighbor[2] = {BondType::Vertical, prev_row, bond.column_};

            top_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                       _lattice.getBond(bond_neighbor[1]).isActive() ||
                       _lattice.getBond(bond_neighbor[2]).isActive());

            // bottom end of bond
            bond_neighbor[3] = {BondType::Horizontal, next_row, bond.column_};
            bond_neighbor[4] = {BondType::Horizontal, next_row, prev_column};

            bottom_site = !(_lattice.getBond(bond_neighbor[3]).isActive() ||
                          _lattice.getBond(bond_neighbor[4]).isActive());
        }
    }
        /* site.x_ > min_index && site.x_ < max_index &&  is not possible anymore*/
    else if (bond.column_ == min_index) {
        // left edge not in the corners
        bond_neighbor.resize(4);
        // top end of bond
        bond_neighbor[0] = {BondType::Horizontal, bond.row_, bond.column_};
        bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_};

        top_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                   _lattice.getBond(bond_neighbor[1]).isActive());

        // bottom end of bond
        bond_neighbor[2] = {BondType::Horizontal, next_row, bond.column_};
        bond_neighbor[3] = {BondType::Vertical, next_row, bond.column_};

        bottom_site = !(_lattice.getBond(bond_neighbor[2]).isActive() ||
                      _lattice.getBond(bond_neighbor[3]).isActive());

    } else if (bond.column_ == max_index) {
        // right edge no corners
        bond_neighbor.resize(4);

        // top end of bond
        bond_neighbor[0] = {BondType::Horizontal, bond.row_, prev_column};
        bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_};

        top_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                   _lattice.getBond(bond_neighbor[1]).isActive());

        // bottom end of bond
        bond_neighbor[2] = {BondType::Vertical, next_row, bond.column_};
        bond_neighbor[3] = {BondType::Horizontal, next_row, prev_column};

        bottom_site = !(_lattice.getBond(bond_neighbor[2]).isActive() ||
                      _lattice.getBond(bond_neighbor[3]).isActive());
    }else{
        // one level inside the lattice
        bond_neighbor.resize(6);

        // top end of bond
        bond_neighbor[0] = {BondType::Horizontal, bond.row_, bond.column_};
        bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column};
        bond_neighbor[2] = {BondType::Vertical, prev_row, bond.column_};

        top_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                   _lattice.getBond(bond_neighbor[1]).isActive() ||
                   _lattice.getBond(bond_neighbor[2]).isActive());

        // bottom end of bond
        bond_neighbor[3] = {BondType::Horizontal, next_row, bond.column_};
        bond_neighbor[4] = {BondType::Horizontal, next_row, prev_column};
        bond_neighbor[5] = {BondType::Vertical, next_row, bond.column_};

        bottom_site = !(_lattice.getBond(bond_neighbor[3]).isActive() ||
                      _lattice.getBond(bond_neighbor[4]).isActive() ||
                      _lattice.getBond(bond_neighbor[5]).isActive());
    }

    site_neighbor.reserve(2);
    if (top_site) {
        site_neighbor.push_back({bond.row_, bond.column_});
    }
    if (bottom_site) {
        site_neighbor.push_back({next_row, bond.column_});
    }
}

void BondPercolation_pb_v1::connection_2_horizontal_no_periodicity(const BondIndex &bond,
                                                                   vector<Index> &site_neighbor,
                                                                   vector<BondIndex> &bond_neighbor,
                                                                   value_type next_column,
                                                                   value_type prev_column,
                                                                   value_type prev_row) {
    bool left_site{false}, right_site{false};

    // without periodicity
    if (bond.row_ == min_index) {
        // top edge including corners
        if (bond.column_ == min_index) {
            // upper left corner

            bond_neighbor.resize(3);
            // left end of bond
            bond_neighbor[0] = {BondType::Vertical, bond.row_, bond.column_};
            // right end bond
            bond_neighbor[1] = {BondType::Vertical, bond.row_, next_column};
            bond_neighbor[2] = {BondType::Horizontal, bond.row_, next_column};

            left_site = !_lattice.getBond(bond_neighbor[0]).isActive();

            right_site = !(_lattice.getBond(bond_neighbor[1]).isActive() ||
                         _lattice.getBond(bond_neighbor[2]).isActive());

        } else if (bond.column_ == max_index) {
            // upper right corner
            bond_neighbor.resize(3);
            // left end of bond
            bond_neighbor[0] = {BondType::Vertical, bond.row_, bond.column_};
            bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column};

            // right end bond
            bond_neighbor[3] = {BondType::Vertical, bond.row_, next_column};

            left_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                        _lattice.getBond(bond_neighbor[1]).isActive());

            right_site = !_lattice.getBond(bond_neighbor[3]).isActive();

        } else {
            // top edge excluding corners

            bond_neighbor.resize(4);
            // left end of bond
            bond_neighbor[0] = {BondType::Vertical, bond.row_, bond.column_};
            bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column}; // left end bond

            left_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                        _lattice.getBond(bond_neighbor[1]).isActive());

            // right end bond
            bond_neighbor[2] = {BondType::Vertical, bond.row_, next_column};
            bond_neighbor[3] = {BondType::Horizontal, bond.row_, next_column}; // right end bond

            right_site = !(_lattice.getBond(bond_neighbor[2]).isActive() ||
                         _lattice.getBond(bond_neighbor[3]).isActive());
        }
    } else if (bond.row_ == max_index) {
        // bottom edge including corners

        if (bond.column_ == min_index) {
            // lower left corner
            bond_neighbor.resize(3);
            // left end of bond
            bond_neighbor[0] = {BondType::Vertical, prev_row, bond.column_};
            left_site = !_lattice.getBond(bond_neighbor[0]).isActive();

            // right end bond
            bond_neighbor[1] = {BondType::Vertical, prev_row, next_column}; // right end bond
            bond_neighbor[2] = {BondType::Horizontal, bond.row_, next_column};

            right_site = !(_lattice.getBond(bond_neighbor[1]).isActive() ||
                         _lattice.getBond(bond_neighbor[2]).isActive());

        } else if (bond.column_ == max_index) {
            // lower right corner
            bond_neighbor.resize(3);
            // left end of bond
            bond_neighbor[0] = {BondType::Vertical, prev_row, bond.column_};
            bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column}; // left end of bond

            left_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                        _lattice.getBond(bond_neighbor[1]).isActive());

            // right end bond
            bond_neighbor[2] = {BondType::Vertical, prev_row, next_column};

            right_site = !_lattice.getBond(bond_neighbor[2]).isActive();

        } else {
            // bottom edge excluding corners
            //  bottom edge
            bond_neighbor.resize(4);

            bond_neighbor[0] = {BondType::Vertical, prev_row, bond.column_}; // left end of bond
            bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column}; // left end of bond

            left_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                        _lattice.getBond(bond_neighbor[1]).isActive());

            bond_neighbor[2] = {BondType::Vertical, prev_row, next_column}; // right end bond
            bond_neighbor[3] = {BondType::Horizontal, bond.row_, next_column}; // right end bond

            right_site = !(_lattice.getBond(bond_neighbor[2]).isActive() ||
                         _lattice.getBond(bond_neighbor[3]).isActive());


        }
    }
        /* site.x_ > min_index && site.x_ < max_index &&  is not possible anymore*/
    else if (bond.column_ == min_index) { // left edge excluding corners
        bond_neighbor.resize(5);

        bond_neighbor[0] = {BondType::Vertical, bond.row_, bond.column_}; // left end of bond
        bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_}; // left end of bond

        left_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                    _lattice.getBond(bond_neighbor[1]).isActive());

        bond_neighbor[2] = {BondType::Vertical, prev_row, next_column}; // right end bond
        bond_neighbor[3] = {BondType::Vertical, bond.row_, next_column}; // right end bond
        bond_neighbor[4] = {BondType::Horizontal, bond.row_, next_column}; // right end bond

        right_site = !(_lattice.getBond(bond_neighbor[2]).isActive() ||
                     _lattice.getBond(bond_neighbor[3]).isActive() ||
                     _lattice.getBond(bond_neighbor[4]).isActive());

    } else if (bond.column_ == max_index) {
        // right edge no corners
        bond_neighbor.resize(5);

        bond_neighbor[0] = {BondType::Vertical, bond.row_, bond.column_}; // left end of bond
        bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_}; // left end of bond
        bond_neighbor[2] = {BondType::Horizontal, bond.row_, prev_column}; // left end of bond

        left_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                    _lattice.getBond(bond_neighbor[1]).isActive() ||
                    _lattice.getBond(bond_neighbor[2]).isActive());

        bond_neighbor[3] = {BondType::Vertical, prev_row, next_column}; // right end bond
        bond_neighbor[4] = {BondType::Vertical, bond.row_, next_column}; // right end bond

        right_site = !(_lattice.getBond(bond_neighbor[3]).isActive() ||
                     _lattice.getBond(bond_neighbor[4]).isActive());

    }else{
        // one level inside the lattice

        bond_neighbor.resize(6);

        // left end of bond
        bond_neighbor[0] = {BondType::Vertical, bond.row_, bond.column_};
        bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_};
        bond_neighbor[2] = {BondType::Horizontal, bond.row_, prev_column};

        left_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                    _lattice.getBond(bond_neighbor[1]).isActive() ||
                    _lattice.getBond(bond_neighbor[2]).isActive());

        // right end bond
        bond_neighbor[3] = {BondType::Vertical, prev_row, next_column};
        bond_neighbor[4] = {BondType::Vertical, bond.row_, next_column};
        bond_neighbor[5] = {BondType::Horizontal, bond.row_, next_column};

        right_site = !(_lattice.getBond(bond_neighbor[3]).isActive() ||
                     _lattice.getBond(bond_neighbor[4]).isActive() ||
                     _lattice.getBond(bond_neighbor[5]).isActive());
    }

    site_neighbor.reserve(2);
    if (left_site) {
        site_neighbor.push_back({bond.row_, bond.column_});
    }
    if (right_site) {
        site_neighbor.push_back({bond.row_, next_column});
    }
}

/**
 * Takes 0.093 sec time when total time is 9.562 sec
 *
 * Find one row from _cluster to place 4 or less new bonds
 * Also remove the matched index values, because they will be inserted later.
 * This gives an advantage, i.e., you don't need to perform a checking.
 * todo takes so much time
 * @param hv_bonds
 * @return a set
 */
set<value_type>
BondPercolation_pb_v1::find_index_for_placing_new_bonds(const vector<BondIndex> &neighbors) {

    set<value_type> found_index_set;    // use set to prevent repeated index

    int id;
    for (auto n: neighbors) {
        id = _lattice.getBond(n).get_groupID();
        if(id >= 0){
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
BondPercolation_pb_v1::find_cluster_index_for_placing_new_bonds(const std::vector<BondIndex> &neighbors,
                                                                std::set<value_type> &found_indices){

    found_indices.clear();
    int base_id{-1};
    int id{0};
    value_type size_base{}, tmp, index;
    for (auto n: neighbors) {
        id = _lattice.getBond(n).get_groupID();
        if(id >= 0){
            index = value_type(id);
            tmp = _clusters[index].numberOfBonds();
            if(tmp > size_base){
                size_base = tmp;
                base_id = id;
                continue;
            }
            if(id == base_id){
                // already a base. harmful if included in found_index_set
                continue;
            }
            found_indices.insert(index);

        }
    }
    return base_id;
}


/**
 *
 * @param neighbors         :
 * @param found_index_set   : index of the clusters that will be merged together.
 *                            Does not contain the base cluster index or id.
 * @return                  : id of the base cluster
 */
int
BondPercolation_pb_v1::find_cluster_index_for_placing_new_bonds_v2(const std::vector<BondIndex> &neighbors,
                                                                std::set<value_type> &found_indices){

    found_indices.clear();
    int base_id{-1};
    int id{0};
    value_type size_base{}, tmp, index, base{};
    for (auto n: neighbors) {
        id = _lattice.getBond(n).get_groupID();
        if(id >= 0){
            index = value_type(id);
            tmp = _clusters[index].numberOfBonds();
            if(tmp > size_base){
                size_base = tmp;
                base_id = id;
                base = index;
            }
            found_indices.insert(index);
        }
    }

    // erase if it exists
    found_indices.erase(base);

    return base_id;
}



void BondPercolation_pb_v1::calculate_spanning_probability() {
    calculate_spanning_probability_by_largest_cluster();
}


/**
 * Use Group_ID to identify Bond and Site in the same cluster
 */
void BondPercolation_pb_v1::calculate_spanning_probability_by_largest_cluster() {
    if (debug_calculate_spanning_probability_by_largest_cluster_3) {
        cout << "calculate_spanning_probability_by_largest_cluster() : line " << __LINE__ << endl;
    }
    // find the largest cluster
    value_type l_largest_cluster{0}; // number of site of the largest cluster
    double b;
    for (auto a: _clusters) {
        b = a.numberOfSites();
        if (b < length()) // no percolation yet
            continue;
        if (b > l_largest_cluster) {
            l_largest_cluster = b;
        }
    }
//    double spanning_prob = l_largest_cluster/_number_of_occupied_sites;

    number_of_bonds_to_span.push_back(l_largest_cluster);

}


/**
 * Counts the number of active sites in the lattice
 */
value_type BondPercolation_pb_v1::count_number_of_active_site() {
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


//
///**
// * Give a site Index and it will give 4 bonds corresponding to that site
// * @param index
// * @return
// */
//vector<Index> BondPercolation_pb_v1::get_Sites_for_bonds(BondIndex index) {
//    if (debug_get_Bond_for_site) {
//        cout << "Entry -> get_Bond_for_site() : line " << __LINE__ << endl;
//    }
//    vector<Index> connected_sites(2);
//
//    if(index.horizontal()) {
//        // for horizontal bond, row remains the same
//        connected_sites[0] = _lattice.getSite({index.x_, index.column_}).set_ID();
//        auto c = (index.column_ + 1) % length();
//        connected_sites[1] = _lattice.getSite({index.x_, c}).set_ID();
//    }else{
//        // for vertical bond, column remains the same
//        connected_sites[0] = _lattice.getSite({index.x_, index.column_}).set_ID();
//        auto r = (index.x_ + 1) % length();
//        connected_sites[1] = _lattice.getSite({r, index.column_}).set_ID();
//    }
//
//
//    if (debug_get_Bond_for_site) {
//        cout << "Bonds for untouched site " << index << endl;
//        for (auto it = connected_sites.begin(); it != connected_sites.end(); ++it) {
//            cout << "set_ID " << (*it) << " = " << *it << endl;
//        }
//    }
//    return connected_sites;
//}



value_type BondPercolation_pb_v1::number_of_site_in_spanning_clusters(unordered_set<int> g_ids) {
    value_type nos{};   // number of sites
    for (auto b: g_ids) {
        for (auto a: _clusters) {
            if (a.get_ID() == b) {
                nos += a.numberOfSites();
            }
        }
    }
    return nos;
}


/***
 *
 * @return
 */
bool BondPercolation_pb_v1::detectSpanning() {
    if (debug_4_detectSpanning) {
        cout << "Entry -> detectSpanning() : line " << __LINE__ << endl;
    }

    // if any of the clusters does not have sites > length() --> no spanning
    value_type x{};
    for (auto a: _clusters) {
        if (x < a.numberOfSites())
            x = a.numberOfSites();
    }
    if (x < length()) {
        if (debug_4_detectSpanning) {
            cout << "not enough sites in any cluster to span : line " << __LINE__ << endl;
        }
        number_of_bonds_to_span.push_back(0);
        return false;
    }

    // vertical spanning
    std::unordered_set<int> first_row_group_ids;
    std::unordered_set<int> last_row_group_ids;
    // horizontal spanning
    std::unordered_set<int> first_column_group_ids;
    std::unordered_set<int> last_column_group_ids;

    for (value_type i{}; i != length(); ++i) {
        // for rows
        if (_lattice.getSite({0, i}).isActive()) {
            first_row_group_ids.insert(_lattice.getSite({0, i}).get_groupID());
        }
        if (_lattice.getSite({length() - 1, i}).isActive()) {
            last_row_group_ids.insert(_lattice.getSite({length() - 1, i}).get_groupID());
        }

        // for columns
        if (_lattice.getSite({i, 0}).isActive()) {
            first_column_group_ids.insert(_lattice.getSite({i, 0}).get_groupID());
        }
        if (_lattice.getSite({i, length() - 1}).isActive()) {
            last_column_group_ids.insert(_lattice.getSite({i, length() - 1}).get_groupID());
        }

    }

//    cout << '{';
//    for(auto a: first_row_group_ids){
//        cout << a << ',';
//    }
//    cout << '}' << endl;
//
//    cout << '{';
//    for(auto b: last_row_group_ids){
//        cout << b << ',';
//    }
//    cout << '}' << endl;

    std::vector<int> vertical_spanning_id;
    std::vector<int> horizontal_spanning_id;




    // to prevent repetitionn, use set
    unordered_set<int> spanning_ids;
    for (auto a: first_row_group_ids) {
        for (auto b: last_row_group_ids) {
            if (a == b) {
                spanning_ids.insert(a);
                vertical_spanning_id.push_back(a);
            }
        }
    }

    for (auto a: first_column_group_ids) {
        for (auto b: last_column_group_ids) {
            if (a == b) {
                spanning_ids.insert(a);
                horizontal_spanning_id.push_back(a);
            }
        }
    }


    value_type number_of_sites_in_spanning_cluster = number_of_site_in_spanning_clusters(spanning_ids);

    number_of_bonds_to_span.push_back(number_of_sites_in_spanning_cluster);

    if (debug_4_detectSpanning) {
        cout << "Number of sites to span " << number_of_sites_in_spanning_cluster << ','
             << "Number of active sites " << _total_number_of_active_bonds << ','
             << "Total number of sites " << maxSites() << endl;
        cout << "Number of cluters " << _clusters.size() << endl;
        cout << "Spanning Info : line " << __LINE__ << endl;
        cout << "Vertical spanning for cluster_ID   : " << vertical_spanning_id << endl;
        cout << "Horizontal spanning for cluster_ID : " << horizontal_spanning_id << endl;
        viewClusterExtended();
        viewLatticeExtended();
    }

    return !spanning_ids.empty();
}



/***********************************
 * Wrapping Detection
 **********************************/

/**
 *
 * @param bond  : last placed bond
 * @return
 */
/**
 * Wrapping is detected here. using the last placed bond.
 * @return bool. True if wrapping occured.
 */
bool BondPercolation_pb_v1::detectWrapping() {
    BondIndex bond = lastPlacedBond();
    // only possible if the cluster containing 'site' has bonds >= length of the lattice
    if(_number_of_occupied_bonds < length()){
        return false;
    }

    // check if it is already a wrapping site
    int id = _lattice.getBond(bond).get_groupID();
    int tmp_id{};
    for (auto i: _wrapping_indices){
        tmp_id = _lattice.getSite(i).get_groupID();
        if(id == tmp_id ){
//            cout << "Already a wrappig cluster : line " << __LINE__ << endl;
            return true;
        }
    }

    auto sites = SqLattice::get_neighbor_indices(length(), bond);
    IndexRelative r1 = _lattice.getSite(sites[0]).relativeIndex();
    IndexRelative r2 = _lattice.getSite(sites[1]).relativeIndex();

    if(abs(r1.x_ - r2.x_) > 1 || abs(r1.y_ - r2.y_) > 1){
//                cout << "Wrapping : line " << __LINE__ << endl;
        _wrapping_indices.push_back(sites[0]);
        return true;
    }


//    cout << "wrapping sites " << _wrapping_indices << endl;
    // if %_wrapping_indices is not empty but wrapping is not detected for the current site (%site)
    // that means there is wrapping but not for the %site
    return !_wrapping_indices.empty();
}

std::vector<double> BondPercolation_pb_v1::spanningProbability() const {
    vector<double> x(number_of_bonds_to_span.size());
    for (value_type i{}; i != x.size(); ++i) {
        x[i] = number_of_bonds_to_span[i] / double(maxSites());
    }
    return x;
}


double BondPercolation_pb_v1::occupationProbability() const {
    return _number_of_occupied_bonds / double(maxBonds());
}


/**
 *
 */
void BondPercolation_pb_v1::periodicity_status() {
    cout << "Periodicity = (" << _periodicity << ")" << endl;
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
void BondPercolation_pb_v1::subtract_entropy_for_site(const set<value_type>& found_index){
    double nos, mu_site;
    double H{};
    for(auto x : found_index){
        nos = _clusters[x].numberOfSites();
        mu_site = nos / maxSites();
        H += log(mu_site) * mu_site;
    }
    H *= -1;
    _entropy -= H;
}

/**
 * Must be called after merging the clusters
 * Cluster length is measured by sites
 * @param index
 */
void BondPercolation_pb_v1::add_entropy_for_site(value_type index){
    double nos = _clusters[index].numberOfSites();
    double mu = nos / maxSites();
    double H = -log(mu) * mu;
    _entropy += H;   // keeps track of entropy
}

/**
 * Very efficient approach for calculating entropy
 * @return
 */
double BondPercolation_pb_v1::entropy() {

    double mu = 1.0/ maxSites();
    double number_of_cluster_with_size_one = (maxSites() - sites_in_cluster_with_size_greater_than_one);
    double H = number_of_cluster_with_size_one * mu * log(mu);
    H *= -1;
    _entropy_current = _entropy + H;
    return _entropy_current;
}





