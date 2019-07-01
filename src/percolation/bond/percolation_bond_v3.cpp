//
// Created by shahnoor on 10/31/2017.
//

#include <cstdlib>
#include <climits>
#include <unordered_set>
#include <algorithm>
#include "../percolation.h"
#include "../../util/printer.h"
#include "percolation_bond_v3.h"


using namespace std;



/**
 *
 * @param length
 */
BondPercolation_pb_v3::BondPercolation_pb_v3(value_type length, bool periodicity)
        : SqLatticePercolation_v2(length) {
    std::cout << "Constructing BondPercolation_pb_v3 object : line " << __LINE__ << endl;
    SqLatticePercolation_v2::set_type('b');
    _periodicity = periodicity;
    _index_sequence_position = 0;
    _number_of_occupied_bonds = 0;
    _max_iteration_limit = maxBonds();
    randomized_index.resize(maxIterationLimit());
    for(value_type i{}; i < maxIterationLimit(); ++i){randomized_index[i] = i;}
    _clusters = vector<Cluster_v4>();

//    init();
}


/**
 * some methods that are directly affected by the seed of random number generator
 * + methods that requires a bit computaion (so not in constructor)
 *
 */
void BondPercolation_pb_v3::init(bool random_seed) {
//    SqLatticePercolation::init();
    _lattice = SqLattice_v2(length());
    _lattice.init(false, true, true, true);
    if(random_seed) setRandomState(0, random_seed);
    initialize_index_sequence(); // not need to put in reset method
    initialize();
    initialize_cluster();
    randomize_v2();  // randomize the untouched_site_indices
}

/**
 * Initialize the Class Object
 */
void BondPercolation_pb_v3::initialize() {

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
void BondPercolation_pb_v3::initialize_index_sequence() {
    bond_index_sequence = _lattice.getBonds();
    site_index_sequence = _lattice.getSites();

    if (!_periodicity) {
        bond_index_sequence.clear();
        bond_index_sequence.reserve(maxBonds());
        value_type row{}, col{};
        value_type limit = maxSites() - 1; // since if row==4 and col ==4 we don't want that bond
        for (value_type i{}; i < limit; ++i) {

            if (row == max_index) {
                bond_index_sequence.push_back({BondType::Horizontal, row, col});
            } else if (col == max_index) {
                bond_index_sequence.push_back({BondType::Vertical, row, col});
            } else {
                bond_index_sequence.push_back({BondType::Horizontal, row, col});
                bond_index_sequence.push_back({BondType::Vertical, row, col});
            }

            ++col;
            if (col == length()) {
                col = 0;
                ++row;
            }
        }
    }
}



/**
 * Randomize the site placing order
 * Takes 3.031000 sec for 1000 x 1000 sites
 */
void BondPercolation_pb_v3::randomize_v2() {
//    cout << "default_seed " << _random.default_seed << endl;
    std::shuffle(randomized_index.begin(), randomized_index.end(), _random);
//    cout << "Index sequence : " << randomized_index_sequence << endl;
}


/**
 * Reset all calculated values and then call initiate()
 * to initiallize for reuse
 * must call reset before actually running the program for data generation.
 * Because random state is not assigned by the constructor rather by the user
 */
void BondPercolation_pb_v3::reset() {
    SqLatticePercolation_v2::reset();
    _total_number_of_active_bonds = 0;
    _index_sequence_position = 0;
    _number_of_occupied_bonds = 0;
    number_of_bonds_to_span.clear();

    sites_in_cluster_with_size_greater_than_one = 0;
    _wrapping_indices.clear();

    initialize();
    initialize_cluster();

//    cout << "use suffle : line " << __LINE__ << endl;
    randomize_v2();
}



/**
 *
 */
double BondPercolation_pb_v3::entropy_slow() {
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
void BondPercolation_pb_v3::relabel_sites(const Cluster_v4 &clstr, int id) {

    const vector<uint> sites = clstr.getSiteIndices();
    for (auto a: sites) {
        _lattice.setGroupIDSite(a, id);
    }

}

/**
 * Before this method is called, current bond must have a suitable group id.
 *
 * @param bond     : last placed bond
 * @param clstr_b  : given cluster
 * @param bond_pos : position of the bonds from where bonds will be relabeled
 *                   according to bond_pos-1 bond ?
 * @param site_pos : position of the sites from where sites will be relabeled
 *                   according to site_pos-1 site ?
 */
void BondPercolation_pb_v3::relabel_cluster(BondIndex bond, const Cluster_v4& clstr_b, size_t bond_pos, size_t site_pos){
    const vector<uint> & sites = clstr_b.getSiteIndices();
    const vector<uint> & bonds = clstr_b.getBondIndices();

    int id = _lattice.getBond(bond).get_groupID();
    auto neighbor_sites = _lattice.get_neighbor_sites(bond);
    Index site_a = neighbor_sites[0];
    Index site_b = neighbor_sites[1];

//    int id_a = _lattice.getSite(site_a).get_groupID();
//    int id_b = _lattice.getSite(site_b).get_groupID();

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
    cout << "site_index_sequence to be relabeled :" ;
    copy(site_index_sequence.begin()+site_pos, site_index_sequence.end(), std::ostream_iterator<Index>(cout));
    cout << endl;
    cout << "site_b " << site_b << endl;
    // END debugging*/

    // relabeling bonds
    for(size_t i{bond_pos}; i < bonds.size(); ++i){
        _lattice.setGroupIDBond(bonds[i], id);
    }

    // relabeling site_index_sequence
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
    uint a;
    IndexRelative relative_site__a;
    for (value_type i{site_pos}; i < sites.size(); ++i) {
        a = sites[i];
        _lattice.setGroupIDSite(sites[i], id);
        relative_site__a = _lattice.getRelativeIndex(a);
//        cout << "index " << a << " relative index " << relative_site__a << " to ";
        x = relative_site__a.x_ + delta_x_ab;
        y = relative_site__a.y_ + delta_y_ab;
        _lattice.setRelativeIndex(a, {x, y});
//        cout << _lattice.getSite(a).relativeIndex() << endl;
    }
}

void BondPercolation_pb_v3::relabel_sites(const vector<uint> &sites, int id_a, int delta_x_ab, int delta_y_ab)  {
    int x, y;
    uint a;
    IndexRelative relative_site__a;
    for (value_type i = 0; i < sites.size(); ++i) {
        a = sites[i];
        _lattice.setGroupIDSite(a, id_a);
        relative_site__a = _lattice.getRelativeIndex(a);
        x = relative_site__a.x_ + delta_x_ab;
        y = relative_site__a.y_ + delta_y_ab;
        _lattice.setRelativeIndex(a, {x, y});
    }
}


/**
 * All site placing method in one place
 *
 * @return true if operation is successfull
 */
bool BondPercolation_pb_v3::occupy() {
    if (_index_sequence_position >= maxBonds()) {
        return false;
    }

    value_type v = placeBond_v1();
    _occuption_probability = occupationProbability(); // for super class
//    auto s = SqLattice::get_neighbor_indices(length(), _last_placed_bond);
//    if(_lattice.getSite(s[0]).get_groupID() != _lattice.getSite(s[1]).get_groupID()){
//        cout << "relabeling is not perfect for site " << endl;
//        for(auto a : s) {
//            cout << a << " id " << _lattice.getSite(a).get_groupID() << endl;
//        }
//        return false; // to reminate iteration
//    }
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
value_type BondPercolation_pb_v3::placeBond_v1() {

    if (_number_of_occupied_bonds >= maxBonds()){
        return ULONG_MAX;// unsigned long int maximum value
    }
    value_type index = randomized_index[_index_sequence_position];
    ++_index_sequence_position;
    BondIndex current_bond = bond_index_sequence[index];
    _last_placed_bond = current_bond;
//    cout << "placing bond " << current_bond << " : " << __LINE__ << endl;

    _lattice.activateBond(current_bond);

    ++_number_of_occupied_bonds;


    // merging clusters

    // one bond connects site
    auto sites = _lattice.get_neighbor_sites(current_bond);
//    cout << "neighbor " << sites << endl;
    if(sites.size() != 2){
        cerr << "How can a bond can link more than two sites " << __LINE__ << endl;
        return 0;
    }
    int id_1 = _lattice.getGroupID(sites[0]);
    int id_2 = _lattice.getGroupID(sites[1]);
//    cout << "id1 " << id_1 << ", id2 " << id_2 << endl;
    if(id_1 < 0 || id_2 < 0){
        cerr << "negative group id!!! " << __LINE__ << endl;
        return 0;
    }
    int base_index = 0;
    if(id_1 == id_2){
        _clusters[id_1].addBondIndex(length(), current_bond);
        _lattice.setGroupID(current_bond, id_1);
//        cout << "belongs to the same cluster" << endl;
    }else {
        size_t size_1 = _clusters[id_1].numberOfSites();
        size_t size_2 = _clusters[id_2].numberOfSites();

        if(size_1 > size_2){
            base_index = id_1;
        }else{
            base_index = id_2;
            id_2 = id_1;
            // caution : both id_1 and id_2 are equal.
        }
        _clusters[base_index].addBondIndex(length(),current_bond);
        _lattice.setGroupID(current_bond, base_index);
        _clusters[base_index].insert_v2(_clusters[id_2]);
        relabel_cluster(current_bond, _clusters[id_2], 0, 0);
        _clusters[id_2].clear(); // clear the cluster
    }



    // end merging
//    subtract_entropy_for_site(found_index_set, base_id);  // tracking entropy change
//    value_type base_index = manageClusters(
//            found_index_set, sites, current_bond, base_id
//    );
//    add_entropy_for_site(base_index); // tracking entropy change

    // running tracker
//    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster
    return base_index;
}



void BondPercolation_pb_v3::calculate_spanning_probability() {
    calculate_spanning_probability_by_largest_cluster();
}


/**
 * Use Group_ID to identify Bond and Site in the same cluster
 */
void BondPercolation_pb_v3::calculate_spanning_probability_by_largest_cluster() {
//    cout << "calculate_spanning_probability_by_largest_cluster() : line " << __LINE__ << endl;
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
value_type BondPercolation_pb_v3::count_number_of_active_site() {
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
//vector<Index> BondPercolation_pb_v3::get_Sites_for_bonds(BondIndex index) {
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



value_type BondPercolation_pb_v3::number_of_site_in_spanning_clusters(unordered_set<int> g_ids) {
    value_type nos{};   // number of site_index_sequence
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
bool BondPercolation_pb_v3::detectSpanning() {
//    cout << "Entry -> detectSpanning() : line " << __LINE__ << endl;

    // if any of the clusters does not have site_index_sequence > length() --> no spanning
    value_type x{};
    for (auto a: _clusters) {
        if (x < a.numberOfSites())
            x = a.numberOfSites();
    }
    if (x < length()) {
        cout << "not enough site_index_sequence in any cluster to span : line " << __LINE__ << endl;
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
bool BondPercolation_pb_v3::detectWrapping() {
    bool a;
    a = detect_wrapping_v1();
//    a = detect_wrapping_v2();
    return a;
}

std::vector<double> BondPercolation_pb_v3::spanningProbability() const {
    vector<double> x(number_of_bonds_to_span.size());
    for (value_type i{}; i != x.size(); ++i) {
        x[i] = number_of_bonds_to_span[i] / double(maxSites());
    }
    return x;
}


double BondPercolation_pb_v3::occupationProbability() const {
    return _number_of_occupied_bonds / double(maxBonds());
}


/**
 *
 */
void BondPercolation_pb_v3::periodicity_status() {
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
void BondPercolation_pb_v3::subtract_entropy_for_site(const set<value_type>& found_index, int base){
    double nos, mu_site;
    double H{};
    if(base >= 0){
        nos = _clusters[base].numberOfSites();
        mu_site = nos / maxSites();
        H += log(mu_site) * mu_site;
    }
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
void BondPercolation_pb_v3::add_entropy_for_site(value_type index){
    double nos = _clusters[index].numberOfSites();
    double mu = nos / maxSites();
    double H = -log(mu) * mu;
    _entropy += H;   // keeps track of entropy
}

/**
 * Very efficient approach for calculating entropy
 * @return
 */
double BondPercolation_pb_v3::entropy() {

    double mu = 1.0/ maxSites();
    double number_of_cluster_with_size_one = (maxSites() - sites_in_cluster_with_size_greater_than_one);
    double H = number_of_cluster_with_size_one * mu * log(mu);
    H *= -1;
    _entropy_current = _entropy + H;
    return _entropy_current;
}


value_type BondPercolation_pb_v3::numberOfBondsInTheWrappingClusters(){
    value_type nob{};
    int id{};
    for(auto i: _wrapping_indices){
        id = _lattice.getSite(i).get_groupID();
        if(id >= 0) {
            nob += _clusters[id].numberOfBonds();
        }
    }
    return nob;
}

value_type BondPercolation_pb_v3::numberOfSitesInTheWrappingClusters(){
    value_type nos{};
    int id{};
    for(auto i: _wrapping_indices){
        id = _lattice.getSite(i).get_groupID();
        if(id >= 0) {
            nos += _clusters[id].numberOfSites();
        }
    }
    return nos;
}

/**
 * Only applicable if the number of bonds in the largest cluster is calculated when occupying the lattice.
 * Significantly efficient than the previous version numberOfBondsInTheLargestCluster()
 * @return
 */
value_type BondPercolation_pb_v3::numberOfBondsInTheLargestCluster() {
    return _number_of_bonds_in_the_largest_cluster;
}

/**
 * Only applicable if the number of bonds in the largest cluster is calculated when occupying the lattice.
 * Significantly efficient than the previous version numberOfBondsInTheLargestCluster()
 * @return
 */
value_type BondPercolation_pb_v3::numberOfSitesInTheLargestCluster() {
    return _number_of_bonds_in_the_largest_cluster;
}

/**
 * Condition: must be called each time a site is placed
 */
void BondPercolation_pb_v3::track_numberOfBondsInLargestCluster() {

    // calculating number of bonds in the largest cluster // by cluster index
    // checking number of bonds
    if(_clusters[_index_last_modified_cluster].numberOfBonds() > _number_of_bonds_in_the_largest_cluster){
        _number_of_bonds_in_the_largest_cluster = _clusters[_index_last_modified_cluster].numberOfBonds();
    }

}

/**
 *
 */
void BondPercolation_pb_v3::track_numberOfSitesInLargestCluster(){

    // calculating number of bonds in the largest cluster // by cluster index
    // checking number of bonds
    if(_clusters[_index_last_modified_cluster].numberOfSites() > _number_of_sites_in_the_largest_cluster){
        _number_of_sites_in_the_largest_cluster = _clusters[_index_last_modified_cluster].numberOfSites();
    }
}

const std::vector<double> BondPercolation_pb_v3::clusterSizeDistribution() const {
    double area=0; // area under curve should be 1
    vector<double> cluster_counts;
    size_t n, sz;
//    cout << "number of max clusters " << _clusters.size() << endl;
    size_t mx_n{};
    size_t total_sites=0;
    for(size_t i{}; i < _clusters.size(); ++i){
        if(_clusters[i].empty())   continue;
//        cout << i << " th cluster ";
        n = _clusters[i].numberOfSites();
        total_sites += n;
        if(n > mx_n){
            mx_n = n;
        }
        sz = cluster_counts.size();
        if (sz <= n){
            cluster_counts.resize(n+1);
        }
//        cout << cluster_counts[n];
        ++cluster_counts[n]; // increase cluster count

//        cout << endl;
    }
    cluster_counts[1] += length()*length() - total_sites; // isolated site_index_sequence forms custer of unit size
//    if(total_sites != length()*length()){
//        cout << "total_sites != length()*length() " << total_sites << endl;
//    }
//    cout << "before returning " << mx_n << endl;
    for(size_t i{}; i < cluster_counts.size(); ++i){
        area += cluster_counts[i] * i;
    }
//    cout << "area under curve " << area /(length()*length()) << endl;
    return cluster_counts;
}

bool BondPercolation_pb_v3::detect_wrapping_v1() {
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

    auto sites = _lattice.get_neighbor_sites(bond);
    IndexRelative r1 = _lattice.getSite(sites[0]).relativeIndex();
    IndexRelative r2 = _lattice.getSite(sites[1]).relativeIndex();

    if(abs(r1.x_ - r2.x_) > 1 || abs(r1.y_ - r2.y_) > 1){
//                cout << "Wrapping : line " << __LINE__ << endl;
        _wrapping_indices.push_back(sites[0]);
        return true;
    }


//    cout << "wrapping site_index_sequence " << _wrapping_indices << endl;
    // if %_wrapping_indices is not empty but wrapping is not detected for the current site (%site)
    // that means there is wrapping but not for the %site
    return !_wrapping_indices.empty();
}

void BondPercolation_pb_v3::initialize_cluster() {
    // all site_index_sequence are a unit cluster initially

    _clusters.resize(maxSites()); // since reset might change sice
    for(int i{}; i < site_index_sequence.size(); ++i){
        _clusters[i].clear(); // first clear the cluster
        _clusters[i].addSiteIndex(length(), site_index_sequence[i]);
        _clusters[i].set_ID(i);
        _lattice.setGroupID(site_index_sequence[i], i);
    }
}


