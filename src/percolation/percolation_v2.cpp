//
// Created by shahnoor on 6/28/19.
//

#include "percolation_v2.h"

using namespace std;


/**
 *
 * @param length
 */
SqLatticePercolation_v2::SqLatticePercolation_v2(value_type length) {
    if (length <= 2) {
        /*
         * Because if _length=2
         * there are total of 4 distinct bond. But it should have been 8, i.e, (2 * _length * _length = 8)
         */
        cerr << "_length <= 2 does not satisfy lattice_ref properties for percolation : line" << __LINE__ << endl;
        exit(1);
    }
    _length = length;
    value_type _length_squared = length * length;
    _max_number_of_bonds = 2*_length_squared;
    _max_number_of_sites = _length_squared;
    _clusters = vector<Cluster_v4>();
    min_index = 0;
    max_index = _length - 1;
}


/**
 *
 */
void SqLatticePercolation_v2::viewCluster() {
    cout << "clusters with numberOfBonds greater than 1" << endl;
    value_type total_bonds{}, total_sites{};

    for (value_type i{}; i != _clusters.size(); ++i) {
        if(_clusters[i].empty()){
//            cout << "Empty cluster : line " << endl;
            continue;
        }
        cout << "cluster [" << i << "] : " << '{' << endl;
        //////>>>>>>>>
        cout << "Sites : size (" << _clusters[i].numberOfSites() << ") : ";
        cout << '{';
        for(auto a: _clusters[i].getSiteIndices()){
            cout << _lattice.translate1DToSite(a) << ',';
        }
        cout << '}' << endl;

        cout << "Bonds : size (" << _clusters[i].numberOfBonds() <<") : ";
        cout << '{';
        for(auto a: _clusters[i].getBondIndices()){
            cout << _lattice.translate1DToBond(a) << ',';
        }
        cout << '}';
        /////<<<<<<<
        cout << '}' << endl;
        total_bonds += _clusters[i].numberOfBonds();
        total_sites += _clusters[i].numberOfSites();
    }
    cout << "Total bonds " << total_bonds << endl;
    cout << "Total site_index_sequence " << total_sites << endl;
}



/**
 * Extended version of view_cluster
 */
void SqLatticePercolation_v2::viewClusterExtended() {
//    cout << "clusters with numberOfBonds greater than 1" << endl;
    value_type total_bonds{}, total_sites{};

    std::vector<uint> sites;
    std::vector<uint> bonds;
    size_t cluster_count{};
    for (value_type i{}; i != _clusters.size(); ++i) {
        if(_clusters[i].empty()){
//            cout << "Empty cluster : line " << endl;
            continue;
        }
        cout << "cluster [" << i << "] : ID (" << _clusters[i].get_ID() << "){" << endl;
        // printing site_index_sequence
        sites = _clusters[i].getSiteIndices();
        cout << "Sites : size (" << sites.size() << ") : ";
        cout << '{';
        for (auto a: sites) {
            cout << _lattice.translate1DToSite(a) << ',';
        }
        cout << '}' << endl;

        bonds = _clusters[i].getBondIndices();
        cout << "Bonds : size (" << bonds.size() << ") : ";
        cout << '{';
        for (auto a: bonds) {
            cout << _lattice.translate1DToBond(a) << ',';
        }
        cout << '}';

        cout << endl;

        total_bonds += _clusters[i].numberOfBonds();
        total_sites += _clusters[i].numberOfSites();

        cout << '}' << endl;
        ++cluster_count;
    }
    cout << "Total bonds " << total_bonds << endl;
    cout << "Total sites " << total_sites << endl;
    cout << "Total clusters " << cluster_count << endl;
}

/**
 *
 * @param site
 * @param bond
 * @param total_site
 * @param total_bond
 */
void
SqLatticePercolation_v2::get_cluster_info(
        vector<value_type> &site,
        vector<value_type> &bond
) {
    value_type total_site{}, total_bond{};
    site.clear();
    bond.clear();

    unsigned long size = _clusters.size();
    site.reserve(size);
    bond.reserve(size);

    value_type s, b;

    for(value_type i{}; i < size; ++i){
        if(_clusters[i].empty()){
//            cout << "Empty cluster : line " << endl;
            continue;
        }
        s = _clusters[i].numberOfSites();
        b = _clusters[i].numberOfBonds();
        site.push_back(s);
        bond.push_back(b);
        total_site += s;
        total_bond += b;
    }
    if(site.size() != bond.size()){
        cout << "Size mismatched : line " << __LINE__ << endl;
    }
//    cout << "total bonds " << total_bond << endl;
//    cout << "tatal site_index_sequence " << total_site << endl;
    if(type == 's'){
        for(value_type j{total_bond}; j < maxBonds(); ++j){
            bond.push_back(1); // cluster of length 1
            total_bond += 1;
        }
    }
    if(type == 'b'){
        for(value_type j{total_site}; j < maxSites(); ++j){
            total_site += 1;
            site.push_back(1); // cluster of length 1
        }
    }
//    if(site.size() != bond.size()){
//        cout << "Size mismatched : line " << __LINE__ << endl;
//    }
//    cout << "total bonds " << total_bond << endl;
//    cout << "tatal site_index_sequence " << total_site << endl;

}

void SqLatticePercolation_v2::reset() {
    _lattice.reset();
//    cout << _clusters.size() << " : line " << __LINE__ << endl;
    _clusters.clear(); // SIGABRT for SitePercolation_ps_v9 in cluster size distribution
//    cout << " line " << __LINE__ << endl;
    _index_sequence_position = 0;
    _cluster_id = 0;
    _occuption_probability = 0;
    // entropy
    _entropy=0;
    _entropy_current=0;
    _entropy_previous=0;
    _largest_jump_entropy = 0;
    _entropy_jump_pc = 0;
    _total_relabeling = 0;
    time_relabel = 0;
    _cluster_count = 0;
    _number_of_bonds_in_the_largest_cluster=0;
    _number_of_sites_in_the_largest_cluster=0;
    _index_last_modified_cluster=0;
}

void SqLatticePercolation_v2::jump() {

    double delta_H{};
    if(_index_sequence_position > 1) {
        delta_H = _entropy_current - _entropy_previous;
    }
    if(abs(delta_H) > abs(_largest_jump_entropy)){
        _largest_jump_entropy = delta_H;
        _entropy_jump_pc = _occuption_probability;
    }
    _entropy_previous = _entropy_current; // be ready for next step
}


/**
 * Entropy calculation is performed here. The fastest method possible.
 * Cluster size is measured by site.
 * @return current entropy of the lattice
 */
double SqLatticePercolation_v2::entropy_by_site() {
    double H{}, mu ;

    for(size_t i{}; i < _clusters.size(); ++i){
        if(!_clusters[i].empty()){
            mu = _clusters[i].numberOfSites() / double(_number_of_occupied_sites);
            H += mu*log(mu);
        }
    }

    return -H;
}

/**
 * Entropy calculation is performed here. The fastest method possible.
 * Cluster size is measured by site.
 * @return current entropy of the lattice
 */
double SqLatticePercolation_v2::entropy_by_bond() {
    double H{}, mu ;

    for(size_t i{}; i < _clusters.size(); ++i){
        if(!_clusters[i].empty()){
            mu = _clusters[i].numberOfBonds() / double(maxBonds());
            H += mu*log(mu);
        }
    }

    double number_of_cluster_with_size_one = maxBonds() - _bonds_in_cluster_with_size_two_or_more;
//    cout << " _bonds_in_cluster_with_size_two_or_more " << _bonds_in_cluster_with_size_two_or_more << " : line " << __LINE__ << endl;
    mu = 1.0/double(maxBonds());
    H += number_of_cluster_with_size_one * log(mu) * mu;

    return -H;
}

void SqLatticePercolation_v2::setRandomState(size_t seed, bool generate_seed) {
//    size_t seed = 0;
//    cerr << "automatic seeding is commented : line " << __LINE__ << endl;
    _random_state = seed;
    if(generate_seed) {
        std::random_device _rd;
        _random_state = _rd();
    }else{
        cerr << "generate_seed = false : line " << __LINE__ << endl;
    }
    _random.seed(_random_state); // seeding
    cout << "seeding with " << _random_state << endl;
}


value_type SqLatticePercolation_v2::getRandomState() {
//    size_t seed = 0;
//    cerr << "automatic seeding is commented : line " << __LINE__ << endl;
//    std::random_device _rd;
//    auto seed = _rd();
//    _random.seed(seed); // seeding
//    cout << "seeding with " << seed << endl;
    return _random_state;
}

void SqLatticePercolation_v2::ckeckCluster() {
//    cout << "clusters with numberOfBonds greater than 1" << endl;
    value_type total_bonds{}, total_sites{};

    std::vector<uint> sites;
    std::vector<uint> bonds;
    value_type count_bonds{},    count_sites{};

    for (value_type i{}; i != _clusters.size(); ++i) {
        if(_clusters[i].empty()){
//            cout << "Empty cluster : line " << endl;
            continue;
        }


        sites = _clusters[i].getSiteIndices();
        bonds = _clusters[i].getBondIndices();


        count_bonds = _clusters[i].numberOfBonds();
        count_sites = _clusters[i].numberOfSites();

        if(count_bonds != bonds.size()){
            cout << "count_bonds != bonds.size()" << endl;
        }

        if(count_sites != sites.size()){
            cout << "count_sites != sites.size()" << endl;
        }

        total_bonds += _clusters[i].numberOfBonds();
        total_sites += _clusters[i].numberOfSites();


    }
    cout << "Total bonds " << total_bonds << endl;
    cout << "Total sites " << total_sites << endl;
}

void SqLatticePercolation_v2::init() {
//    setRandomState(0, true);
}

IndexRelative SqLatticePercolation_v2::getRelativeIndex(Index root, Index site_new){
//    cout << "Entry \"SitePercolation_ps_v10::getRelativeIndex\" : line " << __LINE__ << endl;
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
 * Take a bond index only if the corresponding site is active
 * takes longer? time than version 1?, i.e.,  connection()
 * @param site
 * @param site_neighbor
 * @param bond_neighbor
 */
void SqLatticePercolation_v2::connection_periodic(Index site, vector<Index> &site_neighbor, vector<BondIndex> &bond_neighbor)
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

//    if(!lattice_ref.getSite(site_neighbor[0]).isActive()) {
//        bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
//    }
//    if(!lattice_ref.getSite(site_neighbor[1]).isActive()){
//        bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
//    }
//    if(!lattice_ref.getSite(site_neighbor[2]).isActive()){
//        bond_neighbor.push_back({BondType::Vertical,    site.row_, site.column_});
//    }
//    if(!lattice_ref.getSite(site_neighbor[3]).isActive()) {
//        bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
//    }
    bond_neighbor.resize(4);
    bond_neighbor[0] = {BondType::Horizontal, site.row_, site.column_};
    bond_neighbor[1] = {BondType::Horizontal, site.row_, prev_column};
    bond_neighbor[2] = {BondType::Vertical,    site.row_, site.column_};
    bond_neighbor[3] = {BondType::Vertical, prev_row, site.column_};

}



