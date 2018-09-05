//
// Created by shahnoor on 10/5/2017.
//

#include "percolation.h"

/**
 *
 * @param length
 */
SqLatticePercolation::SqLatticePercolation(value_type length) {
    if (length <= 2) {
        /*
         * Because if _length=2
         * there are total of 4 distinct bond. But it should have been 8, i.e, (2 * _length * _length = 8)
         */
        cerr << "_length <= 2 does not satisfy _lattice properties for percolation : line" << __LINE__ << endl;
        exit(1);
    }
    _length = length;
    value_type _length_squared = length * length;
    _max_number_of_bonds = 2*_length_squared;
    _max_number_of_sites = _length_squared;
    _clusters = vector<Cluster_v2>();
    min_index = 0;
    max_index = _length - 1;

//    size_t seed = 0;
//    cerr << "automatic seeding is commented : line " << __LINE__ << endl;

    auto seed = _random_device();

    _g.seed(seed); // seeding
    cout << "seeding with " << seed << endl;
}


/**
 *
 */
void SqLatticePercolation::viewCluster() {
    cout << "clusters with numberOfBonds greater than 1" << endl;
    value_type total_bonds{}, total_sites{};

    for (value_type i{}; i != _clusters.size(); ++i) {
        if(_clusters[i].empty()){
//            cout << "Empty cluster : line " << endl;
            continue;
        }
        cout << "cluster [" << i << "] : " << '{' << endl;
        cout << _clusters[i];
        total_bonds += _clusters[i].numberOfBonds();
        total_sites += _clusters[i].numberOfSites();

        cout << '}' << endl;
    }
    cout << "Total bonds " << total_bonds << endl;
    cout << "Total sites " << total_sites << endl;
}



/**
 * Extended version of view_cluster
 */
void SqLatticePercolation::viewClusterExtended() {
    cout << "clusters with numberOfBonds greater than 1" << endl;
    value_type total_bonds{}, total_sites{};

    std::vector<Index> sites;
    std::vector<BondIndex> bonds;
    for (value_type i{}; i != _clusters.size(); ++i) {
        if(_clusters[i].empty()){
//            cout << "Empty cluster : line " << endl;
            continue;
        }
        cout << "cluster [" << i << "] : ID (" << _clusters[i].get_ID() << "){" << endl;
        // printing sites
        sites = _clusters[i].getSiteIndices();
        cout << "Sites : size (" << sites.size() << ") : ";
        cout << '{';
        for (auto a: sites) {
            cout << a << ',';
        }
        cout << '}' << endl;

        bonds = _clusters[i].getBondIndices();
        cout << "Bonds : size (" << bonds.size() << ") : ";
        cout << '{';
        for (auto a: bonds) {
            if (a.horizontal()) {
                // horizontal bond
                cout << _lattice.get_h_bond({a.row_, a.column_}) << ',';
            } else if (a.vertical()) {
                // vertical bond
                cout << _lattice.get_v_bond({a.row_, a.column_}) << ',';
            } else {
                cout << '!' << a << '!' << ','; // bond is not valid
            }
        }
        cout << '}';

        cout << endl;

        total_bonds += _clusters[i].numberOfBonds();
        total_sites += _clusters[i].numberOfSites();

        cout << '}' << endl;
    }
    cout << "Total bonds " << total_bonds << endl;
    cout << "Total sites " << total_sites << endl;
}

/**
 *
 * @param site
 * @param bond
 * @param total_site
 * @param total_bond
 */
void
SqLatticePercolation::get_cluster_info(
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
//    cout << "tatal sites " << total_site << endl;
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
//    cout << "tatal sites " << total_site << endl;

}

void SqLatticePercolation::reset() {
    _lattice.reset();
    _cluster_index_from_id.clear();
    _clusters.clear();
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
}

void SqLatticePercolation::jump() {

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



//
//const vector<double> &SqLatticePercolation::get_occupation_probabilities() const {
//    return _occupation_probabilities;
//}
//
//const vector<double> &SqLatticePercolation::get_nob_spanning() const {
//    return _nob_spanning;
//}
//
//const vector<double> &SqLatticePercolation::get_nob_largest() const {
//    return _nob_largest;
//}
//
//const vector<double> &SqLatticePercolation::get_nos_spanning() const {
//    return _nos_spanning;
//}
//
//const vector<double> &SqLatticePercolation::get_nos_largest() const {
//    return _nos_largest;
//}
//
//const vector<double> &SqLatticePercolation::get_entropy_sites() const {
//    return _entropy_sites;
//}
//
//const vector<double> &SqLatticePercolation::get_entropy_bonds() const {
//    return _entropy_bonds;
//}
//
//
//const vector<double> &SqLatticePercolation::get_pcs() const {
//    return _pcs;
//}
//
//const vector<double> &SqLatticePercolation::get_spanning_cluster_size_sites() const {
//    return _spanning_cluster_size_sites;
//}
//
//const vector<double> &SqLatticePercolation::get_spanning_cluster_size_bonds() const {
//    return _spanning_cluster_size_bonds;
//}

