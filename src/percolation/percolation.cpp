//
// Created by shahnoor on 10/5/2017.
//

#include "percolation.h"

/**
 *
 * @param length
 */
SqLatticePercolation::SqLatticePercolation(value_type length) {
    _length = length;
    _length_squared = length * length;
    _max_number_of_bonds = 2*_length_squared;
    _max_number_of_sites = _length_squared;
    _clusters = vector<Cluster_v2>();

    min_index = 0;
    max_index = _length - 1;
}


/**
 *
 */
void SqLatticePercolation::viewCluster() {
    cout << "clusters with numberOfBonds greater than 1" << endl;
    value_type total_bonds{}, total_sites{};

    for (value_type i{}; i != _clusters.size(); ++i) {
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

