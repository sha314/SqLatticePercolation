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
    _max_bonds = 2*_length_squared;
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
        cout << "cluster [" << i << "] : ID (" << _clusters[i].ID() << "){" << endl;
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
 */
void SqLatticePercolation::viewSiteByRelativeIndex(){
    _lattice.view_sites_by_relative_index();
}

