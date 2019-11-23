//
// Created by shahnoor on 11/22/19.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_V12_H
#define SQLATTICEPERCOLATION_PERCOLATION_V12_H

#include <random>
#include "../lattice/lattice.h"

class Percolation{
    int  _length;
    int _max_number_of_bonds;
    int _max_number_of_sites;
    char type{'0'}; // percolation type. 's' -> site percolation. 'b' -> bond percolation
protected:

    // structural variables of lattice
    SqLattice _lattice;

    value_type _index_sequence_position{};
    // cluster

    std::vector<Cluster_v3> _clusters;   // check and remove reapeted index manually
    // every birthTime we create a cluster we assign an set_ID for them
    int _cluster_id{};
    size_t number_of_clusters{};
    value_type min_index{};
    value_type max_index{};

    double _occuption_probability {};
    // entropy
    double _entropy{};
    double _entropy_current{};
    double _entropy_previous{};
    double _largest_jump_entropy{}; // lrgest jump in entropy
    double _entropy_jump_pc{}; // at what pc there have been the largest jump
    size_t _cluster_count{};
    value_type _bonds_in_cluster_with_size_two_or_more{0};   // total number of bonds in the clusters. all cluster has bonds > 1


    value_type _total_relabeling{};
    double time_relabel{};
    value_type _number_of_occupied_sites{};
    value_type _number_of_occupied_bonds{};
    value_type _max_iteration_limit{};
    std::vector<value_type> randomized_index;
    // order parameter calculation ingradiants
    // id of the cluster which has maximum number of bonds. used to calculate order parameter
    value_type _number_of_bonds_in_the_largest_cluster{};
    value_type _number_of_sites_in_the_largest_cluster{};   // might be useful later
    value_type _index_last_modified_cluster{};  // id of the last modified cluster
    std::mt19937 _random;
    value_type _random_state;
};

class SqLatticeRegularSite :  public Percolation{

};

class SqLatticeRSBDSite :  public SqLatticeRegularSite{

};

class SqLatticeExplosiveSite :  public SqLatticeRegularSite{

};

class SqLatticeRegularBond :  public Percolation{

};

class SqLatticeExplosiveBond :  public SqLatticeRegularBond{

};


#endif //SQLATTICEPERCOLATION_PERCOLATION_V12_H
