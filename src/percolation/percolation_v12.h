//
// Created by shahnoor on 11/22/19.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_V12_H
#define SQLATTICEPERCOLATION_PERCOLATION_V12_H

#include <random>
#include "../lattice/lattice.h"

/*****
 * Percolation class
 * version 12
 * Bond and Site percolation class will be extended from here
 */
class Percolation_v12{
    int  _length;
    size_t _max_number_of_bonds;
    size_t _max_number_of_sites;
    char type{'0'}; // percolation type. 's' -> site percolation. 'b' -> bond percolation
protected:

    // structural variables of lattice
    Lattice_v12 _lattice;

    value_type _index_sequence_position{};
    // cluster

    std::vector<Cluster_v12> _clusters;   // check and remove reapeted index manually
    // every birthTime we create a cluster we assign an set_ID for them
    int _cluster_id{};

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
    size_t _random_state;

public:
    ~Percolation_v12() = default;
    Percolation_v12() = default;

    explicit Percolation_v12(int length);

    void setRandomState(size_t seed, bool generate_seed=true);
    size_t getRandomState();

    size_t maxSites() const {return _max_number_of_sites;}
    size_t maxBonds() const {return _max_number_of_bonds;}

};

class SqLatticeRegularSite :  public Percolation_v12{

public:
    SqLatticeRegularSite() = default;
    explicit SqLatticeRegularSite(int length):Percolation_v12(length){}

    void init();

};

class SqLatticeRSBDSite :  public SqLatticeRegularSite{

};

class SqLatticeExplosiveSite :  public SqLatticeRegularSite{

};

class SqLatticeRegularBond :  public Percolation_v12{

};

class SqLatticeExplosiveBond :  public SqLatticeRegularBond{

};


#endif //SQLATTICEPERCOLATION_PERCOLATION_V12_H
