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

    value_type min_index{};
    value_type max_index{};

    double _occuption_probability {};
    // entropy
    /*
     * _entropy_list contains value of entropy for different clusters.
     * id of a cluster is the index of this array that will give the entropy for that cluster
     */
    std::vector<long double> _entropy_list;
    long double _entropy_{}; // for jump and specific heat caclculation
    long double _entropy_site{},_entropy_bond{};
    long double _entropy_current{};
    long double _entropy_previous{};
    double _largest_jump_entropy{}; // lrgest jump in entropy
    double _entropy_jump_pc{}; // at what pc there have been the largest jump
    size_t _cluster_count{};


    value_type _total_relabeling{};
    double time_relabel{};
    value_type _number_of_occupied_sites{};
    value_type _number_of_occupied_bonds{};
    size_t _max_iteration_limit{};
    std::vector<int> index_sequence;
    std::vector<int> randomized_index;
    // order parameter calculation ingradiants
    // id of the cluster which has maximum number of bonds. used to calculate order parameter
    value_type _number_of_bonds_in_the_largest_cluster{};
    value_type _number_of_sites_in_the_largest_cluster{};   // might be useful later
    int bond_largest_cluster_index{};
    int site_largest_cluster_index{};
    value_type _index_last_modified_cluster{};  // id of the last modified cluster
    std::mt19937 _random;
    size_t _random_state;

    // self cluster jump flag
    bool _self_cluster_jump_bond{false};
    bool _self_cluster_jump_site{false};
public:
    ~Percolation_v12() = default;
    Percolation_v12() = default;

    explicit Percolation_v12(int length);

    virtual void reset();

    void setRandomState(size_t seed, bool generate_seed=true);
    size_t getRandomState();

    size_t maxSites() const {return _max_number_of_sites;}
    size_t maxBonds() const {return _max_number_of_bonds;}

    void viewCluster();
    void viewLattice_original(){_lattice.view_as_assigned();};
    void viewLattice(){ _lattice.view_all_v2();};
    void viewLattice_by_id(){ _lattice.view_by_id_v2();};
    void viewLattice_by_gid(){ _lattice.view_by_gid_v2();};
    void viewLattice_by_relative_index(){ _lattice.view_by_relative_index_v2();};
    void viewLattice_sites_by_gid(){ _lattice.view_sites_by_gid();}
    int length(){ return  _lattice.length();}

    size_t numberOfBondsInTheLargestCluster();
    size_t numberOfSitesInTheLargestCluster();

    // must define methods
    virtual size_t maxIterationLimit() = 0;
    // numerical properties
    virtual long double entropy() = 0;
    virtual double occupationProbability() = 0;
    virtual size_t wrappingClusterSize() = 0;

    virtual size_t numberOfSitesInTheWrappingClusters() = 0;
    virtual size_t numberOfBondsInTheWrappingClusters() = 0;

    void track_clusters(int root);
};

//
//
//class SqLatticeRSBDSite :  public SitePercolation_ps_v12{
//
//};
//
//class SqLatticeExplosiveSite :  public SitePercolation_ps_v12{
//
//};
//
//class SqLatticeRegularBond :  public Percolation_v12{
//
//};
//
//class SqLatticeExplosiveBond :  public SqLatticeRegularBond{
//
//};


#endif //SQLATTICEPERCOLATION_PERCOLATION_V12_H
