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
    long double _entropy{};
    long double _entropy_current{};
    long double _entropy_previous{};
    double _largest_jump_entropy{}; // lrgest jump in entropy
    double _entropy_jump_pc{}; // at what pc there have been the largest jump
    size_t _cluster_count{};


    value_type _total_relabeling{};
    double time_relabel{};
    value_type _number_of_occupied_sites{};
    value_type _number_of_occupied_bonds{};
    value_type _max_iteration_limit{};
    std::vector<int> index_sequence;
    std::vector<int> randomized_index;
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

    void viewCluster();
    void viewLattice_original(){_lattice.view_as_assigned();};
    void viewLattice(){_lattice.view_all();};
    void viewLattice_by_id(){_lattice.view_by_id();};
    void viewLattice_by_gid(){_lattice.view_by_gid();};
    void viewLattice_by_relative_index(){_lattice.view_by_relative_index();};
    int length(){ return  _lattice.length();}

};

class SqLatticeRegularSite :  public Percolation_v12{
    size_t index_counter{};
    int id_last_site{};
    std::vector<int> _wrapping_site_ids;
public:
    SqLatticeRegularSite() = default;
    explicit SqLatticeRegularSite(int length);
    void reset();

    void init();

    bool occupy();

    IndexRelative  relabel_new_site(int id_current);
    void relabel(Cluster_v12& clstr, int id_current);
    void relabel_v2(Cluster_v12& clstr, int id_current, IndexRelative dx_dy);
    void relabel_v3(int id_current_a, std::vector<Index>& neighbors_a, Cluster_v12 &clstr_b);

    IndexRelative getRelativeIndexDX(Index root, Index site_new);
    IndexRelative getRelativeIndexDX_v2(Index root, Index site_new);

        int lastSite() const {return id_last_site;}
    Index lastSiteIndex()  {return _lattice.getSite(id_last_site).get_index();}

    int sign(int a);

    bool detectWrapping();
    int wrappingSite_id();
    Index wrappingSite();

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
