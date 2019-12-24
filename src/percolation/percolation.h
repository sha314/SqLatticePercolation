//
// Created by shahnoor on 10/5/2017.
//

#ifndef SITEPERCOLATION_PERCOLATION_H
#define SITEPERCOLATION_PERCOLATION_H


#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <climits>
#include <fstream>
#include <random>
#include "../flags.h"
#include "../types.h"
#include "../lattice/lattice.h"
#include "../index/index.h"
#include "calculation.h"
#include "site_position.h"
#include "id_index.h"
#include "boundary.h"


/**
 * class : SqLatticePercolation
 * used classed : SqLattice, Cluster_v3
 *
 * The Square Lattice Percolation class
 */
class SqLatticePercolation{
    // constants
    value_type  _length;
    value_type _max_number_of_bonds;
    value_type _max_number_of_sites;
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

//    double _occuption_probability {};
    // entropy
    long double _entropy{}, delta_H;
    long double _entropy_current{};
    long double _entropy_previous{};
    long double _largest_jump_entropy{}; // lrgest jump in entropy
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
    value_type delta_P; // jump of order parameter, i.e., number of bonds (or sites) in wrapping cluster
    value_type _wrapping_cluster_size, _wrapping_cluster_size_previous;
    value_type _index_last_modified_cluster{};  // id of the last modified cluster
    std::mt19937 _random_engine;
    value_type _random_state;

    void set_type(char t){type = t;} // setting percolation type
public:
    static constexpr const char* signature = "SqLatticePercolation";

    virtual ~SqLatticePercolation() = default;
    explicit SqLatticePercolation(value_type length);
    virtual void init();
    void setRandomState(size_t seed, bool generate_seed=true);
    value_type getRandomState();
    void reset();
    void set_cluster_measuring_unit(int i){
        std::cout << "Cluster_v3 measuring unit = " << ((i==0) ? "bond" : "site")
                  << " : line " << __LINE__ << std::endl;
    }

    bool occupy();
    value_type length() const { return _length;}
    value_type maxSites() const {return _max_number_of_sites;}
    value_type maxBonds() const { return _max_number_of_bonds;}

    /*********
     *  I/O functions
     */
    virtual void viewCluster();
    virtual void viewClusterExtended();
    void viewCluster(Index i);
    void viewCluster(BondIndex i);
    void viewClusterByIndex(int i);
    virtual void ckeckCluster();
    virtual void view_bonds(){
        _lattice.view_bonds();
    }
    virtual void viewLattice(){
        _lattice.view_sites();

    }

    /**
 * Also shows the cluster index of the sites
 */
    virtual void viewLatticeExtended(){
        _lattice.view_sites_extended();

    }

    /**
    * Displays group ids of sites in a matrix form
    */
    virtual void viewLatticeByID(){
        _lattice.view_sites_by_id();
        _lattice.view_bonds_by_id();
    }

    virtual void viewSiteByID(){
        _lattice.view_sites_by_id();
    }

    virtual void viewBondByID(){
        _lattice.view_bonds_by_id();
    }

    virtual void viewSiteByRelativeIndex(){
        _lattice.view_sites_by_relative_index();
    }
    virtual void viewBondByRelativeIndex(){
        _lattice.view_bonds_by_relative_index_v4();
    }

    virtual void viewByRelativeIndex(){
        _lattice.view_by_relative_index();
    }

    virtual void view(){
        _lattice.view();
    }

    void viewRandomized();


    virtual long double entropy() { return _entropy_current;}
    double entropy_by_site(); // for future convenience. // the shannon entropy. the full calculations. time consuming
    double entropy_by_bond(); // for future convenience. // the shannon entropy. the full calculations. time consuming
    double orderParameter();
    size_t numberOfcluster() const {return _cluster_count;}
    void jump();
    long double jump_entropy();
    value_type jump_wrapping_cluster();
    double largestEntropyJump()const { return _largest_jump_entropy;}
    double largestEntropyJump_pc()const { return _entropy_jump_pc;}


    void get_cluster_info(
            std::vector<value_type> &site,
            std::vector<value_type> &bond
    );

    char get_type() const {return type;} // get percolation type
    virtual value_type maxIterationLimit() {return _max_iteration_limit;};

    double get_relabeling_time() const {return time_relabel;}
    value_type relabeling_count() const {return _total_relabeling;}

    // tracking
//    virtual void track_largestCluster();
//    void track_numberOfBondsInLargestCluster(){}; // delete in future
//    virtual void track_numberOfSitesInLargestCluster(){}; // delete in future
    virtual void track_largestCluster(value_type base);


    virtual const std::vector<double> clusterSizeDistribution() const {std::cout << "not defined in SqLatticePercolation" << std::endl;return {};};  // 2019.06.17

    IndexRelative getRelativeIndex(Index root, Index site_new);

    size_t numberOfClusters() const {return number_of_clusters;}

    // must implement virtual functions
    virtual double occupationProbability() const = 0;

};




#endif //SITEPERCOLATION_PERCOLATION_H

