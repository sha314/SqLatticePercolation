//
// Created by shahnoor on 2/15/21.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_H
#define SQLATTICEPERCOLATION_PERCOLATION_H


#include "lattice_v13.h"
#include "cluster.h"
#include "../types.h"
#include <set>
#include <algorithm>
#include <random>

class Percolation_v13{

public:
    int _length = 0;
    int _seed = 0;
    value_type _random_state = 0;
    std::mt19937 _random_engine;

    Lattice_v13 lattice_ref;
    ClusterPool cluster_pool_ref;

    Percolation_v13(int length, int seed=-1){
        _length = length;
        _seed = seed;
        lattice_ref = Lattice_v13(_length);
        cluster_pool_ref = ClusterPool();
    }

    virtual void reset() {
        lattice_ref.reset();
        cluster_pool_ref.reset();
    }

    void viewLattice(int formatt=0);
    void viewCluster(int view_mode) {
        cluster_pool_ref.view(view_mode);
    }

    std::vector<int>  get_bond_gids(std::vector<int>& bond_ids);
    std::vector<int>  get_site_gids(std::vector<int>& site_ids);
    Index_v13 wrapping_correction_relative_index(Index_v13 delta_X);
    RelativeIndex_v13 wrapping_correction_relative_index(RelativeIndex_v13 delta_X);

    RelativeIndex_v13 get_relative_index(int central_site_id, int neighbor_site_id);
    RelativeIndex_v13 get_change_in_relative_index(RelativeIndex_v13 old_relative_index, RelativeIndex_v13 new_relative_index);

    void setRandomState(size_t seed, bool generate_seed=true);
    value_type getRandomState() { return _random_state;}
};

class SitePercolation_v13: public Percolation_v13{
    std::string signature = "SitePercolation";
    std::vector<int> site_ids_indices;
    int current_idx = 0;
    Site_v13 current_site;
    int selected_id =-1;
    long cluster_count = -1;
    long largest_cluster_sz = 0;
    int largest_cluster_id = -1;
    double max_entropy = 0;
    double entropy_value = 0;
    bool after_wrapping = false;
    int wrapping_cluster_id = -1;

public:
    SitePercolation_v13() = default;
    explicit SitePercolation_v13(int length, int seed=-1);

    virtual std::string get_signature(){return signature;}
    void init_clusters();
    void shuffle_indices() {
        std::cout << "shuffle_indices not set up" << std::endl;
        std::shuffle(site_ids_indices.begin(), site_ids_indices.end(), _random_engine);
    }
    int get_length(){ return lattice_ref.get_length();}
    void reset() override;
    int get_neighbor_site(int central_id, int connecting_bond_id);
    std::vector<int> get_connected_sites(Site_v13 site, std::vector<int>& bond_neighbors);

    bool select_site();
    bool place_one_site();
    void track_largest_cluster(int new_cluster);
    void entropy_subtract(std::vector<int>& bond_neighbors);
    void entropy_add(int new_cluster_id);

    double occupation_prob();

    double entropy();
    double entropy_v1();
    double entropy_v2();
    long largest_cluster(){return largest_cluster_sz;}
    double order_param_largest_clstr(){ return double(largest_cluster_sz) / lattice_ref.get_bond_count();}
    double order_param_wrapping();

    int merge_clusters(std::vector<int>& bond_neighbors);
    std::vector<int> uniqe_gid_bond_neighbors(std::vector<int>& bond_neighbors);
    int merge_clusters_v2(std::vector<int>& bond_neighbors);
    void relabel_relative_indices(int connecting_bond);
    bool detect_wrapping();
};


class SitePercolationL0_v13: public SitePercolation_v13{
    bool first_run;
    std::vector<double> occupation_prob_list;
    std::vector<double> entropy_list;
    std::vector<double> order_wrapping_list;
    std::vector<double> order_largest_list;
public:

    SitePercolationL0_v13(int length, int seed=-1, bool generate_seed=true);
    SitePercolationL0_v13(SitePercolationL0_v13&&) = default;
    std::string get_signature() override { return SitePercolation_v13::get_signature() + "L0_";};
    void reset() override;

    std::vector<double> get_entropy_array(){ return entropy_list;}
    std::vector<double> get_occupation_prob_array(){ return occupation_prob_list;}
    std::vector<double> get_order_param_wrapping_array(){ return order_wrapping_list;}
    std::vector<double> get_order_param_largest_array(){ return order_largest_list;}

    void run_once();
};


template<class tType>
void remove_vector_element(std::vector<tType> &vec, tType elem, bool once=true) {
    for(int i=0; i < vec.size();++ i){
        if(vec[i] == elem){
            vec.erase(vec.begin() + i);
            if(once) break;
        }
    }
}

template<class tType>
bool find_elm(std::vector<tType> &vec, tType elem) {
    for(int i=0; i < vec.size();++ i){
        if(vec[i] == elem){
            return true;
        }
    }
    return false;
}

template<class tType>
void print_vectors(std::vector<tType> &vec, const std::string &end_str) {
    for(auto b: vec){
        std::cout << b << ",";
    }
    std::cout << end_str;


}


#endif //SQLATTICEPERCOLATION_PERCOLATION_H
