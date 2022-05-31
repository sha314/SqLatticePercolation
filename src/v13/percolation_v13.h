//
// Created by shahnoor on 2/15/21.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_H
#define SQLATTICEPERCOLATION_PERCOLATION_H


#include "lattice_v13.h"
#include "lattice_v14.h"
#include "cluster.h"
#include "../types.h"
#include <set>
#include <algorithm>
#include <random>
#include "status.h"
#include "../util/printer.h"

class Percolation_v13{
#ifdef DESTRUCTOR_CHECK
private:
    int random_ID = 0;
public:
    ~Percolation_v13(){
        std::cout << "~Percolation_v13() " << random_ID << std::endl;
    }
    Percolation_v13(Percolation_v13&)= default;
    Percolation_v13(Percolation_v13&&)= default;
    Percolation_v13& operator=(Percolation_v13 &)= default;
    Percolation_v13& operator=(Percolation_v13&&)= default;
#endif

public:
    int _length = 0;
    int _seed = 0;
    size_t wrapping_cluster_site_count_pc = 0;
    size_t wrapping_cluster_bond_count_pc = 0;

    value_type _random_state = 0;
    std::mt19937 _random_engine;

    Lattice_v14 lattice_ref;
    ClusterPool cluster_pool_ref;

    Percolation_v13(int length, int seed=-1){
        _length = length;
        _seed = seed;
//        lattice_ref = Lattice_v14(_length);
        lattice_ref.create_lattice(_length);
        cluster_pool_ref = ClusterPool();

#ifdef DESTRUCTOR_CHECK
        random_ID = rand();
        std::cout << "Percolation_v13() " << random_ID << std::endl;
#endif
    }

    virtual void reset() {
        lattice_ref.reset();
        cluster_pool_ref.reset();
        wrapping_cluster_site_count_pc = 0;
        wrapping_cluster_bond_count_pc = 0;
    }

    void viewLattice(int formatt=0);
    void viewLattice_clsuter_k(int cluster_gid){
        // view only the cluster of gid k in the lattice
        lattice_ref.viewLattice_clsuter_k(cluster_gid);
    }
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

    void test_lattice(){lattice_ref.test_id_index_conversion();}
};

class SitePercolation_v13: public Percolation_v13{
#ifdef DESTRUCTOR_CHECK
private:
    int random_ID = 0;
public:
    ~SitePercolation_v13(){
        std::cout << "~SitePercolation_v13() " << random_ID << std::endl;
    }
    SitePercolation_v13(SitePercolation_v13&)= default;
    SitePercolation_v13(SitePercolation_v13&&)= default;
    SitePercolation_v13& operator=(SitePercolation_v13 &)= default;
    SitePercolation_v13& operator=(SitePercolation_v13&&)= default;
#endif

    std::string signature = "SitePercolation";



    size_t site_count_pc = 0;
//    Site_v13 current_site;

    long cluster_count = -1;
    long largest_cluster_sz = 0;
    int largest_cluster_id = -1;
    std::vector<int> merging_cluster_ids;

// 0-th element = largest cluster
// 1-st element = 2nd largest cluster
// 2-rd element = 3rd largest cluster
    // std::vector<long> cluster_sizes_sorted;
    // std::vector<int> cluster_sizes_sorted_ids;
    

    // std::set<int> largest_2nd_cluster_ids; // 2nd largest cluster ids. As there could be more than one 2nd largest clusters
    // long largest_2nd_cluster_size=0; // 2nd largest clsuter

    
    
    std::vector<int> cluster_size_n; // list of cluster sizes
    std::vector<std::set<int>> cluster_size_n_ids; // corresponding ids

protected:
    P_STATUS status; // holds the current status of percolation.
    // if site occupation was successful or not... Helpful for L1, L2 percolation RSBD

    int occupied_site_count = 0;
    int current_idx = 0;
    std::vector<int> site_ids_indices;
    int selected_id =-1;
    double max_entropy = 0;
    double entropy_value = 0;
    bool after_wrapping = false;
    int wrapping_cluster_id = -1;

public:
    SitePercolation_v13() = default;
    explicit SitePercolation_v13(int length, value_type seed=0);

    virtual std::string get_signature(){return signature;}
    void init_clusters();
    void shuffle_indices() {
//        std::cout << "shuffle_indices not set up" << std::endl;
//        std::cout << "Before shuffle" << std::endl;
//        print_list(site_ids_indices);
        std::shuffle(site_ids_indices.begin(), site_ids_indices.end(), _random_engine);
//        std::cout << "after shuffle " << __LINE__<< std::endl;
//        print_list(site_ids_indices);
    }
    int get_length(){ return lattice_ref.get_length();}
    void reset() override;
    int get_neighbor_site(int central_site_id, int connecting_bond_id);
    std::vector<int> get_connected_sites(Site_v13 site, std::vector<int>& bond_neighbors);

    virtual P_STATUS select_site();
    bool place_one_site();
    void track_largest_cluster(int new_cluster);
    void track_largest_clusters_v4(int new_cluster);
    
    void entropy_subtract(std::vector<int>& bond_neighbors);
    void entropy_add(int new_cluster_id);

    double occupation_prob();
    double get_pc(){ return double(site_count_pc)/(_length*_length);}
    double get_wrapping_cluster_site_count_at_pc(){ return wrapping_cluster_site_count_pc;}
    double get_wrapping_cluster_bond_count_at_pc(){ return wrapping_cluster_bond_count_pc;}

    double entropy();
    double entropy_v1();
    double entropy_v2();
    long largest_cluster(){return largest_cluster_sz;}
    long largest_2nd_cluster(){return cluster_size_n[1];} // second largest clsuter
    long largest_3rd_cluster(){return cluster_size_n[2];} // 3rd largest clsuter
    double order_param_largest_clstr(){ return double(largest_cluster_sz) / lattice_ref.get_bond_count();}
    double order_param_wrapping();

    int merge_clusters(std::vector<int>& bond_neighbors);
    std::vector<int> uniqe_gid_bond_neighbors(std::vector<int>& bond_neighbors_ids);
    int merge_clusters_v2(std::vector<int>& bond_neighbor_ids);
    int merge_clusters_v3(std::vector<int>& bond_neighbor_ids);
    int merge_clusters_v4(std::vector<int>& bond_neighbor_ids);
    void relabel_relative_indices(int connecting_bond_id);
    bool detect_wrapping();
    Site_v14& get_current_site(); // TODO
    int find_root_cluster(const std::vector<int> &bond_gids) ;
    int find_bonds_connected_to_root_cluster(const std::vector<int> &bond_ids) ;

    value_type random_int(int min, int max){
        if (min <0 || max <0 || max <= min){
            std::cout << "Error. random_int. wrong range " << __FILE__ << " :" << __LINE__ << std::endl;
            std::cout << "min, max = " << min << "," << max << std::endl;
            return size_t(min);
        }
        return min + _random_engine() % (max - min);
    }

    size_t maxIterationLimit(){ return _length*_length;}
};


class SitePercolationL0_v13: public SitePercolation_v13{
#ifdef DESTRUCTOR_CHECK
private:
    int random_ID = 0;
public:
    ~SitePercolationL0_v13(){
        std::cout << "~SitePercolationL0_v13() " << random_ID << std::endl;
    }
    SitePercolationL0_v13(SitePercolationL0_v13&)= default;
//    SitePercolationL0_v13(SitePercolationL0_v13&&)= default;
    SitePercolationL0_v13& operator=(SitePercolationL0_v13 &)= default;
    SitePercolationL0_v13& operator=(SitePercolationL0_v13&&)= default;
#endif

protected:
    bool first_run;
    std::vector<double> occupation_prob_list;
    std::vector<double> entropy_list;
    std::vector<double> order_wrapping_list;
    std::vector<double> order_largest_list;

    std::vector<double> larest_1st_list;
    std::vector<double> larest_2nd_list;
    std::vector<double> larest_3rd_list;
public:

    SitePercolationL0_v13(int length, value_type seed=0, bool generate_seed=true);
    SitePercolationL0_v13(SitePercolationL0_v13&&) = default;
    std::string get_signature() override { return SitePercolation_v13::get_signature() + "L0_";};
    void reset() override;

    std::vector<double> get_entropy_array(){ return entropy_list;}
    std::vector<double> get_occupation_prob_array(){ return occupation_prob_list;}
    std::vector<double> get_order_param_wrapping_array(){ return order_wrapping_list;}
    std::vector<double> get_order_param_largest_array(){ return order_largest_list;}
    std::vector<double> get_1st_largest_cluster_array() {return larest_1st_list;}
    std::vector<double> get_2nd_largest_cluster_array() {return larest_2nd_list;}
    std::vector<double> get_3rd_largest_cluster_array() {return larest_3rd_list;}

    void run_once();
    void run_once_v2();

    void test_cluster();

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
