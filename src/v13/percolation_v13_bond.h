//
// Created by shahnoor on 2/15/21.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_BOND_H
#define SQLATTICEPERCOLATION_PERCOLATION_BOND_H


#include "lattice_v13.h"
#include "cluster.h"
#include "../types.h"
#include <set>
#include <algorithm>
#include <random>
#include "status.h"
#include "../util/printer.h"
#include "percolation_v13.h"



class BondPercolation_v13: public Percolation_v13{
    std::string signature = "BondPercolation_v13";

    size_t bond_count_pc = 0;
//    Site_v13 current_site;

    long cluster_count = -1;
    long largest_cluster_sz = 0;
    int largest_cluster_id = -1;

protected:
    P_STATUS status; // holds the current status of percolation.
    // if site occupation was successful or not... Helpful for L1, L2 percolation RSBD

    int occupied_bond_count = 0;
    int current_idx = 0;
    // std::vector<int> site_ids_indices;
    std::vector<int> bond_ids_indices;

    int selected_id =-1;
    double max_entropy = 0;
    double entropy_value = 0;
    bool after_wrapping = false;
    int wrapping_cluster_id = -1;

    double sum_cluster_size_squared=0;
    double sum_cluster_size=0;


    bool first_run;
    std::vector<double> occupation_prob_list;
    std::vector<double> entropy_list;
    std::vector<double> order_wrapping_list;
    std::vector<double> order_largest_list;
    std::vector<double> mean_cluster_sz_list;

public:
    BondPercolation_v13() = default;
    ~BondPercolation_v13();
    explicit BondPercolation_v13(int length, value_type seed=0,  bool generate_seed=true);

    virtual std::string get_signature(){return signature;}
    void init_clusters();
    void shuffle_indices() {
//        std::cout << "shuffle_indices not set up" << std::endl;
    //    std::cout << "Before shuffle" << std::endl;
    //    print_list(site_ids_indices);
        std::shuffle(bond_ids_indices.begin(), bond_ids_indices.end(), _random_engine);
    //    std::cout << "after shuffle " << __LINE__<< std::endl;
    //    print_list(site_ids_indices);
    }
    int get_length(){ return lattice_ref.get_length();}
    void reset() override;
    int get_neighbor_site(int central_site_id, int connecting_bond_id);
    std::vector<int> get_connected_sites(Site_v13 site, std::vector<int>& bond_neighbors);

    virtual P_STATUS select_bond();
    bool place_one_bond();
    void track_largest_cluster(int new_cluster);
    void entropy_subtract(std::vector<int>& bond_neighbors);
    void entropy_add(int new_cluster_id);

    double occupation_prob();
    double get_pc(){ return double(bond_count_pc)/(_length*_length);}
    double get_wrapping_cluster_site_count_at_pc(){ return wrapping_cluster_site_count_pc;}
    double get_wrapping_cluster_bond_count_at_pc(){ return wrapping_cluster_bond_count_pc;}

    double entropy();
    double entropy_v1();
    double entropy_v2();
    long largest_cluster(){return largest_cluster_sz;}
    double order_param_largest_clstr(){ return double(largest_cluster_sz) / lattice_ref.get_bond_count();}
    double order_param_wrapping();

    double get_mean_cluster_size();
    double get_mean_cluster_size_v2();

    int merge_clusters(std::vector<int>& bond_neighbors);
    std::vector<int> uniqe_gid_site_neighbors(std::vector<int>& bond_neighbors_ids);
    int merge_clusters_v2(std::vector<int>& bond_neighbor_ids);
    int merge_clusters_v3(std::vector<int>& bond_neighbor_ids);
    int merge_clusters_v4(std::vector<int>& bond_neighbor_ids);
    void relabel_relative_indices(int connecting_bond_id);
    bool detect_wrapping();
    Bond_v13& get_current_bond(); // TODO
    int find_root_cluster(const std::vector<int> &bond_gids) ;
    int find_sites_connected_to_root_cluster(const std::vector<int> &bond_ids) ;

    value_type random_int(int min, int max){
        if (min <0 || max <0 || max <= min){
            std::cout << "Error. random_int. wrong range " << __FILE__ << " :" << __LINE__ << std::endl;
            std::cout << "min, max = " << min << "," << max << std::endl;
            return size_t(min);
        }
        return min + _random_engine() % (max - min);
    }

    size_t maxIterationLimit(){ return _length*_length;}

    /**
     * @brief view site id list in curly brace {}.
     *  The unoccupied sites are inside square bracket [].
     * The site with angle bracket is the currenlty selected site <>.
     * 
     */
    void view_site_id_list();
    // {
    //     std::cout << "view_site_id_list() {";
    //     auto got_selected_site_flag = false;
    //     for(int i=0; i< site_ids_indices.size(); ++i){
    //         if (i == current_idx){
    //             if(got_selected_site_flag){
    //                 // 0 for background Color(Black)
    //                 // 4 for text color(Red)
    //                 system("Color 04");
    //                 std::cerr << "How can the current site be outside the unoccupied sites?" << std::endl;
    //                 std::cout << "Probably you are calling view_site_id_list without occupying any new sites" << std::endl;
    //                 // exit(0);
    //             }
    //             std::cout << "[";
    //         }
    //         auto a = site_ids_indices[i];
    //         if(a == selected_id){
    //             std::cout <<"<" << a << ">,";
    //             got_selected_site_flag = true;
    //             continue;
    //         }
    //         std::cout << a << ",";
    //     }
    //     // for(auto a: site_ids_indices){
    //     //     if(a == selected_id){
    //     //         std::cout <<"<" << a << ">,";
    //     //         continue;
    //     //     }
    //     //     std::cout << a << ",";
    //     // }
    //     std::cout << "]}" << std::endl;
    // }

    /**
 * Not normalized
 * @return an array : index of the array is the cluster size.
 *                    value of the array is the number of times that cluster appears
 */
const std::vector<double> clusterSizeDistribution();
};

#endif //SQLATTICEPERCOLATION_PERCOLATION_H
