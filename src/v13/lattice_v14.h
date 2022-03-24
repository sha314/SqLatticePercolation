//
// Created by shahnoor on 3/21/22.
//

#ifndef SQLATTICEPERCOLATION_LATTICE_V14_H
#define SQLATTICEPERCOLATION_LATTICE_V14_H


#include "site_v13.h"
#include "bond_v13.h"

class Lattice_v14{
    int length = 0;
    long bond_count = 0;
    long site_count = 0;

    std::vector<Site_v14> site_matrix;
    std::vector<Bond_v14> bond_matrix;

    std::vector<int> bond_ids;
    std::vector<int> site_ids;

public:
    explicit Lattice_v14(int length=3){
        if (length < 0){
            std::cout << "length cannot be negative" << std::endl;
        }
        this->length = length;
        site_count = length*length;
        bond_count = 2*site_count;
//        std::cout << "reached " << __LINE__ << std::endl;
        site_matrix = std::vector<Site_v14>(site_count);
//        site_ids = std::vector<int>(site_count);
        bond_matrix = std::vector<Bond_v14>(bond_count);
//        bond_ids = std::vector<int>(bond_count);
//        std::cout << "reached " << __LINE__ << std::endl;
//        site_matrix.resize(site_count);
//        _site_ids.resize(site_count);
//        bond_matrix.resize(bond_count);
//        bond_ids.resize(bond_count);

        init_lattice();
//        std::cout << "reached " << __LINE__ << std::endl;
//        init_neighbor_ids();
//        std::cout << "reached " << __LINE__ << std::endl;
    }
    int get_length() { return length;}
    void reset();
    void init_lattice();
//    void init_neighbor_ids();
    long get_bond_count(){ return bond_count;}
    long get_site_count(){ return site_count;}
    int bottom_bond_of_site(int s0_index);
    int top_bond_of_site(int s0_index);
    int right_bond_of_site(int s0_index);
    int left_bond_of_site(int s0_index);
    int get_col_from_id(int so_index);
    int get_row_from_id(int so_index);

    std::vector<int> find_neighbor_bonds(int so_index);
    std::vector<int> get_neighbor_bonds(int so_index);
    std::vector<int> get_neighbor_sites(int so_index);
    std::vector<int> get_site_neighbor_of_site(int s0_index);
    std::vector<int> get_sites_for_wrapping_test(int s0_index);

    void print_bonds();
    void print_sites();

    Site_v14& get_site_by_id(int id);

    Bond_v14& get_bond_by_id(int id);

    void set_site_gid_by_id(int id, int gid);

    void set_bond_gid_by_id(int id, int gid);

    int get_site_gid_by_id(int id);


    int get_bond_gid_by_id(int id);

    Site_v14 get_site_by_index(size_t row, size_t col) ;


    /**
     *
     * @param row
     * @param col
     * @param hv_flag : 0 is for horizontal and 1 is for vertical bonds
     * @return
     */
    Bond_v14 get_bond_by_index(int row, int col, int hv_flag=0);

    std::vector<int> get_site_id_list(){ return site_ids;}
    std::vector<int> get_bond_id_list(){ return bond_ids;}

    void view(int formatt=0);
    void view_relative_index();
    void view_site_gids();
    void print_row_separator(int str_sz=10);

    void viewLattice_clsuter_k(int cluster_gid);
    void print_site_bond_list();

    std::vector<std::string> get_row_str(int row, int format=0);
    std::vector<std::string> get_row_v_str(int row, int format=0);

    void init_relative_index(int id);
    void set_relative_index(int id, RelativeIndex_v13 relative_index);

    std::vector<int> get_all_neighbor_sites(int central_site_id);

    void test_id_index_conversion();

    void check_if_it_is_not_equal(int id, int s_index, int line) const;

    bool isRemovable(int site, int type);
};

#endif //SQLATTICEPERCOLATION_LATTICE_V14_H
