//
// Created by shahnoor on 2/15/21.
//

#ifndef SQLATTICEPERCOLATION_CLUSTER_H
#define SQLATTICEPERCOLATION_CLUSTER_H

#include <vector>
#include <string>
#include "lattice_v13.h"

class OneCluster{
    std::vector<int> _site_ids, _bond_ids;
    int _gid = -1;
    int _id = -1;

public:
    OneCluster() = default;

    std::string get_str();
    bool is_empty(){ return _site_ids.size() == 0 and _bond_ids.size() ==0;}
    void set_id(int id);
    void set_gid(int gid);
    int get_id(){ return _id;}
    int get_gid(){ return _gid;}
    void add_site(int site_id);
    void add_sites(std::vector<int>& site_ids);
    void add_bond(int bond_id);
    void add_bonds(std::vector<int>& bond_ids);

    std::vector<int>& get_sites() { return _site_ids;}
    std::vector<int>& get_bonds() { return _bond_ids;}

    int get_site_count() { return _site_ids.size();}
    int get_bond_count() { return _bond_ids.size();}

    void clear(){_site_ids.clear();_bond_ids.clear();_gid=-1;_id=-1;}

    int view(int view_mode=0);

};

class ClusterPool {
    std::vector<OneCluster> _cluster_list;
    int cluster_id = 0;

public:
    ClusterPool() = default;

    void reset() {
        int cluster_id = 0;
        _cluster_list.clear();
    }

    int get_cluster_bond_count(int id) {
        return _cluster_list[id].get_bond_count();
    }

    int get_cluster_site_count(int id) {
        return _cluster_list[id].get_site_count();
    }


    void add_sites(int index, std::vector<int> &site_ids) {
        _cluster_list[index].add_sites(site_ids);
    }

    void add_sites(int index, int site_ids) {
        _cluster_list[index].add_site(site_ids);
    }

    void add_bonds(int index, std::vector<int> &bond_ids){
        _cluster_list[index].add_bonds(bond_ids);
    }

    std::vector<int>  get_sites(int index) {
        return _cluster_list[index].get_sites();
    }

    std::vector<int> get_bonds(int index) {
        return _cluster_list[index].get_bonds();
    }

    OneCluster& get_cluster(int cluster_id) {
        if (cluster_id >= _cluster_list.size()) {
            std::cout << "Cluster does not exists" << std::endl;
        }

        return _cluster_list[cluster_id];
    }

    void clear_cluster(int cluster_id) {
        _cluster_list[cluster_id].clear();
    }

    void create_new_cluster(std::vector<int> &site_ids, std::vector<int> &bond_ids, Lattice_v13 &lattice_ref);
    void create_new_cluster(int site_id, int bond_id, Lattice_v13 &lattice_ref);
    void merge_cluster_with(int cluster_A_id, int cluster_B_id, Lattice_v13& lattice_ref);
    void view(int view_mode=0);

};

#endif //SQLATTICEPERCOLATION_CLUSTER_H
