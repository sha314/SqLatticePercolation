//
// Created by shahnoor on 2/15/21.
//

#include "cluster.h"

using namespace std;

void ClusterPool::view(int view_mode) {
    cout << "View cluster < BEGIN" << endl;
//# print("self.cluster_list ", self.cluster_list)
    int counter = 0;
    for (auto clstr : _cluster_list) {
//# print("clstr ", clstr)
        counter += clstr.view(view_mode);
    }

    cout << "\n View cluster END >" << endl;
    cout << "Number of clusters " << counter << endl;

}

void ClusterPool::create_new_cluster(std::vector<int> &site_ids, std::vector<int> &bond_ids, Lattice_v13 &lattice_ref) {
//# print("method : create_new_cluster")
    auto clsstr = OneCluster();
    clsstr.add_sites(site_ids);
    clsstr.add_bonds(bond_ids);
//# print("setting gid ")
    clsstr.set_gid(cluster_id);
    clsstr.set_id(cluster_id);
//# print("bond_ids ", bond_ids)
    for(auto ss: site_ids){
        lattice_ref.set_site_gid_by_id(ss, cluster_id); // reassign group id
    }


    for(auto bb: bond_ids) {
        lattice_ref.set_bond_gid_by_id(bb, cluster_id); // reassign group id
    }
    cluster_id += 1;
//# print("clsstr ")
//# print(clsstr)
    _cluster_list.push_back(clsstr);
}

void ClusterPool::create_new_cluster(int site_id, int bond_id, Lattice_v13 &lattice_ref) {
//# print("method : create_new_cluster")
    auto clsstr = OneCluster();

//# print("setting gid ")
    clsstr.set_gid(cluster_id);
    clsstr.set_id(cluster_id);
//# print("bond_ids ", bond_ids)

    if(site_id >= 0){
        clsstr.add_site(site_id);
        lattice_ref.set_site_gid_by_id(site_id, cluster_id); // reassign group id
    }
    if(bond_id >= 0) {
        clsstr.add_bond(bond_id);
        lattice_ref.set_bond_gid_by_id(bond_id, cluster_id); // reassign group id
    }

    cluster_id += 1;
//# print("clsstr ")
//# print(clsstr)
    _cluster_list.push_back(clsstr);
}

/**
 *
 * @param cluster_A_id : the cluster. group id of cluster A will persist.
 * @param cluster_B_id : the other cluster that will be merged to this cluster
 * @param lattice_ref  : so that it gid of sites and bonds can be modified here
 */
void ClusterPool::merge_cluster_with(int cluster_A_id, int cluster_B_id, Lattice_v13 &lattice_ref) {
    cout << "merge_cluster_with" << endl;
    int gid = _cluster_list[cluster_A_id].get_gid();
    cout << "cluster " << cluster_A_id << " gid " << gid << endl;
    cout << "cluster " << cluster_B_id << " gid " << endl;
    cout << "cauing error " << __FILE__ << " " << __LINE__ << endl;
    auto bondstmp = _cluster_list[cluster_B_id].get_bonds();
    cout << "bonds size " << bondstmp.size() << " why?" << endl;
    for (auto bb : bondstmp){ //
        lattice_ref.set_bond_gid_by_id(bb, gid);
    }

    for (auto ss : _cluster_list[cluster_B_id].get_sites()){
//# print("todo: site ", ss, " gid => ", gid)
        lattice_ref.set_site_gid_by_id(ss, gid);
//      int tmp = lattice_ref.get_site_by_id(ss).get_gid();
//# print("done: site ", ss, " gid => ", tmp)
    }
//    cout << "SIGSEGV : line " << __LINE__ << endl;
    std::vector<int> bonds = _cluster_list[cluster_B_id].get_bonds();
    _cluster_list[cluster_A_id].add_bonds(bonds);

//# print("before ", self.cluster_list[cluster_A_id].site_ids)
//# print("adding ", self.cluster_list[cluster_B_id].site_ids)
    std::vector<int> sites = _cluster_list[cluster_B_id].get_sites();
    _cluster_list[cluster_A_id].add_sites(sites);
//# print("after ", self.cluster_list[cluster_A_id].site_ids)
    _cluster_list[cluster_B_id].clear(); // clear redundent cluster
}

/**
 *
 * @param view_mode : if 0 then only culsters with non-zero sites and bonds will be shown.
 * @return
 */
int OneCluster::view(int view_mode) {
    if (is_empty()) {
        return 0;
    }
    if(view_mode == 0){
        if(_site_ids.empty() || _bond_ids.empty()) {
            return 1;
        }
    }
    cout <<  "cluster[" << setw(3) <<  _id << "] (gid " << _gid << ") :{" << endl;
    cout << "  sites (" << _site_ids.size() << ") {";
    for(auto aa: _site_ids){
        cout << aa << ",";
    }
    cout << "}" << endl;

    cout << "  bonds (" << _bond_ids.size() << ") {";
    for(auto aa: _bond_ids){
        cout << aa << ",";
    }
    cout << "}" << endl;

    return 1;
}

void OneCluster::add_bonds(std::vector<int> &bond_ids) {
    _bond_ids.insert(_bond_ids.end(), bond_ids.begin(), bond_ids.end());
}

void OneCluster::add_bond(int bond_id) {
    _bond_ids.push_back(bond_id);
}

void OneCluster::add_sites(std::vector<int> &site_ids) {
    _site_ids.insert(_site_ids.end(), site_ids.begin(), site_ids.end());
}

void OneCluster::add_site(int site_id) {
    _site_ids.push_back(site_id);
}

void OneCluster::set_gid(int gid) {
_gid = gid;
}

void OneCluster::set_id(int id) {
_id = id;
}

std::string OneCluster::get_str() {
    stringstream ret_str;
    ret_str << "[" << _id << "] (gid " << _gid << "): ";
    ret_str << _site_ids.size() << " sites and " << _bond_ids.size() << " bonds";
    return ret_str.str();
}
