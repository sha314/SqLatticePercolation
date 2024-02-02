//
// Created by shahnoor on 2/15/21.
//

#include <cmath>
#include <set>
#include <random>
#include "percolation_v13_bond.h"

using namespace std;

BondPercolation_v13::~BondPercolation_v13(){}

BondPercolation_v13::BondPercolation_v13(int length, value_type seed, bool generate_seed) : Percolation_v13(length, seed) {
    cout << "class : BondPercolation_v13" << endl;
    init_clusters();

    bond_ids_indices = lattice_ref.get_bond_id_list();
//    cout << "{";
//    for(auto s: site_ids_indices){
//        cout << s << ",";
//    }
//    cout << "}" << endl;
    int current_idx = 0;
    setRandomState(seed, generate_seed);
    shuffle_indices();
    selected_id = -1;
    cluster_count = lattice_ref.get_bond_count();
    largest_cluster_sz = 0;
    largest_cluster_id = -1;
    max_entropy = log(lattice_ref.get_bond_count());
    entropy_value = max_entropy;
    after_wrapping = false;
    wrapping_cluster_id = -1;
    
}

void BondPercolation_v13::init_clusters() {
//    cluster_pool_ref.reset();
    auto sites = lattice_ref.get_site_id_list();
    for (auto ss: sites){
        cluster_pool_ref.create_new_cluster(ss, -1, lattice_ref);
    }
    // all clusters have length one initially. and no need to exclude largest cluster initially
    sum_cluster_size_squared = _length*_length;
    sum_cluster_size = _length*_length;
}

void BondPercolation_v13::reset() {
    Percolation_v13::reset();

    bond_ids_indices = lattice_ref.get_bond_id_list();;
    current_idx = 0;
    occupied_bond_count = 0;
    bond_count_pc = 0;
    selected_id =-1;
    cluster_count = lattice_ref.get_bond_count();
    largest_cluster_sz = 0;
    largest_cluster_id = -1;
    entropy_value = max_entropy;
    after_wrapping = false;
    wrapping_cluster_id = -1;



//# print("Initial entropy ", self.entropy_value)
    init_clusters();
    shuffle_indices();

    //    int ll = get_length();
//    int l_squared = ll*ll;

    if (first_run) {
        occupation_prob_list.clear();
//        occupation_prob_list.resize(l_squared);
    }

    entropy_list.clear();
    order_wrapping_list.clear();
    order_largest_list.clear();
    mean_cluster_sz_list.clear();


//    entropy_list.resize(l_squared);
//    order_wrapping_list.resize(l_squared);
//    order_largest_list.resize(l_squared);
}

double BondPercolation_v13::entropy_v1() {
    double H = 0;
    int empty_count = 0;
    for(int i=0; i < cluster_count;++i){
        double s_count = cluster_pool_ref.get_cluster_site_count(i);
        if (s_count==0) {
            empty_count += 1;
            continue;
        }
        double mu = s_count / lattice_ref.get_site_count();
        double log_mu = log(mu);
        H += mu * log_mu;
    }
//    cout << "empty cluster count " << empty_count << endl;
#ifdef UNIT_TEST
    if (cluster_count <= 0){
        cout << "Error. Cluster count can't be zero : " << __LINE__ << " file " << __FILE__ << endl;
        cout << "cluster_count " << cluster_count << endl;
    }
    if (empty_count == cluster_count){
        cout << "Error. All clusters are empty?! : " << __LINE__ << " file " << __FILE__ << endl;
    }
#endif
    return -H;
}

/**
 * @brief see ## 2# in Docs/definitions.md
 * 
 * @return double 
 */
double BondPercolation_v13::order_param_wrapping() {
    if (after_wrapping) {
//# print("wrapping cluster id ", self.wrapping_cluster_id)
        double count = cluster_pool_ref.get_cluster_site_count(wrapping_cluster_id);
#ifdef UNIT_TEST
    if(count <= 0){
        cout << "Error . After wrapping, number of bond in the wrapping cluster cannot be zero " << count << endl;
        cout << "wrapping_cluster_id " << wrapping_cluster_id << endl;
        cluster_pool_ref.view(0);
        exit(-1);
    }
#endif
        return count / lattice_ref.get_site_count();
    }
    return 0;
}

/**
 * @brief see ## 1# in Docs/definitions.md
 * 
 * 
 * doing the sum every time is time consuming. use v2 for that
 * @return double 
 */
double BondPercolation_v13::get_mean_cluster_size(){
    double sum_s = 0;
    double sum_s2 = 0; // sum of s squared
    int max_clusters = lattice_ref.get_site_count();
    for(int i=0; i<max_clusters; ++i){
        double s = cluster_pool_ref.get_cluster_site_count(i);
        // sum_s += s;
        sum_s2 += s*s;
    }
    // cout << "total bonds " << sum_s << endl;
    sum_s2 -= largest_cluster_sz*largest_cluster_sz;
    // sum_s -= largest_cluster_sz;
    sum_s = lattice_ref.get_site_count() - largest_cluster_sz;
    return sum_s2/sum_s;
}

/**
 * @brief Largest cluster is not excluded from the private variables. that is done here.
 * 
 * @return double 
 */
double BondPercolation_v13::get_mean_cluster_size_v2(){
    auto a = sum_cluster_size_squared - largest_cluster_sz*largest_cluster_sz;
    auto b = sum_cluster_size - largest_cluster_sz;
    if(b==0) {
        // cout << "b is zero" << endl;
        return 0;
    }
    return a/b;
}

int BondPercolation_v13::get_neighbor_site(int central_site_id, int connecting_bond_id) {
//# central_id = current_site.get_id()
//# print("central site id : ", central_id)

    auto connected_sites = lattice_ref.get_neighbor_sites(connecting_bond_id);
//# print("connected ", connected_sites)
//    for (int i=0; i < connected_sites.size(); ++i){
//        if (connected_sites[i] == central_id){
//            connected_sites.erase(connected_sites.begin() + i);
//            break;
//        }
//    }
    remove_vector_element(connected_sites, central_site_id);
    if(connected_sites.size() != 1){
        cout << "Number of neighbors must be exactly 1 : get_connected_sites()" << endl;
    }
    int found_neighbor = connected_sites[0];

#ifdef UNIT_TEST
// test if the found neighbor is actually a neighbor
    auto siteA = lattice_ref.get_site_by_id(central_site_id).get_index();
    auto siteB = lattice_ref.get_site_by_id(found_neighbor).get_index();
    auto bondAB = lattice_ref.get_bond_by_id(connecting_bond_id);
//    cout << "A, B =connected by bond > AB. " << siteA.get_str()
//         << ", " << siteB.get_str()
//         << " => AB " << bondAB.get_str() << endl;
    int dRow = abs(siteA.row() - siteB.row());
    int dCol = abs(siteA.col() - siteB.col());
    if((dRow+1) == _length){
        dRow = 1;
    }
    if((dCol+1) == _length){
        dCol = 1;
    }
//    cout << "dRow, dCol " << dRow << ", " << dCol << endl;
    if(dRow == 1 && dCol == 0){
//        cout << "Ok" << endl;
    }
    else if(dRow == 0 && dCol == 1){
//        cout << "Ok" << endl;
    }
    else{
        cout << "Error not a neighbor" << endl;
        cout << "lattice length " << _length << endl;
        exit(-1);
    }
#endif

//# print("neighbor site ids ", neighbor_sites)
    return found_neighbor ;
//    pass
}

std::vector<int> BondPercolation_v13::get_connected_sites(Site_v13 site, std::vector<int> &bond_neighbors) {
//# print("central site index : ", site.get_index())
    auto central_id = site.get_id();
//# print("central site id : ", central_id)
    vector<int> neighbor_sites;
    for(auto bb : bond_neighbors) {
        auto sb2 = lattice_ref.get_bond_by_id(bb);
        auto connected_sites = sb2.get_connected_sites();
//# print("connected ", connected_sites)
        remove_vector_element(connected_sites, central_id);
        if (connected_sites.size() > 1) {
            cout << "Number of neighbors cannot exceed 2 : get_connected_sites()" << endl;
        }
        neighbor_sites.push_back(connected_sites[0]);
    }
//# print("neighbor site ids ", neighbor_sites)
    return neighbor_sites;
}

P_STATUS BondPercolation_v13::select_bond() {
    // cout << "SitePercolation_v13::select_site() " << endl;
    if (current_idx >= lattice_ref.get_bond_count()) {
//        cout << "No sites to occupy" << endl;
        return P_STATUS::EMPTY_SITE_LIST;
    }
    selected_id = bond_ids_indices[current_idx];
    current_idx += 1;
//    current_site = lattice_ref.get_site_by_id(selected_id);
//    cout << "selected id " << selected_id << " site " << current_site.get_str() << endl;
//    cout << "selected id " << selected_id << " site " << get_current_site().get_str() << endl;
    return P_STATUS::SUCESS;

}

bool BondPercolation_v13::place_one_bond() {
//    cout << "************************ place_one_site. count " << current_idx << endl;
    status = select_bond();
    if(status == P_STATUS::SUCESS) {
        // view_site_id_list();   // selected site will be inside angle bracket <>

//        cout << "selected site ", self.current_site.get_index(), " id ", self.current_site.get_id())

        /**
         * @brief 
         * Since it's a bond, we can just initialize the relative index of one of the 
         * connected sites
         */
        // lattice_ref.init_relative_index(selected_id);  // initialize        relative index
        auto site_neighbors = get_current_bond().get_connected_sites();
        for( auto a: site_neighbors){
            cout << a << ", ";
        }
        cout << endl;
//# site_neighbors = self.get_connected_sites(self.current_site, bond_neighbors)
        exit(0);
        // entropy_subtract(site_neighbors);

        // auto merged_cluster_index = merge_clusters_v4(site_neighbors);

        // track_largest_cluster(merged_cluster_index);
        // entropy_add(merged_cluster_index);

//# self.lattice_ref.set_site_gid_by_id(selected_id, merged_cluster_index)
//# self.cluster_pool_ref.add_sites(merged_cluster_index, selected_id)
        occupied_bond_count += 1;

#ifdef UNIT_TEST
        if ( site_neighbors.size() != 4){
            cout << "Number of bond neighbor must be 4 : " << __LINE__ << endl;
            exit(-1);
        }
#endif
        return true;
    }
    else if(status == P_STATUS::EMPTY_SITE_LIST){
        return false;
    }

    return true; // to keep the loop running
}

void BondPercolation_v13::track_largest_cluster(int new_cluster) {
    auto new_size = cluster_pool_ref.get_cluster_bond_count(new_cluster);
//# self.cluster_pool_ref.get_cluster_site_count(new_cluster)
    if (new_size > largest_cluster_sz) {
        largest_cluster_id = new_cluster;
        largest_cluster_sz = new_size;
    }
}

void BondPercolation_v13::entropy_subtract(std::vector<int> &sites) {
//# print("entropy_subtract")
    // cout << "{"; for(auto a: bonds) {cout << a << ",";} cout << "} vs " << endl;
    // auto bonds = lattice_ref.get_neighbor_bonds(selected_id);
    // cout << "{"; for(auto a: bonds) {cout << a << ",";} cout << "}" << endl;

//# print(self.current_site, " neighbors => ", sites)
    set<int> gids;
    for(auto ss: sites) {
        auto gid = lattice_ref.get_site_gid_by_id(ss);
        if (gid == -1)        continue;
        gids.insert(gid);
    }
//# print("gids ", gids)
    double H = 0, mu;
    double bc = lattice_ref.get_site_count();
    for (auto gg : gids) {
        auto s_count = cluster_pool_ref.get_cluster_site_count(gg);
        if (s_count == 0)        continue;

        mu = s_count / bc;
        H += mu * log(mu);

        /*****For Mean Cluster size. begin*****/
        
        sum_cluster_size_squared -= s_count*s_count;
        sum_cluster_size -= s_count;

        /*****Ended the quest for mean cluster size*/
    }
//# print("before ", self.entropy_value)
    entropy_value += H;
//# print("after ", self.entropy_value)

}

void BondPercolation_v13::entropy_add(int new_cluster_id) {
//# print("entropy_add")
    double b_count = cluster_pool_ref.get_cluster_site_count(new_cluster_id);
    double mu = b_count / lattice_ref.get_site_count();
//# print("before ", self.entropy_value)
    entropy_value -= mu*log(mu);
//# print("after ", self.entropy_value)


    /*****For Mean Cluster size. begin*****/
    
    sum_cluster_size_squared += b_count*b_count;
    sum_cluster_size += b_count;

    /*****Ended the quest for mean cluster size*/
}

double BondPercolation_v13::entropy() {
//    return entropy_v1();
    return entropy_v2();
}

double BondPercolation_v13::entropy_v2() {
    return entropy_value;
}

/**
 *
 * @param connecting_bond_id
 */
void BondPercolation_v13::relabel_relative_indices(int connecting_bond_id) {
//    cout << "SitePercolation_v13::relabel_relative_indices" << endl;
    auto bond = lattice_ref.get_bond_by_id(connecting_bond_id);
    auto bbg = bond.get_gid();
//    auto central_site = current_site.get_id();
    auto central_site = get_current_bond().get_id();
    auto neighbor_site = get_neighbor_site(central_site, bond.get_id());
//    cout <<"central_site " << central_site << " and neighbor_site " << neighbor_site << endl;
//    cout << "sites of cluster " << bbg << " will be relabeled" << endl;
    auto sites_to_relabel = cluster_pool_ref.get_sites(bbg);
//# print("sites_to_relabel ", sites_to_relabel)
//# relabel neighbor according to central site
    if (sites_to_relabel.empty()) {
//        cout << "sites_to_relabel.empty()" << endl;
        return;
    }
    auto old_relative_idx = lattice_ref.get_site_by_id(neighbor_site).get_relative_index();
    auto new_relative_idx = get_relative_index(central_site, neighbor_site);
//# if the BBB lines are commented then it sould not affect the result. so why extra lines;
//    # self.lattice_ref.set_relative_index(neighbor_site, new_relative_idx)  # BBB
//        # then relabel all sites belonging to the cluster according to the neighbor

    if (lattice_ref.get_site_gid_by_id(central_site) >= 0){
//# old_relative_index = self.get_relative_index(central_site, self.selected_id)
        auto change = get_change_in_relative_index(old_relative_idx, new_relative_idx);
        RelativeIndex_v13 changeR = RelativeIndex_v13(change);
//        cout << "change " << change.get_str() << endl;
//        cout << "old_relative_index "  << old_relative_idx.get_str() << endl;
        for (auto ss : sites_to_relabel) {
//# if ss == neighbor_site:  # BBB
//#     print("already got relabeled") # BBB
//#     continue # BBB
            RelativeIndex_v13 ss_relative_index = lattice_ref.get_site_by_id(ss).get_relative_index();
//# change = self.get_change_in_relative_index(old_relative_index, ss_relative_index)
//# print("change ", change)
//# print("new_relative_index type ", type(ss_relative_index))
//            cout << "relative index before : " << ss_relative_index.get_str() << endl;
            auto temp = ss_relative_index + changeR;
            ss_relative_index = RelativeIndex_v13(temp) ;
//            cout << "relative index after  : " << ss_relative_index.get_str() << endl;
//            cout << "new_relative_index " << new_relative_index << endl;
            lattice_ref.get_site_by_id(ss).set_relative_index(ss_relative_index);
//            cout << "get relative index " << lattice_ref.get_site_by_id(ss).get_relative_index().get_str() << endl;
        }
    }
}

double BondPercolation_v13::occupation_prob() {
#ifdef UNIT_TEST
    if(occupied_site_count > _length*_length){
        cout << "this many sites cannot be occupied : " << __FILE__ << ": " << __LINE__ << endl;
        cout << "occupied_site_count " << occupied_site_count << endl;
        exit(-1);
    }
#endif
    return double(occupied_bond_count) / lattice_ref.get_bond_count();
}

int BondPercolation_v13::merge_clusters(std::vector<int> &bond_neighbors) {
    auto bond_gids = get_bond_gids(bond_neighbors);
//# print("merging clusters ", bond_gids)
    int ref_sz = 0,    root_clstr = 0;
    for (auto bb : bond_gids) {
        int sz = cluster_pool_ref.get_cluster_bond_count(bb);
        if (sz >= ref_sz) {
            root_clstr = bb;
            ref_sz = sz;
        }
    }
    cout << "root cluster is " << root_clstr << endl;
    for (auto bb : bond_gids) {
        if (bb == root_clstr) {
        //# print("bb ", bb, " is a root cluster")
            continue;
        }
    cout << "merging " << bb << " to " << root_clstr << endl;
        cluster_pool_ref.merge_cluster_with(root_clstr, bb, lattice_ref);
    }

    return root_clstr;
}

/**
 * merging with relabeling relative indices
 * @param bond_neighbor_ids
 * @return
 */
int BondPercolation_v13::merge_clusters_v2(std::vector<int> &bond_neighbor_ids) {

    auto bond_neighbors = uniqe_gid_bond_neighbors(bond_neighbor_ids);
    auto bond_gids = get_bond_gids(bond_neighbors);
//# site_gids = self.get_site_gids(site_neighbors)
//# print("site_gids ", site_gids)
//    cout << "bond_neighbors {";
//    print_vectors(bond_neighbors, "}\n");
//    cout << "bond_gids {";
//    print_vectors(bond_gids, "}\n");

//# print("set minus ", set(site_gids) - set(bond_gids))
//# print("merging clusters ", bond_gids)
    int root_clstr = bond_gids[0];
    value_type ref_sz = 0;
    for (auto bbg : bond_gids) {
        if(after_wrapping && bbg == wrapping_cluster_id){
            root_clstr = wrapping_cluster_id;
            break;
        }
        value_type sz = cluster_pool_ref.get_cluster_bond_count(bbg);
        if (sz >= ref_sz) {
            root_clstr = bbg;
            ref_sz = sz;
        }
    }
#ifdef UNIT_TEST
    if(after_wrapping){
        for(auto bbg: bond_gids){
            if (bbg == wrapping_cluster_id && root_clstr != wrapping_cluster_id){
                cout << "Error. Wrapping cluster is connected to one of the neighbors but it's not the root cluster" << endl;

                exit(-1);
            }
        }
    }

    set<int> tmp(bond_gids.begin(), bond_gids.end());
    if(tmp.size() != bond_gids.size()){
        cout << "Error. duplicate gids " << __LINE__ << endl;
        exit(-1);
    }

#endif

//    cout << "root cluster is " << root_clstr << endl;
//    cout << "Assign and relabel currently selected site" << endl;
    for (auto bb : bond_neighbors) {
        int bbg = lattice_ref.get_bond_by_id(bb).get_gid();
        if (bbg == root_clstr) {
//            cout << "bbg " << bbg << " is a root cluster" << endl;
//# relabel and assign the current site here
            lattice_ref.set_site_gid_by_id(selected_id, root_clstr);
            cluster_pool_ref.add_sites(root_clstr, {selected_id});
//# relabeling current site. relative index
//            auto neighbor_site = get_neighbor_site(current_site.get_id(), bb);
            auto neighbor_site = get_neighbor_site(get_current_bond().get_id(), bb);
//            cout << "central site " << selected_id << " neighbor site " << neighbor_site << endl;
            if (lattice_ref.get_site_gid_by_id(neighbor_site) >= 0) {
//# relabel selected site with respect to neighbor site. so neighbor_site is the central site
                auto rri = get_relative_index(neighbor_site, selected_id);
//# rri = self.get_relative_index(self.selected_id, neighbor_site)

//# sitttte = self.lattice_ref.get_site_by_id(self.selected_id)
//                cout << "relative index before " << neighbor_site << endl;
//                cout << selected_id << " => rri " << rri.get_str() << endl;

                lattice_ref.set_relative_index(selected_id, rri);
//                current_site.set_relative_index(rri); // local variable. needs to be updated
            }else{
//                cout << "Does not belong to any cluster yet" << endl;
            }
        }
    }


//    cout << "Relabel all cluster according to root cluster and selected site" << endl;
    int i = 0;
    for (auto bb : bond_neighbors) {
//        cout << "relabeling attempt [" << i << "]" << endl;
        i+= 1;
        int bbg = lattice_ref.get_bond_by_id(bb).get_gid();
        if (bbg == root_clstr) {
//            cout << "bb " << bbg << " is a root cluster" << endl;
            continue;
        }
//        cout << "relabeling relative index of cluster " << bbg << endl;
        relabel_relative_indices(bb);
//        cout << "merging " << bbg << " to " << root_clstr << endl;
        cluster_pool_ref.merge_cluster_with(root_clstr, bbg, lattice_ref);
    }
//    cout << "*************************************DONE relabeling" << endl;

    for (int bbg : bond_gids) {
        if (bbg == root_clstr) {
//            cout << "bb " << bbg << " is a root cluster" << endl;
            continue;
        }
        cluster_pool_ref.clear_cluster(bbg);
    }
    return root_clstr;

}


/**
 * merging with relabeling relative indices
 * @param bond_neighbor_ids
 * @return
 */
int BondPercolation_v13::merge_clusters_v3(std::vector<int> &bond_neighbor_ids) {

    auto bond_neighbors = uniqe_gid_bond_neighbors(bond_neighbor_ids);
    auto bond_gids = get_bond_gids(bond_neighbors);
//# site_gids = self.get_site_gids(site_neighbors)
//# print("site_gids ", site_gids)
//    cout << "bond_neighbors {";
//    print_vectors(bond_neighbors, "}\n");
//    cout << "bond_gids {";
//    print_vectors(bond_gids, "}\n");

//# print("set minus ", set(site_gids) - set(bond_gids))
//# print("merging clusters ", bond_gids)
    int root_clstr = find_root_cluster(bond_gids);
#ifdef UNIT_TEST
    if(after_wrapping){
        for(auto bbg: bond_gids){
            if (bbg == wrapping_cluster_id && root_clstr != wrapping_cluster_id){
                cout << "Error. Wrapping cluster is connected to one of the neighbors but it's not the root cluster" << endl;

                exit(-1);
            }
        }
    }

    set<int> tmp(bond_gids.begin(), bond_gids.end());
    if(tmp.size() != bond_gids.size()){
        cout << "Error. duplicate gids " << __LINE__ << endl;
        exit(-1);
    }

#endif

//    cout << "root cluster is " << root_clstr << endl;
//    cout << "Assign and relabel currently selected site" << endl;
    for (auto bb : bond_neighbors) {
        int bbg = lattice_ref.get_bond_by_id(bb).get_gid();
        if (bbg == root_clstr) {

//# relabel and assign the current site here
            lattice_ref.set_site_gid_by_id(selected_id, root_clstr);
            cluster_pool_ref.add_sites(root_clstr, {selected_id});
//# relabeling current site. relative index
//            auto neighbor_site = get_neighbor_site(current_site.get_id(), bb);
            auto neighbor_site = get_neighbor_site(get_current_bond().get_id(), bb);
//            cout << "central site " << selected_id << " neighbor site " << neighbor_site << endl;
            if (lattice_ref.get_site_gid_by_id(neighbor_site) >= 0) {
//# relabel selected site with respect to neighbor site. so neighbor_site is the central site
                auto rri = get_relative_index(neighbor_site, selected_id);
//# rri = self.get_relative_index(self.selected_id, neighbor_site)

//# sitttte = self.lattice_ref.get_site_by_id(self.selected_id)
//                cout << "relative index before " << neighbor_site << endl;
//                cout << selected_id << " => rri " << rri.get_str() << endl;

                lattice_ref.set_relative_index(selected_id, rri);
//                current_site.set_relative_index(rri); // local variable. needs to be updated
            }else{
//                cout << "Does not belong to any cluster yet" << endl;
            }
        }
    }


//    cout << "Relabel all cluster according to root cluster and selected site" << endl;
    int i = 0;
    for (auto bb : bond_neighbors) {
//        cout << "relabeling attempt [" << i << "]" << endl;
        i+= 1;
        int bbg = lattice_ref.get_bond_by_id(bb).get_gid();
        if (bbg == root_clstr) {
//            cout << "bb " << bbg << " is a root cluster" << endl;
            continue;
        }
//        cout << "relabeling relative index of cluster " << bbg << endl;
        relabel_relative_indices(bb);
//        cout << "merging " << bbg << " to " << root_clstr << endl;
        cluster_pool_ref.merge_cluster_with(root_clstr, bbg, lattice_ref);
    }
//    cout << "*************************************DONE relabeling" << endl;

    for (int bbg : bond_gids) {
        if (bbg == root_clstr) {
//            cout << "bb " << bbg << " is a root cluster" << endl;
            continue;
        }
        cluster_pool_ref.clear_cluster(bbg);
    }
    return root_clstr;

}

/**
 * merging with relabeling relative indices
 * @param bond_neighbor_ids
 * @return
 */
int BondPercolation_v13::merge_clusters_v4(std::vector<int> &bond_neighbor_ids) {

    int rootcounter = 0;
    auto bond_neighbors = uniqe_gid_bond_neighbors(bond_neighbor_ids);
    auto bond_gids = get_bond_gids(bond_neighbors);
//# site_gids = self.get_site_gids(site_neighbors)
//# print("site_gids ", site_gids)
//    cout << "bond_neighbors {";
//    print_vectors(bond_neighbors, "}\n");
//    cout << "bond_gids {";
//    print_vectors(bond_gids, "}\n");

//# print("set minus ", set(site_gids) - set(bond_gids))
//# print("merging clusters ", bond_gids)
//    int root_clstr = find_root_cluster(bond_gids);
    int suitable_bond_id = find_bonds_connected_to_root_cluster(bond_neighbors);
    int root_clstr = lattice_ref.get_bond_by_id(suitable_bond_id).get_gid();

#ifdef UNIT_TEST
    if(after_wrapping){
        for(auto bbg: bond_gids){
            if (bbg == wrapping_cluster_id && root_clstr != wrapping_cluster_id){
                cout << "Error. Wrapping cluster is connected to one of the neighbors but it's not the root cluster" << endl;

                exit(-1);
            }
        }
    }

    set<int> tmp(bond_gids.begin(), bond_gids.end());
    if(tmp.size() != bond_gids.size()){
        cout << "Error. duplicate gids " << __LINE__ << endl;
        exit(-1);
    }

#endif

//    cout << "root cluster is " << root_clstr << endl;
//    cout << "Assign and relabel currently selected site" << endl;
//# relabel and assign the current site here
    lattice_ref.set_site_gid_by_id(selected_id, root_clstr);
    cluster_pool_ref.add_sites(root_clstr, {selected_id});
//# relabeling current site. relative index
//            auto neighbor_site = get_neighbor_site(current_site.get_id(), bb);
    auto neighbor_site = get_neighbor_site(selected_id, suitable_bond_id);
//            cout << "central site " << selected_id << " neighbor site " << neighbor_site << endl;
    if (lattice_ref.get_site_gid_by_id(neighbor_site) >= 0) {
//# relabel selected site with respect to neighbor site. so neighbor_site is the central site
        auto rri = get_relative_index(neighbor_site, selected_id);
//# rri = self.get_relative_index(self.selected_id, neighbor_site)

//# sitttte = self.lattice_ref.get_site_by_id(self.selected_id)
//                cout << "relative index before " << neighbor_site << endl;
//                cout << selected_id << " => rri " << rri.get_str() << endl;

        lattice_ref.set_relative_index(selected_id, rri);
//                current_site.set_relative_index(rri); // local variable. needs to be updated
    }else{
//                cout << "Does not belong to any cluster yet" << endl;
    }


//    cout << "Relabel all cluster according to root cluster and selected site" << endl;
    int i = 0;
    for (auto bb : bond_neighbors) {
//        cout << "relabeling attempt [" << i << "]" << endl;
        i+= 1;
        int bbg = lattice_ref.get_bond_by_id(bb).get_gid();
        if (bbg == root_clstr) {
//            cout << "bb " << bbg << " is a root cluster" << endl;
            continue;
        }
//        cout << "relabeling relative index of cluster " << bbg << endl;
        relabel_relative_indices(bb);
//        cout << "merging " << bbg << " to " << root_clstr << endl;
        cluster_pool_ref.merge_cluster_with(root_clstr, bbg, lattice_ref);
    }
//    cout << "*************************************DONE relabeling" << endl;


#ifdef UNIT_TEST
// check new cluster for errors
    for (auto bb : bond_neighbors) {

        int bbg = lattice_ref.get_bond_by_id(bb).get_gid();
        if (bbg != root_clstr) {
            // these clsuters must be empty
            value_type bonds_c = cluster_pool_ref.get_cluster_bond_count(bbg);
            value_type sites_c = cluster_pool_ref.get_cluster_site_count(bbg);
            if(bonds_c > 0 || sites_c > 0){
                cerr << "Error " << " Cluster wasn't cleared properly" << endl;
                exit(-1);
            }
        }else{
            // root cluster
            // check for gid of bonds and sites. and relative indices of the sites
            auto bonds = cluster_pool_ref.get_bonds(bbg);
            auto sites = cluster_pool_ref.get_sites(bbg);
            for(auto bb: bonds){
                auto gid = lattice_ref.get_bond_gid_by_id(bb);
                if(gid != root_clstr){
                    cerr << "Error. gid of bonds mismatch" << endl;
                    exit(-1);
                }
            }

            for(auto ss: sites){
                auto gid = lattice_ref.get_site_gid_by_id(ss);
                if(gid != root_clstr){
                    cerr << "Error. gid of sites mismatch" << endl;
                    exit(-1);
                }
            }
        }

    }
#endif

    return root_clstr;

}

/**
 * Find suitable id/index for the root cluster among the provided gids.
 * If the critical point is reached then the wrapping cluster is the root clsuter if it's present in `bond_gids`
 * Otherwise the largest one is.
 * @param bond_gids
 * @return
 */
int BondPercolation_v13::find_root_cluster(const vector<int> &bond_gids)  {
    int root_clstr = bond_gids[0];
    value_type ref_sz = 0;
    for (auto bbg : bond_gids) {
        if(after_wrapping && bbg == wrapping_cluster_id){
            root_clstr = wrapping_cluster_id;
            break;
        }
        value_type sz = cluster_pool_ref.get_cluster_bond_count(bbg);
        if (sz >= ref_sz) {
            root_clstr = bbg;
            ref_sz = sz;
        }
    }
    return root_clstr;
}

/**
 * bonds of bond_ids may be connected to different cluster. we need one of them that is connected to
 * he largest cluster or wrapping cluster.
 * wrapping cluster is preferred if after wrapping.
 * @param bond_ids
 * @return
 */
int BondPercolation_v13::find_bonds_connected_to_root_cluster(const vector<int> &bond_ids)  {
    int selected_id = bond_ids[0];
    int bbg = 0;

    value_type ref_sz = 0;
    for (auto bb : bond_ids) {
        bbg = lattice_ref.get_bond_by_id(bb).get_gid();
        if(after_wrapping && bbg == wrapping_cluster_id){
//            selected_id = wrapping_cluster_id;
            selected_id  = bb;
            break;
        }
        value_type sz = cluster_pool_ref.get_cluster_bond_count(bbg);
        if (sz >= ref_sz) {
            selected_id = bb;
            ref_sz = sz;
        }
    }
    return selected_id;
}

/**
 * @brief view site id list in curly brace {}.
 *  The unoccupied sites are inside square bracket [].
 * The site with angle bracket is the currenlty selected site <>.
 * 
 */
void BondPercolation_v13::view_site_id_list(){
    std::cout << "view_site_id_list() {";
    auto got_selected_site_flag = false;
    for(int i=0; i< bond_ids_indices.size(); ++i){
        if (i == current_idx){
            if(got_selected_site_flag){
                // 0 for background Color(Black)
                // 4 for text color(Red)
                system("Color 04");
                std::cerr << "How can the current site be outside the unoccupied sites?" << std::endl;
                std::cout << "Probably you are calling view_site_id_list without occupying any new sites" << std::endl;
                // exit(0);
            }
            std::cout << "[";
        }
        auto a = bond_ids_indices[i];
        if(a == selected_id){
            std::cout <<"<" << a << ">,";
            got_selected_site_flag = true;
            continue;
        }
        std::cout << a << ",";
    }
    // for(auto a: bond_ids_indices){
    //     if(a == selected_id){
    //         std::cout <<"<" << a << ">,";
    //         continue;
    //     }
    //     std::cout << a << ",";
    // }
    std::cout << "]}" << std::endl;
}


/**
 * @brief 
 * 
 * cluster size is measured by number of bonds in it.
 * @return const std::vector<double> 
 */
const std::vector<double> BondPercolation_v13::clusterSizeDistribution(){
    vector<double> cluster_counts;
    size_t n, sz;
    auto cluster_count = cluster_pool_ref.cluster_count_v2();
    // if (cluster_count != cluster_pool_ref.cluster_count()){
    //    cout << "number of clusters " << cluster_pool_ref.cluster_count() << " vs " << cluster_count << endl;
    //    exit(0);
    // }
    size_t mx_n{}, total_bonds{};
    
    for(size_t i{}; i < cluster_count; ++i){
//         if(_clusters[i].empty())   continue;
// //        cout << i << " th cluster ";
        
        n = cluster_pool_ref.get_cluster_bond_count(i);
        total_bonds += n;
        if(n > mx_n){
            mx_n = n;
        }
        sz = cluster_counts.size();
        if (n >= sz){
            cluster_counts.resize(n+1);
        }
//        cout << cluster_counts[n];
        ++cluster_counts[n]; // increase cluster count

//        cout << endl;
    }
    // since isolated bonds forms cluster of unit size and it might not be in the cluster
    cluster_counts[1] +=  lattice_ref.get_bond_count() - total_bonds;
//    cout << "before returning " << mx_n << endl;
    return cluster_counts;
}


/**
 * Return bonds ids with unique gids.
 * @param bond_neighbors_ids : bond ids.
 * @return bonds ids
 */
std::vector<int> BondPercolation_v13::uniqe_gid_bond_neighbors(std::vector<int> &bond_neighbors_ids) {
//    cout << "bond_neighbors {";
//    print_vectors(bond_neighbors, "}\n");
    vector<int> gids, unique_bond_ids;
    for (int bb : bond_neighbors_ids) {
        int bbg = lattice_ref.get_bond_by_id(bb).get_gid();
        if (find_elm(gids, bbg))  continue;
        else {
            gids.push_back(bbg);
            unique_bond_ids.push_back(bb);
        }
    }
    return unique_bond_ids;

}

/**
 * @brief Detect Wrapping
 * 
 * Remember to call this method whn at least one site is occupied
 * @return true 
 * @return false 
 */
bool BondPercolation_v13::detect_wrapping() {
    if (after_wrapping) return true;

//# print("detect_wrapping")
    auto neighbors = lattice_ref.get_sites_for_wrapping_test(selected_id);
//# print("neighbors of self.selected_id with same gid : ", neighbors)
//    auto central_r_index = current_site.get_relative_index();
    // auto central_r_index = get_current_bond().get_relative_index();
    cout << "TODO  : there are no relative index of a bond. find two sites and then find relative index" 
    << __FILE__ << " " << __LINE__ << endl;
    for (auto ss : neighbors) {
        auto rss = lattice_ref.get_site_by_id(ss).get_relative_index();
        // auto delta_x = central_r_index.x_coord() - rss.x_coord();
        // auto delta_y = central_r_index.y_coord() - rss.y_coord();
//         if ((abs(delta_x) > 1) or (abs(delta_y) > 1)) {
//         //    cout << selected_id << " and " << ss << " are connected via wrapping" << endl;
// //            cout << "indices are " << lattice_ref.get_site_by_id(selected_id).get_index().get_str() <<
// //           " and " << lattice_ref.get_site_by_id(ss).get_index().get_str() << endl;
// //            cout << "relative " << central_r_index.get_str() << " - " << rss.get_str() << endl;

//             after_wrapping = true;
//             wrapping_cluster_id = lattice_ref.get_site_by_id(selected_id).get_gid();
//             bond_count_pc = size_t(occupied_bond_count);
//             wrapping_cluster_site_count_pc = cluster_pool_ref.get_cluster_site_count(wrapping_cluster_id);
//             wrapping_cluster_bond_count_pc = cluster_pool_ref.get_cluster_bond_count(wrapping_cluster_id);

//             return true;
//         }
    }
    // cout << "SitePercolation_v13::detect_wrapping() | Got to Line " << __LINE__ << endl;
    return false;
}

Bond_v13 &BondPercolation_v13::get_current_bond() {
    return lattice_ref.get_bond_by_id(selected_id);
}





