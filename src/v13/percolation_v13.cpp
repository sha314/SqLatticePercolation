//
// Created by shahnoor on 2/15/21.
//

#include <cmath>
#include <set>
#include <random>
#include "percolation_v13.h"

using namespace std;

void Percolation_v13::viewLattice(int formatt) {
        if (formatt < 0) {
            cout << "undefined format" << endl;
        }
        else if(formatt == 3) {
            lattice_ref.view_relative_index();
        }
        else if(formatt == 4) {
            lattice_ref.view_site_gids();
        }
        else {
            lattice_ref.view(formatt);
        }
}

RelativeIndex_v13 Percolation_v13::get_change_in_relative_index(RelativeIndex_v13 old_relative_index,
                                                                RelativeIndex_v13 new_relative_index) {
//# print("get_change_in_relative_index")
//# print("old_relative_index - new_relative_index = ", old_relative_index, " - ", new_relative_index)
    auto change = new_relative_index - old_relative_index;
//# print("get_change_in_relative_index. type of change ", type(change))
//# print("change ", change)
    return RelativeIndex_v13(change);
}

/**
 * Get the gid of the bonds given in the argument.
 * @param bond_ids : ids of bonds
 * @return gids of bonds
 */
vector<int> Percolation_v13::get_bond_gids(std::vector<int> &bond_ids) {
//    set<int> gids;
    vector<int> gids;
    for (auto bb : bond_ids) {
        int gid = lattice_ref.get_bond_by_id(bb).get_gid();
//        gids.insert(gid);
        gids.push_back(gid);
    }
//    vector<int> vec(gids.begin(), gids.end());
//    return vec;
    return gids;
}

vector<int> Percolation_v13::get_site_gids(std::vector<int> &site_ids) {
    set<int> gids;
    for (auto ss : site_ids) {
        int gid = lattice_ref.get_site_by_id(ss).get_gid();
//        cout << "site " << ss << " gid => " << gid << endl;
        gids.insert(gid);
    }
    vector<int> vec(gids.begin(), gids.end());
    return vec;
}

/**
 * when delta X of relative indices are greater than 1 and there are in
    the opposite edge of the lattice.

 * @param delta_X
 * @return
 */
RelativeIndex_v13 Percolation_v13::wrapping_correction_relative_index(RelativeIndex_v13 delta_X) {
//# LL = self.lattice_ref.length
    int xx = delta_X.x_coord();
    int yy = delta_X.y_coord();
    if (abs(xx) > 1) {
        xx = -xx / abs(xx);
    }
    if (abs(yy) > 1) {
        yy = -yy / abs(yy);
    }
//# print(type(xx), " and ", type(yy))
    return RelativeIndex_v13(xx, yy);
}

Index_v13 Percolation_v13::wrapping_correction_relative_index(Index_v13 delta_X) {
//# LL = self.lattice_ref.length
    int xx = delta_X.row();
    int yy = delta_X.col();
    if (abs(xx) > 1) {
        xx = -xx / abs(xx);
    }
    if (abs(yy) > 1) {
        yy = -yy / abs(yy);
    }
//# print(type(xx), " and ", type(yy))
    return Index_v13(xx, yy);
}

/**
 * neighbor_site_id will get a new relative index based on central_site_id. only one condition,
                new site must be a neighbor of old site

 * @param central_site_id
 * @param neighbor_site_id
 * @return
 */
RelativeIndex_v13 Percolation_v13::get_relative_index(int central_site_id, int neighbor_site_id) {
//    cout << "central_site_id " << central_site_id << endl;
//    cout << "neighbor_site_id " << neighbor_site_id << endl;
    Index_v13 central_index = lattice_ref.get_site_by_id(central_site_id).get_index();
    Index_v13 neighbor_index = lattice_ref.get_site_by_id(neighbor_site_id).get_index();
//    cout << "central_index " << central_index.get_str() << endl;
//    cout << "neighbor_index " << neighbor_index.get_str() << endl;
    Index_v13 idx = neighbor_index.subtract(central_index);
//    cout << "idx " << idx.get_str() << endl;
    idx = wrapping_correction_relative_index(idx);
//    cout << "after wrapping correction " << idx.get_str() << endl;
    Index_v13 old_relative_index = lattice_ref.get_site_by_id(central_site_id).get_relative_index();
//# new_relative_index = self.lattice_ref.get_site_by_id(new_site_id).get_relative_index()
//# print(old_relative_index, " old_relative_index type ", type(old_relative_index))
    auto new_relative_index = old_relative_index + idx;
//# print("new_relative_index type ", type(new_relative_index))
//# print("new_relative_index type ", type(RelativeIndex(index=new_relative_index)))
    return RelativeIndex_v13(new_relative_index);
}

void Percolation_v13::setRandomState(size_t seed, bool generate_seed) {
//    size_t seed = 0;
//    cerr << "automatic seeding is commented : line " << __LINE__ << endl;
    _random_state = seed;
    if(generate_seed) {
        std::random_device _rd;
        _random_state = _rd();
    }else{
        cerr << "generate_seed = false : line " << __LINE__ << endl;
    }
    _random_engine.seed(_random_state); // seeding
    cout << "seeding with " << _random_state << endl;
}

SitePercolation_v13::SitePercolation_v13(int length, value_type seed) : Percolation_v13(length, seed) {
    cout << "class : SitePercolation_v13" << endl;
    init_clusters();

    site_ids_indices = lattice_ref.get_site_id_list();
//    cout << "{";
//    for(auto s: site_ids_indices){
//        cout << s << ",";
//    }
//    cout << "}" << endl;
    int current_idx = 0;
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

void SitePercolation_v13::init_clusters() {
//    cluster_pool_ref.reset();
    auto bonds = lattice_ref.get_bond_id_list();
    for (auto bb: bonds){
        cluster_pool_ref.create_new_cluster(-1, bb, lattice_ref);
    }
}

void SitePercolation_v13::reset() {
    Percolation_v13::reset();

    site_ids_indices = lattice_ref.get_site_id_list();;
    current_idx = 0;
    occupied_site_count = 0;
    site_count_pc = 0;
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
}

double SitePercolation_v13::entropy_v1() {
    double H = 0;
    int empty_count = 0;
    for(int i=0; i < cluster_count;++i){
        double b_count = cluster_pool_ref.get_cluster_bond_count(i);
        if (b_count==0) {
            empty_count += 1;
            continue;
        }
        double mu = b_count / lattice_ref.get_bond_count();
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

double SitePercolation_v13::order_param_wrapping() {
    if (after_wrapping) {
//# print("wrapping cluster id ", self.wrapping_cluster_id)
        double count = cluster_pool_ref.get_cluster_bond_count(wrapping_cluster_id);
#ifdef UNIT_TEST
    if(count <= 0){
        cout << "Error . After wrapping, number of bond in the wrapping cluster cannot be zero " << count << endl;
        cout << "wrapping_cluster_id " << wrapping_cluster_id << endl;
        cluster_pool_ref.view(0);
        exit(-1);
    }
#endif
        return count / lattice_ref.get_bond_count();
    }
    return 0;
}

int SitePercolation_v13::get_neighbor_site(int central_site_id, int connecting_bond_id) {
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

std::vector<int> SitePercolation_v13::get_connected_sites(Site_v13 site, std::vector<int> &bond_neighbors) {
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

P_STATUS SitePercolation_v13::select_site() {

    if (current_idx >= lattice_ref.get_site_count()) {
//        cout << "No sites to occupy" << endl;
        return P_STATUS::EMPTY_SITE_LIST;
    }
    selected_id = site_ids_indices[current_idx];
    current_idx += 1;
//    current_site = lattice_ref.get_site_by_id(selected_id);
//    cout << "selected id " << selected_id << " site " << current_site.get_str() << endl;
//    cout << "selected id " << selected_id << " site " << get_current_site().get_str() << endl;
    return P_STATUS::SUCESS;

}

bool SitePercolation_v13::place_one_site() {
//    cout << "************************ place_one_site. count " << current_idx << endl;
    status = select_site();
    if(status == P_STATUS::SUCESS) {

//        cout << "selected site ", self.current_site.get_index(), " id ", self.current_site.get_id())
        lattice_ref.init_relative_index(selected_id);  // initialize        relative index
//        auto bond_neighbors = current_site.get_connecting_bonds();

        auto bond_neighbors = get_current_site().get_connecting_bonds();
//# site_neighbors = self.get_connected_sites(self.current_site, bond_neighbors)

        entropy_subtract(bond_neighbors);

        auto merged_cluster_index = merge_clusters_v4(bond_neighbors);

        track_largest_cluster(merged_cluster_index);
        entropy_add(merged_cluster_index);

//# self.lattice_ref.set_site_gid_by_id(selected_id, merged_cluster_index)
//# self.cluster_pool_ref.add_sites(merged_cluster_index, selected_id)
        occupied_site_count += 1;

#ifdef UNIT_TEST
        if ( bond_neighbors.size() != 4){
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

void SitePercolation_v13::track_largest_cluster(int new_cluster) {
    auto new_size = cluster_pool_ref.get_cluster_bond_count(new_cluster);
//# self.cluster_pool_ref.get_cluster_site_count(new_cluster)
    if (new_size > largest_cluster_sz) {
        largest_cluster_id = new_cluster;
        largest_cluster_sz = new_size;
    }
}

void SitePercolation_v13::entropy_subtract(std::vector<int> &bond_neighbors) {
//# print("entropy_subtract")
    auto bonds = lattice_ref.get_neighbor_bonds(selected_id);
//# print(self.current_site, " neighbors => ", sites)
    set<int> gids;
    for(auto bb: bonds) {
        auto gid = lattice_ref.get_bond_gid_by_id(bb);
        if (gid == -1)        continue;
        gids.insert(gid);
    }
//# print("gids ", gids)
    double H = 0, mu;
    double bc = lattice_ref.get_bond_count();
    for (auto gg : gids) {
        auto b_count = cluster_pool_ref.get_cluster_bond_count(gg);
        if (b_count == 0)        continue;

        mu = b_count / bc;
        H += mu * log(mu);
    }
//# print("before ", self.entropy_value)
    entropy_value += H;
//# print("after ", self.entropy_value)

}

void SitePercolation_v13::entropy_add(int new_cluster_id) {
//# print("entropy_add")
    double b_count = cluster_pool_ref.get_cluster_bond_count(new_cluster_id);
    double mu = b_count / lattice_ref.get_bond_count();
//# print("before ", self.entropy_value)
    entropy_value -= mu*log(mu);
//# print("after ", self.entropy_value)
}

double SitePercolation_v13::entropy() {
//    return entropy_v1();
    return entropy_v2();
}

double SitePercolation_v13::entropy_v2() {
    return entropy_value;
}

/**
 *
 * @param connecting_bond_id
 */
void SitePercolation_v13::relabel_relative_indices(int connecting_bond_id) {
//    cout << "SitePercolation_v13::relabel_relative_indices" << endl;
    auto bond = lattice_ref.get_bond_by_id(connecting_bond_id);
    auto bbg = bond.get_gid();
//    auto central_site = current_site.get_id();
    auto central_site = get_current_site().get_id();
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

double SitePercolation_v13::occupation_prob() {
#ifdef UNIT_TEST
    if(occupied_site_count > _length*_length){
        cout << "this many sites cannot be occupied : " << __FILE__ << ": " << __LINE__ << endl;
        cout << "occupied_site_count " << occupied_site_count << endl;
        exit(-1);
    }
#endif
    return double(occupied_site_count) / lattice_ref.get_site_count();
}

int SitePercolation_v13::merge_clusters(std::vector<int> &bond_neighbors) {
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
int SitePercolation_v13::merge_clusters_v2(std::vector<int> &bond_neighbor_ids) {

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
            auto neighbor_site = get_neighbor_site(get_current_site().get_id(), bb);
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
int SitePercolation_v13::merge_clusters_v3(std::vector<int> &bond_neighbor_ids) {

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
            auto neighbor_site = get_neighbor_site(get_current_site().get_id(), bb);
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
int SitePercolation_v13::merge_clusters_v4(std::vector<int> &bond_neighbor_ids) {

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
int SitePercolation_v13::find_root_cluster(const vector<int> &bond_gids)  {
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
int SitePercolation_v13::find_bonds_connected_to_root_cluster(const vector<int> &bond_ids)  {
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
 * Return bonds ids with unique gids.
 * @param bond_neighbors_ids : bond ids.
 * @return bonds ids
 */
std::vector<int> SitePercolation_v13::uniqe_gid_bond_neighbors(std::vector<int> &bond_neighbors_ids) {
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

bool SitePercolation_v13::detect_wrapping() {
    if (after_wrapping) return true;

//# print("detect_wrapping")
    auto neighbors = lattice_ref.get_sites_for_wrapping_test(selected_id);
//# print("neighbors of self.selected_id with same gid : ", neighbors)
//    auto central_r_index = current_site.get_relative_index();
    auto central_r_index = get_current_site().get_relative_index();
    for (auto ss : neighbors) {
        auto rss = lattice_ref.get_site_by_id(ss).get_relative_index();
        auto delta_x = central_r_index.x_coord() - rss.x_coord();
        auto delta_y = central_r_index.y_coord() - rss.y_coord();
        if ((abs(delta_x) > 1) or (abs(delta_y) > 1)) {
//            cout << selected_id << " and " << ss << " are connected via wrapping" << endl;
//            cout << "indices are " << lattice_ref.get_site_by_id(selected_id).get_index().get_str() <<
//           " and " << lattice_ref.get_site_by_id(ss).get_index().get_str() << endl;
//            cout << "relative " << central_r_index.get_str() << " - " << rss.get_str() << endl;

            after_wrapping = true;
            wrapping_cluster_id = lattice_ref.get_site_by_id(selected_id).get_gid();
            site_count_pc = size_t(occupied_site_count);
            wrapping_cluster_site_count_pc = cluster_pool_ref.get_cluster_site_count(wrapping_cluster_id);
            wrapping_cluster_bond_count_pc = cluster_pool_ref.get_cluster_bond_count(wrapping_cluster_id);

            return true;
        }
    }
    return false;
}

Site_v13 &SitePercolation_v13::get_current_site() {
    return lattice_ref.get_site_by_id(selected_id);
}

void SitePercolationL0_v13::reset() {
    SitePercolation_v13::reset();
//    int ll = get_length();
//    int l_squared = ll*ll;

    if (first_run) {
        occupation_prob_list.clear();
//        occupation_prob_list.resize(l_squared);
    }

    entropy_list.clear();
    order_wrapping_list.clear();
    order_largest_list.clear();


//    entropy_list.resize(l_squared);
//    order_wrapping_list.resize(l_squared);
//    order_largest_list.resize(l_squared);
}

SitePercolationL0_v13::SitePercolationL0_v13(int length, value_type seed, bool generate_seed) : SitePercolation_v13(length, seed) {
    cout << "class : SitePercolationL0_v13" << endl;
    first_run = true;

    int ll = get_length();
    int l_squared = ll*ll;
//    occupation_prob_list.resize(l_squared);
//    entropy_list.resize(l_squared);
//    order_wrapping_list.resize(l_squared);
//    order_largest_list.resize(l_squared);
    setRandomState(seed, generate_seed);
}

void SitePercolationL0_v13::run_once() {
//# sq_lattice_p.viewLattice(3)
//# sq_lattice_p.viewCluster()
    double p, H, P1, P2;

    while (place_one_site()) {
        detect_wrapping();
        if (first_run) {
            p = occupation_prob();
            occupation_prob_list.push_back(p);
        }
        H = entropy();
        P1 = order_param_wrapping();
        P2 = order_param_largest_clstr();

        entropy_list.push_back(H);
        order_wrapping_list.push_back(P1);
        order_largest_list.push_back(P2);
#ifdef UNIT_TEST
        double  H1 = entropy_v1();
        double  H2 = entropy_v2();
        if(abs(H1 - H2) > 1e-6){
            cout << "Error : Entropy v1 and v2 are not equal : " << __LINE__ << endl;
            cout << "H1 = " << H1 << endl;
            cout << "H2 = " << H2 << endl;
            cout << "max_entropy  " << max_entropy << endl;
            exit(-1);
        }
#endif
    }

#ifdef UNIT_TEST
//    P1 = order_param_wrapping();
//    P2 = order_param_largest_clstr();

    if (abs(P1-1.0) > 1e-6){
        cout << "Error : order parameter wrapping P1 not equal to 1.0. line " << __LINE__ << endl;
        cout << "P1 = " << P1 << endl;
        exit(-1);
    }

    if (abs(P2-1.0) > 1e-6){
        cout << "Error : order parameter largest P2 not equal to 1.0. line " << __LINE__ << endl;
        cout << "P2 = " << P2 << endl;
        exit(-1);
    }
    p = occupation_prob();

    if (abs(p-1.0) > 1e-6){
        cout << "Error : occupation_prob p not equal to 1.0. line " << __LINE__ << endl;
        exit(-1);
    }

    cout << "entropy_list {" << endl;
    for(auto hh: entropy_list){
        cout << hh << endl;
    }
    cout << "              }" << endl;
#endif

    first_run = false;
}

/**
 * make use of STATUS. Since for L1 and L2 we can have place_one_site() returning true but no site was allocated.
 */
void SitePercolationL0_v13::run_once_v2() {
//# sq_lattice_p.viewLattice(3)
//# sq_lattice_p.viewCluster()
    double p, H, P1, P2;

    while (place_one_site()) {
        if (status != P_STATUS::SUCESS) continue;
        detect_wrapping();
        if (first_run) {
            p = occupation_prob();
            occupation_prob_list.push_back(p);
        }
        H = entropy();
        P1 = order_param_wrapping();
        P2 = order_param_largest_clstr();

        entropy_list.push_back(H);
        order_wrapping_list.push_back(P1);
        order_largest_list.push_back(P2);
#ifdef UNIT_TEST
        double  H1 = entropy_v1();
        double  H2 = entropy_v2();
        if(abs(H1 - H2) > 1e-6){
            cout << "Error : Entropy v1 and v2 are not equal : " << __LINE__ << endl;
            cout << "H1 = " << H1 << endl;
            cout << "H2 = " << H2 << endl;
            cout << "max_entropy  " << max_entropy << endl;
            exit(-1);
        }
#endif
    }

#ifdef UNIT_TEST
//    P1 = order_param_wrapping();
//    P2 = order_param_largest_clstr();

    if (abs(P1-1.0) > 1e-6){
        cout << "Error : order parameter wrapping P1 not equal to 1.0. line " << __LINE__ << endl;
        cout << "P1 = " << P1 << endl;
        exit(-1);
    }

    if (abs(P2-1.0) > 1e-6){
        cout << "Error : order parameter largest P2 not equal to 1.0. line " << __LINE__ << endl;
        cout << "P2 = " << P2 << endl;
        exit(-1);
    }
    p = occupation_prob();

    if (abs(p-1.0) > 1e-6){
        cout << "Error : occupation_prob p not equal to 1.0. line " << __LINE__ << endl;
        exit(-1);
    }

    cout << "entropy_list {" << endl;
    for(auto hh: entropy_list){
        cout << hh << endl;
    }
    cout << "              }" << endl;
#endif

    first_run = false;
}

void SitePercolationL0_v13::test_cluster() {
    size_t count = cluster_pool_ref.cluster_count(true);
    for(size_t i{}; i < count; ++i){
        auto clstr = cluster_pool_ref.get_cluster(i);
        if(clstr.is_empty()){
            if(clstr.get_site_count() != 0 || clstr.get_bond_count() != 0){
                cerr << "Empty cluster contains something : " << __FILE__ << ": " << __LINE__ << endl;
                exit(-1);
            }
        }

        auto gid_cluster = clstr.get_gid();
        auto id_cluster = clstr.get_id();
        if(id_cluster != gid_cluster){
            cout << "id and gid of a cluster does not match. " << __FILE__ << " : " << __LINE__ << endl;
        }
        auto site_ids = clstr.get_sites();
        int root_site_id = -1;
        int root_count = 0;
        for(auto ss: site_ids){
            auto gid = lattice_ref.get_site_gid_by_id(ss);
            if(gid != gid_cluster){
                cerr << "gid mismatch. Site. " << __FILE__ << " : " << __LINE__ << endl;
                exit(-1);
            }
            auto sitesite = lattice_ref.get_site_by_id(ss);
            auto rr = sitesite.get_relative_index();
            if(rr.col() == 0 && rr.row() == 0){
//                cout << "root site found " << endl;
                root_site_id = ss;
                root_count++;
//                cout << "id " << sitesite.get_id() << " index "
//                     << sitesite.get_index().get_str() << " gid " << sitesite.get_gid()
//                     << " relative index " << sitesite.get_relative_index().get_str() << endl;
            }
        }
        if(site_ids.size() > 0) {
            if (root_count > 1) {
                cout << "More than 1 root site. " << __FILE__ << " : " << __LINE__ << endl;
                exit(-1);
            }
            if (root_site_id < 0) {
                cout << "Root site not found. " << __FILE__ << " : " << __LINE__ << endl;
            }
        }
        auto bond_ids = clstr.get_bonds();
        for(auto bb: bond_ids){
            auto gid = lattice_ref.get_bond_gid_by_id(bb);
            if(gid != gid_cluster){
                cerr << "gid mismatch. Bond. " << __FILE__ << " : " << __LINE__ << endl;
                exit(-1);
            }
        }

        // now check the relatice index of all sites given the root site.

    }

}
