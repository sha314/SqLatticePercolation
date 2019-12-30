//
// Created by shahnoor on 12/30/19.
//

#include <algorithm>
#include "percolation_bond_v12.h"

using namespace std;


BondPercolation_pb_v12::BondPercolation_pb_v12(int length) : Percolation_v12(length) {
    auto bonds = _lattice.getBonds();
    for(auto b: bonds) {
        auto id = _lattice.getBondID(b);
        index_sequence.emplace_back(id);
    }
    _max_iteration_limit = size_t(2*length*length);
}

void BondPercolation_pb_v12::init() {
    //
    _entropy = logl(maxSites());

    // activate bonds and initialize cluster
    _clusters.resize(maxSites());
//    cout << "_clusters.size() " << _clusters.size() << endl;
    auto sites = _lattice.getSites();
    for(int i{}; i < _clusters.size(); ++i){

        auto id = _lattice.getSiteID(sites[i]);
        _lattice.activateSite(sites[i]);
        _lattice.setGroupIDSite(sites[i], i);
        _clusters[i].addSite(id);
        _clusters[i].setGroupID(i);
    }

    randomized_index = index_sequence;

    std::shuffle(randomized_index .begin(), randomized_index.end(), _random);

}

bool BondPercolation_pb_v12::occupy() {
    if(index_counter >= randomized_index.size()) return false;
    // select a site randomly
//    size_t i = _random() % randomized_index.size();
    id_last_bond = randomized_index[index_counter];
    index_counter++;

    // activate the site
    _lattice.activateSite(id_last_bond);
    ++_number_of_occupied_sites;
    manageCluster();


    return true;
}

/**
 * Steps :
 *  1. find the largest cluster among bond neighbors of selected site.
 *  2. get largest cluster group id as root.
 *  3. relabel selected site.
 *  4.
 */
void BondPercolation_pb_v12::manageCluster() {// find it's neighbors. sites and then bonds
    auto bond_connected = _lattice.getSite(id_last_bond).connectedBondIDs();
    // find which clusters the bonds belong
    // find a suitble root cluster
    int root{-1};
    size_t root_size{}, tmp{};
    set<int> gids; // to prevent repetition
    for(auto b: bond_connected){
        auto gid = _lattice.getGroupIDBond(b);
        gids.emplace(gid);
        tmp = _clusters[gid].numberOfSites();
        if(tmp >= root_size){
            root_size = tmp;
            root = gid;
        }
    }

    _lattice.setGroupIDSite(id_last_bond, root);
    _clusters[root].addSite(id_last_bond);

    auto site = _lattice.getSite(id_last_bond);
    auto coordinate_new = site.get_index();
//    auto relative_new = site.relativeIndex();

    auto sites = _lattice.get_neighbor_sites_of_site(id_last_bond);
    IndexRelative dx_dy;
    for(auto n: sites){
        if(_lattice.getGroupIDSite(n) == root){
            // find relative index with respect to this site
            dx_dy = getRelativeIndexDX_v2(n, coordinate_new);
            break; // since first time r is set running loop is doing no good
        }
    }
    _lattice.getRelativeIndex(id_last_bond).add(dx_dy);
    // subtract entropy
    subtract_entropy(gids);
    // insert all to it
    for(auto g: gids){
        if(g == root) continue;
//        cout << "g = " << g << endl;
        _clusters[root].insert(_clusters[g]);

        // relabel
//        relabel(_clusters[g], id_last_site);
        relabel_v3(id_last_bond, sites, _clusters[g]);

        _clusters[g].clear();
    }
    // add entropy
    add_entropy(root);
    // track cluster
    track_clusters(root);
}


/**
 *
 * @param id_current_a : last site id
 * @param neighbors_a  : neighbors sites of `id_current` site
 * @param clstr_b        : cluster to be relabeled
 */
void BondPercolation_pb_v12::relabel_v3(int id_current_a, std::vector<Index>& neighbors_a, Cluster_v12 &clstr_b){
    int gid_current = _lattice.getGroupIDSite(id_current_a);
    auto site = _lattice.getSite(id_current_a);
    auto coordinate_a = site.get_index();
    auto relative_a = site.relativeIndex();
//    cout << "coordinate index of new site " << coordinate_a << endl;
//    cout << "relative index of new site   " << relative_a << endl;
    auto bonds = clstr_b.getBondIDs();
    for(auto b: bonds){
        // relabel bond group id
        _lattice.setGroupIDBond(b, gid_current);
    }

    int gid_b = clstr_b.getGroupID();
    IndexRelative dx_dy;
    for(auto s: neighbors_a){
        auto tmp_id = _lattice.getGroupIDSite(s);
        if(gid_b == tmp_id){
            dx_dy = getRelativeIndexDX_v2(coordinate_a, s);
        }
    }


    auto sites = clstr_b.getSiteIDs();
    for(auto s: sites){
//        cout << "relabeling " << s << endl;
        // relabel site group id
        _lattice.setGroupIDSite(s, gid_current);
        // relabel relative index
        _lattice.getRelativeIndex(s).add(dx_dy);
    }

}

void BondPercolation_pb_v12::reset() {
    Percolation_v12::reset();
    index_counter = 0;
    _number_of_occupied_sites = 0;
    _wrapping_site_ids.clear();
    init();
}

/**
 *
 * @param root     : relative index is known
 * @param site_new : relative index is unknown
 * @return what needs to be added to `site_new` to get it's relative index right
 */
IndexRelative BondPercolation_pb_v12::getRelativeIndexDX(Index root, Index site_new){
//    cout << "Entry \"SitePercolation_ps_v10::getRelativeIndex\" : line " << __LINE__ << endl;
    int delta_x = -int(root.column_) + int(site_new.column_); // if +1 then root is on the right ??
    int delta_y = int(root.row_) - int(site_new.row_); // if +1 then root is on the top ??


    cout << "Delta x,y " << delta_x << ", " << delta_y << endl;
    // normalizing delta_x
    if(delta_x > 1){
        delta_x /= -delta_x;
    }
    else if(delta_x < -1){
        delta_x /= delta_x;
    }

    // normalizing delta_y
    if(delta_y > 1){
        delta_y /= -delta_y;
    }else if(delta_y < -1){
        delta_y /= delta_y;
    }

//    IndexRelative indexRelative_root = _lattice.getSite(root).relativeIndex();
////    cout << "Relative index of root " << indexRelative_root << endl;
//    cout << "Delta x,y " << delta_x << ", " << delta_y << endl;
//    IndexRelative r =  {indexRelative_root.x_ + delta_x, indexRelative_root.y_ + delta_y};
////    cout << "Relative index of site_new " << r << endl;
//    return r;



    return {delta_x, delta_y};
}


/**
 * Just adding the value returned by this function to the relative index of `site_new` will
 * transform it accurately
 * @param root     : relative index is known
 * @param site_new : relative index is unknown
 * @return what needs to be added to `site_new` to get it's relative index right
 */
IndexRelative BondPercolation_pb_v12::getRelativeIndexDX_v2(Index root, Index site_new){
    auto relative_old = _lattice.getRelativeIndex(root);
    auto relative_new = _lattice.getRelativeIndex(site_new);

    int dx_r = relative_old.x_ - relative_new.x_;
    int dy_r = relative_old.y_ - relative_new.y_;

    // since row column is interchanged in natural view
    // row will correspond to x here
    int dx_c = int(site_new.row_) - int(root.row_); // since it can be negative
    int dy_c = int(site_new.column_)- int(root.column_);

    // normalizing delta_x
    if(dx_c > 1){
        // new site is actually below by wrapping property and boundary condition
        // far sites are acually neighbors
        dx_c /= -dx_c;
    }
    else if(dx_c < -1){
        dx_c /= dx_c;
    }

    // normalizing delta_y
    if(dy_c > 1){
        dy_c /= -dy_c;
    }else if(dy_c < -1){
        dy_c /= dy_c;
    }


    int dx = dx_r + dx_c;
    int dy = dy_r + dy_c;

#ifdef DEBUG_FLAG
    cout << "dx_r, dy_r         = " << dx_r << ", " << dy_r << endl;
    cout << "dx_c, dy_c         = " << dx_c << ", " << dy_c << endl;
    cout << "dx, dy             = " << dx << ", " << dy << endl;
#endif
    return {dx, dy};
}

/**
 * New site is inserted in root cluster ?
 * @param id_current
 */
IndexRelative  BondPercolation_pb_v12::relabel_new_site(int id_current) {
    int gid_current = _lattice.getGroupIDSite(id_current);
    auto site = _lattice.getSite(id_current);
    auto coordinate_new = site.get_index();
    auto relative_new = site.relativeIndex();
    cout << "coordinate index of new site " << coordinate_new << endl;
    cout << "relative index of new site   " << relative_new << endl;


    // calculate relative index and relabel sites
    auto sites_neighbor = _lattice.get_neighbor_sites_of_site(id_current);
    // if any of these neighboring sites belong to the same cluster that current site does then
    // we can use that site's relative index and coordinate index to relabel the relative index of merging cluster
    Index coordinate_old ;
    IndexRelative relative_old;

    cout << "gid of current " << gid_current << endl;
    for(auto s: sites_neighbor){
        auto gid = _lattice.getGroupIDSite(s);
        cout <<"neighbor =" << s << "  gid =" << gid << endl;
        if(gid == gid_current){
            relative_old = _lattice.getRelativeIndex(s);
            coordinate_old = s;
            cout << "relative index of old site " << relative_old << endl;
            cout << "coordinate index of old site " << coordinate_old << endl;
            break;
        }
    }

    int dx_r = relative_old.x_ - relative_new.x_;
    int dy_r = relative_old.y_ - relative_new.y_;

    int dy_c = int(coordinate_new.row_) - int(coordinate_old.row_); // since it can be negative
    int dx_c = int(coordinate_new.column_)- int(coordinate_old.column_);

    int dx = dx_r + sign(dx_c);
    int dy = dy_r + sign(dy_c);

    cout << "dx_r, dy_r = " << dx_r << ", " << dy_r << endl;
    cout << "dx_c, dy_c = " << sign(dx_c) << ", " << sign(dy_c) << endl;
    cout << "dx, dy = " << dx << ", " << dy << endl;

    _lattice.getSite(id_current).addToRelativeIndex(dx, dy);
    return IndexRelative{dx, dy};
}

int BondPercolation_pb_v12::sign(int a) {
    if(a == 0) return 0;
    return (a > 0) ? 1 : -1;
}

Index BondPercolation_pb_v12::wrappingSite(){
    if (_wrapping_site_ids.empty()) {
        cout << "No wrapping site" << endl;
        return {};
    }
    return _lattice.getSiteLocation(_wrapping_site_ids[0]);
}

int BondPercolation_pb_v12::wrappingSite_id(){
    if (_wrapping_site_ids.empty()) return -1;
    return _wrapping_site_ids[0];
}

bool BondPercolation_pb_v12::detectWrapping() {
    auto site = _lattice.getSite(id_last_bond);
    auto index = site.get_index();
    // only possible if the cluster containing 'site' has site_index_sequence >= length of the lattice
    if(_number_of_occupied_sites < length()){
        return false;
    }

    // check if it is already a wrapping site
    int gid = _lattice.getGroupIDSite(id_last_bond);
    int tmp_id{};
    for (auto i: _wrapping_site_ids){
        tmp_id = _lattice.getSite(i).get_groupID();
        if(gid == tmp_id ){
//            cout << "Already a wrappig cluster : line " << __LINE__ << endl;
            return true;
        }
    }

    // get four neighbors of site always. since wrapping is valid if periodicity is implied
    vector<Index> sites = _lattice.get_neighbor_sites_of_site(index);


    IndexRelative irel = _lattice.getRelativeIndex(index);
//        cout << "pivot's " << site << " relative " << irel << endl;
    IndexRelative b;
    for (auto a:sites){
        if(_lattice.getSite(a).get_groupID() != _lattice.getGroupIDSite(index)){
            // different cluster
            continue;
        }
//            cout << "belongs to the same cluster : line " << __LINE__ << endl;

        b = _lattice.getSite(a).relativeIndex();
//            cout << "neibhbor " << a << " relative " << b << endl;
        if(abs(irel.x_ - b.x_) > 1 || abs(irel.y_ - b.y_) > 1){
//                cout << "Wrapping : line " << __LINE__ << endl;
            _wrapping_site_ids.push_back(id_last_bond);
            return true;
        }
    }

//    cout << "wrapping site_index_sequence " << _wrapping_indices << endl;
    // if %_wrapping_indices is not empty but wrapping is not detected for the current site (%site)
    // that means there is wrapping but not for the %site
    return !_wrapping_site_ids.empty();
}

size_t BondPercolation_pb_v12::wrappingClusterSize() {


    return 0;
}

long double BondPercolation_pb_v12::entropy() {
//    return entropy_v1();
    return entropy_v2();
}

long double BondPercolation_pb_v12::entropy_v1() {
    long double H{}, mu{};
    double n{};
    for(const auto &clstr: _clusters){
        if(clstr.empty()) continue;
        n = clstr.numberOfSites();
        mu = n/maxSites();
        H += logl(mu) * mu;
    }
//    return -H;
    _entropy = -H;
    return -H;
}

long double BondPercolation_pb_v12::entropy_v2() {
    return _entropy;
}

double BondPercolation_pb_v12::occupationProbability() {
    return double(_number_of_occupied_bonds)/maxBonds();
}

std::string BondPercolation_pb_v12::getSignature() {
    string s = "sq_lattice_site_percolation";
    s += "_periodic";
    return s;
}

size_t BondPercolation_pb_v12::numberOfBondsInTheWrappingClusters() {
    int gid{-1};
    for(auto id: _wrapping_site_ids){
        gid = _lattice.getGroupIDSite(id);
        break;
    }
    if(gid == -1) return 0;
    return _clusters[gid].numberOfBonds();
}

size_t BondPercolation_pb_v12::numberOfSitesInTheWrappingClusters() {
    int gid{-1};
    for(auto id: _wrapping_site_ids){
        gid = _lattice.getGroupIDSite(id);
        break;
    }
    return _clusters[gid].numberOfSites();
}

void BondPercolation_pb_v12::subtract_entropy(const std::set<int> &gids) {
    long double H{}, nb{}, mu{};
    for(auto g: gids){
        nb = _clusters[g].numberOfSites();
        if(nb <= 0) continue; // empty cluster
        mu = nb/maxSites();
        H += logl(mu)*mu;
    }
    // H is negative. so adding is subtracting
    _entropy += H;
}

void BondPercolation_pb_v12::add_entropy(int root) {
    long double H{}, nb{}, mu{};
    nb = _clusters[root].numberOfSites();
    if(nb <= 0) { // empty cluster
        cerr << "root cluster cannot be empty" << endl;
    }
    mu = nb/maxSites();
    H += logl(mu)*mu;
    // H is negative. so subtracting  is adding
    _entropy -= H;
}

void BondPercolation_pb_v12::track_clusters(int root) {
    auto nb = _clusters[root].numberOfBonds();
    auto ns = _clusters[root].numberOfSites();
    if(nb > _number_of_bonds_in_the_largest_cluster){
        _number_of_bonds_in_the_largest_cluster = nb;
    }
    if(ns > _number_of_sites_in_the_largest_cluster){
        _number_of_sites_in_the_largest_cluster = ns;
    }
}

