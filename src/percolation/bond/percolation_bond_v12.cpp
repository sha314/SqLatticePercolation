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
    _entropy_site = logl(maxSites());

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
    _lattice.activateBond(id_last_bond);
    ++_number_of_occupied_bonds;
    manageCluster();


    return true;
}

/**
 * Each bond can have only two site neighbors.
 * Select suitable root cluster and other cluster is the one to be merged
 */
void BondPercolation_pb_v12::manageCluster() {// find it's neighbors. sites and then bonds
    auto sites_connected = _lattice.getBond(id_last_bond).connectedSites();
//    cout << "neighbor count " << sites_connected.size() << endl;
    // find which clusters the bonds belong
    // find a suitble root cluster
    int root_gid{-1};
    size_t root_size{}, tmp{};
    set<int> gids; // to prevent repetition
    for(auto ss: sites_connected){
        auto gid = _lattice.getGroupIDSite(ss);
        if(gids.count(gid) > 0) continue; // already taken care of
        gids.emplace(gid);
        tmp = _clusters[gid].numberOfSites();// since sites are always present use this
//        tmp = _clusters[gid].numberOfBonds();
        if(tmp >= root_size){
            root_size = tmp;
            root_gid = gid;
        }
    }

    // subtract entropy before new sites or bonds are added to any clusters. and before merging any clusters
    subtract_entropy(gids);



//    cout << "root gid " << root_gid << endl;
    _lattice.setGroupIDBond(id_last_bond, root_gid);
    _clusters[root_gid].addBond(id_last_bond);


//    int other_id{-1};
//    if(sites_connected[0] == root_id){
//        other_id = sites_connected[1];
//    }else{
//        other_id = sites_connected[0];
//    }
//    cout << "root id " << root_id  << " other_id " << other_id << endl;
//    auto site_a = _lattice.getSite(root_id);
//    auto coordinate_a = site_a.get_index();
//    auto relative_a = site_a.relativeIndex();
//
//    auto site_b = _lattice.getSite(other_id);
//    auto coordinate_b = site_b.get_index();
//    auto relative_b = site_b.relativeIndex();
//
//    auto sites = _lattice.get_neighbor_sites_of_site(id_last_bond);
//    IndexRelative dx_dy = getRelativeIndexDX_v2(coordinate_a, coordinate_b);
//
//    _lattice.getRelativeIndex(coordinate_b).add(dx_dy);


    // insert all to it
    for(auto g: gids){
        if(g == root_gid) continue;
//        cout << "g = " << g << endl;
        _clusters[root_gid].insert(_clusters[g]);

        // relabel
        relabel_v3(id_last_bond, sites_connected, _clusters[g]);

        _clusters[g].clear();
    }
    // add entropy
    add_entropy(root_gid);
    // track cluster
    track_clusters(root_gid);
}


/**
 *
 * @param id_current_bond : last site id
 * @param clstr_b        : cluster to be relabeled
 */
void BondPercolation_pb_v12::relabel_v3(int id_current_bond, vector<int>& id_sites_connected, Cluster_v12 &clstr_b){
    int root_gid = _lattice.getGroupIDBond(id_current_bond);
//    auto bond = _lattice.getBond(id_current_bond);

//    auto coordinate_a = bond.get_index();
//    auto relative_a = bond.relativeIndex();
//    cout << "coordinate index of new site " << coordinate_a << endl;
//    cout << "relative index of new site   " << relative_a << endl;
    auto bonds = clstr_b.getBondIDs();
    for(auto b: bonds){
        // relabel bond group id
        _lattice.setGroupIDBond(b, root_gid);
    }

    int gid_b = clstr_b.getGroupID();
    IndexRelative dx_dy;

    Index root_a{}, other_b{};
    for(auto ss: id_sites_connected){
        if(_lattice.getGroupIDSite(ss) == root_gid){
            root_a = _lattice.getSite(ss).get_index();
        }
        if(_lattice.getGroupIDSite(ss) == gid_b){
            other_b = _lattice.getSite(ss).get_index();
        }
    }

#ifdef DEBUG_FLAG
    cout << "connected sites ids {";
    for(auto ss: id_sites_connected){
        cout << ss << ",";
    }
    cout << "}" << endl;
    cout << "root index " << root_a << " , new site " << other_b << endl;
#endif

    dx_dy = getRelativeIndexDX_v2(root_a, other_b);

#ifdef DEBUG_FLAG
    cout << dx_dy << " will be added as relative index " << endl;
#endif


    auto sites = clstr_b.getSiteIDs();
    for(auto s: sites){
//        cout << "relabeling " << s << endl;
        // relabel site group id
        _lattice.setGroupIDSite(s, root_gid);
        // relabel relative index
        _lattice.getRelativeIndex(s).add(dx_dy);
    }

}

void BondPercolation_pb_v12::reset() {
    Percolation_v12::reset();
    index_counter = 0;
    _number_of_occupied_bonds = 0;
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

    auto relative_root = _lattice.getRelativeIndex(root);
    auto relative_new = _lattice.getRelativeIndex(site_new);



    int dx_r = relative_root.x_ - relative_new.x_;
    int dy_r = relative_root.y_ - relative_new.y_;

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
    cout << "root index " << root << " , new site " << site_new << endl;
    cout << "root relative " << relative_root << " , new relative " << relative_new << endl;
    cout << "dx_r, dy_r         = " << dx_r << ", " << dy_r << endl;
    cout << "dx_c, dy_c         = " << dx_c << ", " << dy_c << endl;
    cout << "dx, dy             = " << dx << ", " << dy << endl;

    cout << IndexRelative{dx, dy} << " will be added as relative index " << endl;
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
    auto bond = _lattice.getBond(id_last_bond);
    auto sites = bond.connectedSites();// tow connected site id

    if(_lattice.getGroupIDSite(sites[0]) == _lattice.getGroupIDSite(sites[1])) {
        auto rel_1 = _lattice.getRelativeIndex(sites[0]);
        auto rel_2 = _lattice.getRelativeIndex(sites[1]);

        if(abs(rel_1.x_ - rel_2.x_) > 1 || abs(rel_1.y_ - rel_2.y_) > 1){
//                cout << "Wrapping : line " << __LINE__ << endl;
            _wrapping_site_ids.push_back(sites[0]);
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
    _entropy_site = -H;
    return -H;
}

long double BondPercolation_pb_v12::entropy_v2() {
    return _entropy_site;
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
    _entropy_site += H;
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
    _entropy_site -= H;
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

