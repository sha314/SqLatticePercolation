//
// Created by shahnoor on 10/11/2017.
//


#include "cluster.h"

using namespace std;

// add Site index
void Cluster_v2::addSiteIndex(Index index) {
    _site_index.push_back(index);
}

void Cluster_v2::addBondIndex(BondIndex bondIndex) {
    _bond_index.push_back(bondIndex);
}


bool Cluster_v2::isPresent(BondIndex bond) const {
    for (auto a: _bond_index) {
        if (a == bond)
            return true;
    }
    return false;
}

bool Cluster_v2::isPresent(Index site) const {
    for (auto a: _site_index) {
        if (a == site)
            return true;
    }
    return false;
}

bool Cluster_v2::checkPresenceAndErase(BondIndex bond) {
    for (auto it = _bond_index.begin(); it != _bond_index.end(); ++it) {
        if (*it == bond) {
            _bond_index.erase(it);
            return true;
        }
    }
    return false;
}

bool Cluster_v2::checkPresenceAndEraseIf(BondIndex bond, bool flag) {
    if(!flag)
        return false;
    for (auto it = _bond_index.begin(); it != _bond_index.end(); ++it) {
        if (*it == bond) {
            _bond_index.erase(it);
            return true;
        }
    }
    return false;
}

bool Cluster_v2::checkPresenceAndEraseIf(Index bond, bool flag) {
    if(!flag)
        return false;
    for (auto it = _site_index.begin(); it != _site_index.end(); ++it) {
        if (*it == bond) {
            _site_index.erase(it);
            return true;
        }
    }
    return false;
}


bool Cluster_v2::checkPresenceAndErase(Index bond) {
    for (auto it = _site_index.begin(); it != _site_index.end(); ++it) {
        if (*it == bond) {
            _site_index.erase(it);
            return true;
        }
    }
    return false;
}


void Cluster_v2::eraseSite(value_type index) {
    auto it = _site_index.begin();
    it += index;
    _site_index.erase(it);
}

void Cluster_v2::eraseBond(value_type index) {
    auto it = _bond_index.begin();
    it += index;
    _bond_index.erase(it);
}

void Cluster_v2::insert(const std::vector<Bond> &bonds) {
    _bond_index.reserve(bonds.size());
    for(value_type i{} ; i != bonds.size() ; ++i){
        _bond_index.push_back(bonds[i].ID());
    }
}

void Cluster_v2::insert(const std::vector<Site> &sites) {
    _site_index.reserve(sites.size());
    for(value_type i{} ; i != sites.size() ; ++i){
        _site_index.push_back(sites[i].ID());
    }
}

void Cluster_v2::insert(const std::vector<BondIndex>& bonds){
    _bond_index.reserve(bonds.size());
    for(value_type i{} ; i != bonds.size() ; ++i){
        _bond_index.push_back(bonds[i]);
    }
}

void Cluster_v2::insert(const std::vector<Index>& sites){
    _site_index.reserve(sites.size());
    for(value_type i{} ; i != sites.size() ; ++i){
        _site_index.push_back(sites[i]);
    }
}

/**
 * Merge two cluster as one
 * All intrinsic property should be considered, e.g., creation time of a cluster must be recalculated
 * @param cluster
 */
void Cluster_v2::insert(const Cluster_v2 &cluster) {
    if(_id > cluster._id){
        cout << "_id > cluster._id : line " << __LINE__ << endl;
        _id = cluster._id;
    }
    // older time or smaller time is the creation birthTime of the cluster
//    cout << "Comparing " << _creation_time << " and " << cluster._creation_time;
    _creation_time = _creation_time < cluster._creation_time ? _creation_time : cluster._creation_time;
//    cout << " Keeping " << _creation_time << endl;
    _bond_index.insert(_bond_index.end(), cluster._bond_index.begin(), cluster._bond_index.end());
    _site_index.insert(_site_index.end(), cluster._site_index.begin(), cluster._site_index.end());
}


/**
 * Merge two cluster as one
 * All intrinsic property should be considered, e.g., creation time of a cluster must be recalculated
 * @param cluster
 */
void Cluster_v2::insert_v2(const Cluster_v2 &cluster) {
//    if(_id > cluster._id){
//        cout << "_id > cluster._id : line " << __LINE__ << endl;
//    }
    // older time or smaller time is the creation birthTime of the cluster
//    cout << "Comparing " << _creation_time << " and " << cluster._creation_time;
    _creation_time = _creation_time < cluster._creation_time ? _creation_time : cluster._creation_time;
//    cout << " Keeping " << _creation_time << endl;
    _bond_index.insert(_bond_index.end(), cluster._bond_index.begin(), cluster._bond_index.end());
    _site_index.insert(_site_index.end(), cluster._site_index.begin(), cluster._site_index.end());
}


void Cluster_v2::insert_with_id_v2(const Cluster_v2 &cluster, int id) {
    _id = id;
    // older time or smaller time is the creation birthTime of the cluster
//    cout << "Comparing " << _creation_time << " and " << cluster._creation_time;
    _creation_time = _creation_time < cluster._creation_time ? _creation_time : cluster._creation_time;
//    cout << " Keeping " << _creation_time << endl;
    _bond_index.insert(_bond_index.end(), cluster._bond_index.begin(), cluster._bond_index.end());
    _site_index.insert(_site_index.end(), cluster._site_index.begin(), cluster._site_index.end());
}


std::ostream &operator<<(std::ostream &os, const Cluster_v2 &cluster) {
    os << "Sites : size (" << cluster._site_index.size() << ") : ";
    os << '{';
    for(auto a: cluster._site_index){
        os << a << ',';
    }
    os << '}' << endl;

    os << "Bonds : size (" << cluster._bond_index.size() <<") : ";
    os << '{';
    for(auto a: cluster._bond_index){
        os << a << ',';
    }
    os << '}';

    return os << endl;
}

/**************************************
 *
 * Version 3
 **********/
void Cluster_v3::addSites(value_type count) {
//    if(count > 2){ // for bond and site percolation
//        cerr << "At a time more than 2 site can't place in a cluster" << endl;
//    }

    _number_of_sites += count;
}

void Cluster_v3::addBonds(value_type count) {
//    if(count > 4){ // for bond and site percolation
//        cerr << "At a time more than 4 site can't place in a cluster" << endl;
//    }
    _number_of_bonds += count;
}


std::ostream &operator<<(std::ostream &os, const Cluster_v3 &cluster) {
    os << "set_ID " << cluster._id << endl;
    os << "Root site " << cluster._root_site_index << endl;
    os << "Root bond " << cluster._root_bond_index << endl;
    os << "Number of sites " << cluster._number_of_sites << endl;
    os << "Number of bonds " << cluster._number_of_bonds << endl;
    os << "Birth time " << cluster._creation_time << endl;

    return os << endl;
}