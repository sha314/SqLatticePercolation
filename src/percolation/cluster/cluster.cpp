//
// Created by shahnoor on 10/2/2017.
//


//
// Created by shahnoor on 10/11/2017.
//


#include "cluster.h"

using namespace std;

// add Site index
void Cluster_v3::addSiteIndex(Index index) {
    _site_index.push_back(index);
}

void Cluster_v3::addBondIndex(BondIndex bondIndex) {
    _bond_index.push_back(bondIndex);
}


bool Cluster_v3::isPresent(BondIndex bond) const {
    for (auto a: _bond_index) {
        if (a == bond)
            return true;
    }
    return false;
}

bool Cluster_v3::isPresent(Index site) const {
    for (auto a: _site_index) {
        if (a == site)
            return true;
    }
    return false;
}

bool Cluster_v3::checkPresenceAndErase(BondIndex bond) {
    for (auto it = _bond_index.begin(); it != _bond_index.end(); ++it) {
        if (*it == bond) {
            _bond_index.erase(it);
            return true;
        }
    }
    return false;
}

bool Cluster_v3::checkPresenceAndEraseIf(BondIndex bond, bool flag) {
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

bool Cluster_v3::checkPresenceAndEraseIf(Index bond, bool flag) {
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


bool Cluster_v3::checkPresenceAndErase(Index bond) {
    for (auto it = _site_index.begin(); it != _site_index.end(); ++it) {
        if (*it == bond) {
            _site_index.erase(it);
            return true;
        }
    }
    return false;
}


void Cluster_v3::eraseSite(value_type index) {
    auto it = _site_index.begin();
    it += index;
    _site_index.erase(it);
}

void Cluster_v3::eraseBond(value_type index) {
    auto it = _bond_index.begin();
    it += index;
    _bond_index.erase(it);
}


void Cluster_v3::insert(const std::vector<BondIndex>& bonds){
    _bond_index.reserve(bonds.size());
    for(value_type i{} ; i != bonds.size() ; ++i){
        _bond_index.push_back(bonds[i]);
    }
}

void Cluster_v3::insert(const std::vector<Index>& sites){
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
void Cluster_v3::insert(const Cluster_v3 &cluster) {
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
void Cluster_v3::insert_v2(const Cluster_v3 &cluster) {
    // older time or smaller time is the creation birthTime of the cluster
//    cout << "Comparing " << _creation_time << " and " << cluster._creation_time;
//    _creation_time = _creation_time < cluster._creation_time ? _creation_time : cluster._creation_time;
//    cout << " Keeping " << _creation_time << endl;
    _bond_index.insert(_bond_index.end(), cluster._bond_index.begin(), cluster._bond_index.end());
    _site_index.insert(_site_index.end(), cluster._site_index.begin(), cluster._site_index.end());
}


void Cluster_v3::insert_with_id_v2(const Cluster_v3 &cluster, int id) {
    _id = id;
    // older time or smaller time is the creation birthTime of the cluster
//    cout << "Comparing " << _creation_time << " and " << cluster._creation_time;
    _creation_time = _creation_time < cluster._creation_time ? _creation_time : cluster._creation_time;
//    cout << " Keeping " << _creation_time << endl;
    _bond_index.insert(_bond_index.end(), cluster._bond_index.begin(), cluster._bond_index.end());
    _site_index.insert(_site_index.end(), cluster._site_index.begin(), cluster._site_index.end());
}


std::ostream &operator<<(std::ostream &os, const Cluster_v3 &cluster) {
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

void Cluster_v12::addBond(int id) {
    _bond_ids.emplace_back(id);
}

void Cluster_v12::addSite(int id) {
    _site_ids.emplace_back(id);
}

void Cluster_v12::insert(const Cluster_v12 &cluster) {
    if(_gid == cluster._gid){
        cerr << "self merging : Cluster_v12::insert" << endl;
        return;
    }
    _site_ids.insert(_site_ids.end(), cluster._site_ids.begin(), cluster._site_ids.end());
    _bond_ids.insert(_bond_ids.end(), cluster._bond_ids.begin(), cluster._bond_ids.end());
}
