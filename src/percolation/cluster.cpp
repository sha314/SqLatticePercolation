//
// Created by shahnoor on 10/2/2017.
//


//
// Created by shahnoor on 10/11/2017.
//


#include "cluster.h"

using namespace std;

// add Site index
void Cluster::addSiteIndex(Index index) {
    _site_index.push_back(index);
}

void Cluster::addBondIndex(BondIndex bondIndex) {
    _bond_index.push_back(bondIndex);
}


bool Cluster::isPresent(BondIndex bond) const {
    for (auto a: _bond_index) {
        if (a == bond)
            return true;
    }
    return false;
}

bool Cluster::isPresent(Index site) const {
    for (auto a: _site_index) {
        if (a == site)
            return true;
    }
    return false;
}

bool Cluster::checkPresenceAndErase(BondIndex bond) {
    for (auto it = _bond_index.begin(); it != _bond_index.end(); ++it) {
        if (*it == bond) {
            _bond_index.erase(it);
            return true;
        }
    }
    return false;
}

bool Cluster::checkPresenceAndEraseIf(BondIndex bond, bool flag) {
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

bool Cluster::checkPresenceAndEraseIf(Index bond, bool flag) {
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


bool Cluster::checkPresenceAndErase(Index bond) {
    for (auto it = _site_index.begin(); it != _site_index.end(); ++it) {
        if (*it == bond) {
            _site_index.erase(it);
            return true;
        }
    }
    return false;
}


void Cluster::eraseSite(value_type index) {
    auto it = _site_index.begin();
    it += index;
    _site_index.erase(it);
}

void Cluster::eraseBond(value_type index) {
    auto it = _bond_index.begin();
    it += index;
    _bond_index.erase(it);
}


void Cluster::insert(const std::vector<BondIndex>& bonds){
    _bond_index.reserve(bonds.size());
    for(value_type i{} ; i != bonds.size() ; ++i){
        _bond_index.push_back(bonds[i]);
    }
}

void Cluster::insert(const std::vector<Index>& sites){
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
void Cluster::insert(const Cluster &cluster) {
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
void Cluster::insert_v2(const Cluster &cluster) {
    // older time or smaller time is the creation birthTime of the cluster
//    cout << "Comparing " << _creation_time << " and " << cluster._creation_time;
    _creation_time = _creation_time < cluster._creation_time ? _creation_time : cluster._creation_time;
//    cout << " Keeping " << _creation_time << endl;
    _bond_index.insert(_bond_index.end(), cluster._bond_index.begin(), cluster._bond_index.end());
    _site_index.insert(_site_index.end(), cluster._site_index.begin(), cluster._site_index.end());
}


void Cluster::insert_with_id_v2(const Cluster &cluster, int id) {
    _id = id;
    // older time or smaller time is the creation birthTime of the cluster
//    cout << "Comparing " << _creation_time << " and " << cluster._creation_time;
    _creation_time = _creation_time < cluster._creation_time ? _creation_time : cluster._creation_time;
//    cout << " Keeping " << _creation_time << endl;
    _bond_index.insert(_bond_index.end(), cluster._bond_index.begin(), cluster._bond_index.end());
    _site_index.insert(_site_index.end(), cluster._site_index.begin(), cluster._site_index.end());
}


std::ostream &operator<<(std::ostream &os, const Cluster &cluster) {
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
