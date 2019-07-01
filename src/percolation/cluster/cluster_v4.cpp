//
// Created by shahnoor on 6/28/19.
//


#include "cluster.h"

using namespace std;

// add Site index
void Cluster_v4::addSiteIndex(size_t length, Index index) {
    index_type i = IndexTranslator::translateSiteTo1D(length, index);
    _site_index.push_back(i);
}

void Cluster_v4::addBondIndex(size_t length, BondIndex bondIndex) {
    index_type i = IndexTranslator::translateBondTo1D(length, bondIndex);
    _bond_index.push_back(i);
}

//
//bool Cluster_v4::isPresent(BondIndex bond) const {
//    for (auto a: _bond_index) {
//        if (a == bond)
//            return true;
//    }
//    return false;
//}
//
//bool Cluster_v4::isPresent(Index site) const {
//    for (auto a: _site_index) {
//        if (a == site)
//            return true;
//    }
//    return false;
//}
//
//bool Cluster_v4::checkPresenceAndErase(BondIndex bond) {
//    for (auto it = _bond_index.begin(); it != _bond_index.end(); ++it) {
//        if (*it == bond) {
//            _bond_index.erase(it);
//            return true;
//        }
//    }
//    return false;
//}
//
//bool Cluster_v4::checkPresenceAndEraseIf(BondIndex bond, bool flag) {
//    if(!flag)
//        return false;
//    for (auto it = _bond_index.begin(); it != _bond_index.end(); ++it) {
//        if (*it == bond) {
//            _bond_index.erase(it);
//            return true;
//        }
//    }
//    return false;
//}
//
//bool Cluster_v4::checkPresenceAndEraseIf(Index bond, bool flag) {
//    if(!flag)
//        return false;
//    for (auto it = _site_index.begin(); it != _site_index.end(); ++it) {
//        if (*it == bond) {
//            _site_index.erase(it);
//            return true;
//        }
//    }
//    return false;
//}
//
//
//bool Cluster_v4::checkPresenceAndErase(Index bond) {
//    for (auto it = _site_index.begin(); it != _site_index.end(); ++it) {
//        if (*it == bond) {
//            _site_index.erase(it);
//            return true;
//        }
//    }
//    return false;
//}
//
//
//void Cluster_v4::eraseSite(value_type index) {
//    auto it = _site_index.begin();
//    it += index;
//    _site_index.erase(it);
//}
//
//void Cluster_v4::eraseBond(value_type index) {
//    auto it = _bond_index.begin();
//    it += index;
//    _bond_index.erase(it);
//}



/**
 * Merge two cluster as one
 * All intrinsic property should be considered, e.g., creation time of a cluster must be recalculated
 * @param cluster
 */
void Cluster_v4::insert_v2(const Cluster_v4 &cluster) {
    // older time or smaller time is the creation birthTime of the cluster
//    cout << "Comparing " << _creation_time << " and " << cluster._creation_time;
//    _creation_time = _creation_time < cluster._creation_time ? _creation_time : cluster._creation_time;
//    cout << " Keeping " << _creation_time << endl;
    _bond_index.insert(_bond_index.end(), cluster._bond_index.begin(), cluster._bond_index.end());
    _site_index.insert(_site_index.end(), cluster._site_index.begin(), cluster._site_index.end());
}


void Cluster_v4::insert_with_id_v2(const Cluster_v4 &cluster, int id) {
    _id = id;
    // older time or smaller time is the creation birthTime of the cluster
//    cout << "Comparing " << _creation_time << " and " << cluster._creation_time;
    _creation_time = _creation_time < cluster._creation_time ? _creation_time : cluster._creation_time;
//    cout << " Keeping " << _creation_time << endl;
    _bond_index.insert(_bond_index.end(), cluster._bond_index.begin(), cluster._bond_index.end());
    _site_index.insert(_site_index.end(), cluster._site_index.begin(), cluster._site_index.end());
}


std::ostream &operator<<(std::ostream &os, const Cluster_v4 &cluster) {
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

