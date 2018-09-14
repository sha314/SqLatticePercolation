//
// Created by shahnoor on 10/3/2017.
//

#ifndef SITEPERCOLATION_CLUSTER_H
#define SITEPERCOLATION_CLUSTER_H

#include <vector>
#include <set>
#include "../lattice/bond.h"
#include "../types.h"
#include "../lattice/site.h"


/**
 * Cluster of bonds and sites
 * version 3
 * final goal -> make a template cluster. so that we can use it for Bond cluster or Site cluster
 * root site (bond) is the first site (bond) of the cluster. nedeed for (wrapping) site percolation
 */
class Cluster{
    // contains bond and site
    std::vector<BondIndex>  _bond_index; // BondIndex for indexing bonds
    std::vector<Index>      _site_index; // Site index


    int _creation_time;  // holds the creation birthTime of a cluster object
    int _id;
public:
//    using iterator = std::vector<Bond>::iterator;

    ~Cluster()                           = default;
    Cluster()                            = default;
    Cluster(Cluster&)                 = default;
    Cluster(Cluster&&)                = default;
    Cluster& operator=(const Cluster&)      = default;
    Cluster& operator=(Cluster&&)     = default;

    explicit Cluster(int id){

        _id = id;       // may be modified in the program

        /*
         * Only readable, not modifiable.
         * when time = 0 => only lattice exists and bonds in site percolation, not any sites
         * When id = 0, time = 1 => we have placed the first site, hence created a cluster with size greater than 1
         *      Only then Cluster constructor is called.
         *
         */
        _creation_time = id + 1;       // only readable, not modifiable
    }


    void addSiteIndex(Index );
    void addBondIndex(BondIndex );

    Index lastAddedSite(){return _site_index.back();}
    BondIndex lastAddedBond(){return _bond_index.back();}

    bool isPresent(BondIndex bond) const ;
    bool isPresent(Index site) const ;

    bool checkPresenceAndErase(BondIndex bond);
    bool checkPresenceAndErase(Index bond);

    bool checkPresenceAndEraseIf(BondIndex bond, bool flag);
    bool checkPresenceAndEraseIf(Index bond, bool flag);

    void eraseSite(value_type index);
    void eraseBond(value_type index);


    void insert(const std::vector<BondIndex>& bonds);
    void insert(const std::vector<Index>& sites);

    void insert(const Cluster& cluster);
    void insert_v2(const Cluster& cluster);
    void insert_with_id_v2(const Cluster& cluster, int id);


    friend std::ostream& operator<<(std::ostream& os, const Cluster& cluster);

    const std::vector<BondIndex>&   getBondIndices()    {return _bond_index;}
    const std::vector<Index>&       getSiteIndices()    {return _site_index;}

    const std::vector<BondIndex>&   getBondIndices()  const  {return _bond_index;}
    const std::vector<Index>&       getSiteIndices()  const  {return _site_index;}

    value_type numberOfBonds() const { return _bond_index.size();}
    value_type numberOfSites() const { return _site_index.size();}
    int get_ID() const { return _id;}
    void set_ID(int id) { _id = id;}

    int birthTime() const {return _creation_time;}

    Index getRootSite()const{return _site_index[0];} // for site percolation
    BondIndex getRootBond()const{return _bond_index[0];} // for bond percolation
    bool empty() const { return _bond_index.empty() && _site_index.empty();}
    void clear() {_bond_index.clear(); _site_index.clear(); }
};



#endif //SITEPERCOLATION_CLUSTER_H
