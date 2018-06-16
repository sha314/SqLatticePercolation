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
 */
class Cluster_v2{
    // contains bond and site ??
    std::vector<BondIndex>  _bond_index; // BondIndex for indexing bonds
    std::vector<Index>      _site_index; // Site index


//    Index _root_site; // root site is the first site of the cluster. nedeed for (wrapping) site percolation
//    BondIndex _root_bond; // root bond is the first bond of the cluster.

    int _creation_time;  // holds the creation birthTime of a cluster object
    int _id;
public:
//    using iterator = std::vector<Bond>::iterator;

    ~Cluster_v2()                           = default;
    Cluster_v2()                            = default;
    Cluster_v2(Cluster_v2&)                 = default;
    Cluster_v2(Cluster_v2&&)                = default;
    Cluster_v2& operator=(Cluster_v2&)      = default;
    Cluster_v2& operator=(Cluster_v2&&)     = default;

    Cluster_v2(int id){
//        _bond_index.reserve(max_size);
        /*
         *
         */
        _id = id;       // may be modified in the program

        /*
         * Only readable, not modifiable.
         * when time = 0 => only lattice exists and bonds in site percolation, not any sites
         * When id = 0, time = 1 => we have placed the first site, hence created a cluster with size greater than 1
         *      Only then Cluster_v2 constructor is called.
         *
         */
        _creation_time = id + 1;       // only readable, not modifiable
    }


//    void addElement(Bond ); // add an element to the cluster
//    void addElement(Site ); // add an element to the cluster

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

    void insert(const std::vector<Bond>& bonds);
    void insert(const std::vector<Site>& sites);

    void insert(const std::vector<BondIndex>& bonds);
    void insert(const std::vector<Index>& sites);

    void insert(const Cluster_v2& cluster);
    void insert_v2(const Cluster_v2& cluster);

//    void insert_weighted_relabling(std::vector<Bond>& bonds);
//    void insert_weighted_relabling(std::vector<Site>& sites);
//
//    void insert_weighted_relabling(Cluster_v2& cluster);


    friend std::ostream& operator<<(std::ostream& os, const Cluster_v2& cluster);

    const std::vector<BondIndex>&   getBondIndices()    {return _bond_index;}
    const std::vector<Index>&       getSiteIndices()    {return _site_index;}

    const std::vector<BondIndex>&   getBondIndices()  const  {return _bond_index;}
    const std::vector<Index>&       getSiteIndices()  const  {return _site_index;}

    value_type numberOfBonds() const { return _bond_index.size();}
    value_type numberOfSites() const { return _site_index.size();}
    int ID() const { return _id;}
    void ID(int id) { _id = id;}

    int birthTime() const {return _creation_time;}

    Index getRootSite()const{return _site_index[0];} // for site percolation
    BondIndex getRootBond()const{return _bond_index[0];} // for bond percolation
};



/**
 * Does not keep track of sites or bonds that has the same id.
 * So not really useful when relabeling sites or bonds
 *
 * version 3
 */
class Cluster_v3{
    // contains bond and site ??
    BondIndex  _root_bond_index; // BondIndex for indexing bonds
    Index      _root_site_index; // Site index
    value_type _number_of_bonds{}, _number_of_sites{};

    int _creation_time;  // holds the creation birthTime of a cluster object
    int _id;
public:
//    using iterator = std::vector<Bond>::iterator;

    ~Cluster_v3()                           = default;
    Cluster_v3()                            = default;
    Cluster_v3(Cluster_v3&)                 = default;
    Cluster_v3(Cluster_v3&&)                = default;
    Cluster_v3& operator=(Cluster_v3&)      = default;
    Cluster_v3& operator=(Cluster_v3&&)     = default;

    Cluster_v3(int id){
//        _bond_index.reserve(max_size);
        /*
         *
         */
        _id = id;       // may be modified in the program

        /*
         * Only readable, not modifiable.
         * when time = 0 => only lattice exists and bonds in site percolation, not any sites
         * When id = 0, time = 1 => we have placed the first site, hence created a cluster with size greater than 1
         *      Only then Cluster_v2 constructor is called.
         *
         */
        _creation_time = id + 1;       // only readable, not modifiable
    }

    void root(Index site, BondIndex bond){
        _root_site_index = site;
        _root_bond_index = bond;
    }

    void addSites(value_type count);
    void addBonds(value_type count);

    friend std::ostream& operator<<(std::ostream& os, const Cluster_v3& cluster);

    value_type numberOfBonds() const { return _number_of_bonds;}
    value_type numberOfSites() const { return _number_of_sites;}

    int ID() const { return _id;}
    void ID(int id) { _id = id;}

    int birthTime() const {return _creation_time;}
};



#endif //SITEPERCOLATION_CLUSTER_H
