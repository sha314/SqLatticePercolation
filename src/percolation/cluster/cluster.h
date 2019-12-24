//
// Created by shahnoor on 10/3/2017.
//

#ifndef SITEPERCOLATION_CLUSTER_H
#define SITEPERCOLATION_CLUSTER_H

#include <vector>
#include <set>
#include "../../lattice/bond.h"
#include "../../types.h"
#include "../../lattice/site.h"


/**
 * Cluster of bonds and sites
 * version 3
 * final goal -> make a template cluster. so that we can use it for Bond cluster or Site cluster
 * root site (bond) is the first site (bond) of the cluster. nedeed for (wrapping) site percolation
 */
class Cluster_v3{
    // contains bond and site
    std::vector<BondIndex>  _bond_index; // BondIndex for indexing bonds
    std::vector<Index>      _site_index; // Site index


    int _creation_time;  // holds the creation birthTime of a cluster object
    int _id;
public:
//    using iterator = std::vector<Bond>::iterator;

    ~Cluster_v3()                           = default;
    Cluster_v3()                            = default;
    Cluster_v3(Cluster_v3&)                 = default;
    Cluster_v3(Cluster_v3&&)                = default;
    Cluster_v3& operator=(const Cluster_v3&)      = default;
    Cluster_v3& operator=(Cluster_v3&&)     = default;

    explicit Cluster_v3(int id){

        _id = id;       // may be modified in the program

        /*
         * Only readable, not modifiable.
         * when time = 0 => only lattice exists and bonds in site percolation, not any site_index_sequence
         * When id = 0, time = 1 => we have placed the first site, hence created a cluster with size greater than 1
         *      Only then Cluster_v3 constructor is called.
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

    void insert(const Cluster_v3& cluster);
    void insert_v2(const Cluster_v3& cluster);
    void insert_with_id_v2(const Cluster_v3& cluster, int id);


    friend std::ostream& operator<<(std::ostream& os, const Cluster_v3& cluster);

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
    bool emptySite() const { return _site_index.empty();}
    bool emptyBond() const { return _bond_index.empty();}
    void clear() {_bond_index.clear(); _site_index.clear(); }
};



/**
 * Cluster of bonds and sites
 * version 4
 * instead of storing Index and BondIndex, just store integers.
 * But keep in mind that these integers are the indices of std::vector<Index>
 * and  std::vector<BondIndex> which must be present in SitePercolation and BondPercolation class.
 * This way the program might be a little faster and memory efficient.
 * (3/5 % less memory use): previously, each Index got 2 variables. each BondIndex got  3 variables.
 *      assuming they are all equal size (for simplifying calculations). now for each site there is
 *      one index and for each bond there is one index, total of 2 variables.
 *      Therefore ((5-2)/5 %) less memory required
 *
 *
 */
class Cluster_v4{
    // contains bond and site
    std::vector<uint>  _bond_index; // BondIndex for indexing bonds
    std::vector<uint>  _site_index; // Site index


    int _creation_time;  // holds the creation birthTime of a cluster object
    int _id;
public:
//    using iterator = std::vector<Bond>::iterator;

    ~Cluster_v4()                           = default;
    Cluster_v4()                            = default;
    Cluster_v4(Cluster_v4&)                 = default;
    Cluster_v4(Cluster_v4&&)                = default;
    Cluster_v4& operator=(const Cluster_v4&)      = default;
    Cluster_v4& operator=(Cluster_v4&&)     = default;

    explicit Cluster_v4(int id){

        _id = id;       // may be modified in the program

        /*
         * Only readable, not modifiable.
         * when time = 0 => only lattice exists and bonds in site percolation, not any site_index_sequence
         * When id = 0, time = 1 => we have placed the first site, hence created a cluster with size greater than 1
         *      Only then Cluster_v3 constructor is called.
         *
         */
        _creation_time = id + 1;       // only readable, not modifiable
    }


    void addSiteIndex(size_t length, Index index );
    void addBondIndex(size_t length, BondIndex index);

    void addSiteIndex(uint i){_site_index.push_back(i);};
    void addBondIndex(uint i){_bond_index.push_back(i);};

    uint lastAddedSite(){return _site_index.back();}
    uint lastAddedBond(){return _bond_index.back();}

    void insert_v2(const Cluster_v4& cluster);
    void insert_with_id_v2(const Cluster_v4& cluster, int id);


    friend std::ostream& operator<<(std::ostream& os, const Cluster_v3& cluster);

    const std::vector<uint>  getBondIndices()  const  {return _bond_index;}
    const std::vector<uint>  getSiteIndices()  const  {return _site_index;}

    value_type numberOfBonds() const { return _bond_index.size();}
    value_type numberOfSites() const { return _site_index.size();}

    int get_ID() const { return _id;}
    void set_ID(int id) { _id = id;}

    int birthTime() const {return _creation_time;}

    uint getRootSite()const{return _site_index[0];} // for site percolation
    uint getRootBond()const{return _bond_index[0];} // for bond percolation
    bool empty() const { return _bond_index.empty() && _site_index.empty();}
    bool emptySite() const { return _site_index.empty();}
    bool emptyBond() const { return _bond_index.empty();}
    void clear() {_bond_index.clear(); _site_index.clear(); }
    friend std::ostream &operator<<(std::ostream &os, const Cluster_v4 &cluster);
};


/**
 * Cluster of bonds and sites
 * compitable for v12
 *
 * root site (bond) is the first site (bond) of the cluster. nedeed for (wrapping) site percolation
 */
class Cluster_v12{
    // contains bond and site
    std::vector<int>  _bond_ids; // BondIndex for indexing bonds
    std::vector<int>  _site_ids; // Site index


    int _gid; // group id of the cluster. All sites and bonds must have this same group id
public:
//    using iterator = std::vector<Bond>::iterator;

    ~Cluster_v12()                           = default;
    Cluster_v12()                            = default;
    Cluster_v12(Cluster_v12&)                 = default;
    Cluster_v12(Cluster_v12&&)                = default;
    Cluster_v12& operator=(const Cluster_v12&)      = default;
    Cluster_v12& operator=(Cluster_v12&&)     = default;

    explicit Cluster_v12(int id){

        _gid = id;       // may be modified in the program

    }


    void addSite(Index index);
    void addBond(Index index);

    void addSite(int id);
    void addBond(int id);

//    int lastAddedSite(){return _site_ids.back();}
//    int lastAddedBond(){return _bond_ids.back();}

    void insert(const Cluster_v12& cluster);

    friend std::ostream& operator<<(std::ostream& os, const Cluster_v12& cluster);

    const std::vector<int>&   getBondIDs()    {return _bond_ids;}
    const std::vector<int>&   getSiteIDs()    {return _site_ids;}

//    const std::vector<int>&   getBondIDs()    {return std::vector<int>(_bond_ids.begin(), _bond_ids.end());}
//    const std::vector<int>&   getSiteIDs()    {return std::vector<int>(_site_ids.begin(), _site_ids.end());}

    value_type numberOfBonds() const { return _bond_ids.size();}
    value_type numberOfSites() const { return _site_ids.size();}

    int getGroupID() const { return _gid;}
    void setGroupID(int id) { _gid = id;}


    int getRootSite()const{return _site_ids[0];} // for site percolation
    int getRootBond()const{return _bond_ids[0];} // for bond percolation
    bool empty() const { return _bond_ids.empty() && _site_ids.empty();}
    bool emptySite() const { return _site_ids.empty();}
    bool emptyBond() const { return _bond_ids.empty();}
    void clear() {_bond_ids.clear(); _site_ids.clear(); }
};


#endif //SITEPERCOLATION_CLUSTER_H
