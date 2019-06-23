//
// Created by shahnoor on 10/2/2017.
//

#ifndef SITEPERCOLATION_LATTICE_H
#define SITEPERCOLATION_LATTICE_H

#include <vector>
#include <cmath>

#include "../percolation/cluster.h"
#include "../types.h"
#include "../flags.h"
#include "site.h"
#include "bond.h"


/**
 * The square Lattice
 * Site and Bonds are always present But they will not be counted unless they are activated
 * always return by references, so that values in the class itself is modified
 */
class SqLattice {
//    std::vector<std::vector<Index>> _clusters;  // only store index in the cluster
    std::vector<std::vector<Site>> _sites;  // holds all the site_index_sequence
    std::vector<std::vector<Bond>> _h_bonds;  // holds all horizontal bonds
    std::vector<std::vector<Bond>> _v_bonds;  // holds all vertical bonds

    bool _bond_resetting_flag=true; // so that we can reset all bonds
    bool _site_resetting_flag=true; // and all site_index_sequence

    value_type _length{};

private:
    void reset_sites();
    void reset_bonds();
public:
    ~SqLattice() = default;
    SqLattice() = default;
    SqLattice(SqLattice&) = default;
    SqLattice(SqLattice&&) = default;
    SqLattice& operator=(const SqLattice&) = default;
    SqLattice& operator=(SqLattice&&) = default;

    SqLattice(value_type length, bool activate_bonds, bool activate_sites, bool bond_reset, bool site_reset);

    void reset(bool reset_all=false);


    /***************************************
     * I/O functions
     **************************************/
    void view_sites();
    void view_sites_extended();
    void view_sites_by_id();
    void view_sites_by_relative_index();
    void view_bonds_by_relative_index();
    void view_bonds_by_relative_index_v2();
    void view_bonds_by_relative_index_v3();
    void view_bonds_by_relative_index_v4();
    void view_by_relative_index();
    void view(); // view lattice bonds and site_index_sequence together

    void view_h_bonds();
    void view_v_bonds();

    void view_bonds(){
        view_h_bonds();
        view_v_bonds();
    }

    void view_h_bonds_extended();
    void view_v_bonds_extended();

    void view_bonds_by_id();

    /************************************
     * Activation functions
     ***********************************/
    void activateAllSite();
    void activateAllBond();
    void activate_site(Index index);
    void activateBond(BondIndex bond);

    void deactivate_site(Index index);
    void deactivate_bond(Bond bond);


    value_type length() const { return  _length;}

//    Site getSite(Index index);
//    Bond get_h_bond(Index set_ID);
//    Bond get_v_bond(Index set_ID);

    Site& getSite(Index index);
//    Site&& getSiteR(Index index);
    Bond& get_h_bond(Index id);
    Bond& get_v_bond(Index id);
    Bond& getBond(BondIndex);

    const Site& getSite(Index index) const ;

    void setGroupID(Index index, int group_id);
    void setGroupID(BondIndex index, int group_id);
    int getGroupID(Index index);
    int getGroupID(BondIndex index);

    std::vector<Index> getSites();
    std::vector<BondIndex> getBonds();

/******************************************************************************
 * Get Neighbor from given index
 ******************************************************************************/
    std::vector<Index> get_neighbor_sites(Index site);   // site neighbor of site
    std::vector<BondIndex> get_neighbor_bonds(Index site);   // bond neighbor of site
    std::vector<BondIndex> get_neighbor_bonds(BondIndex site); // bond neighbor of bond
    std::vector<Index> get_neighbor_sites(BondIndex bond);   // sites neighbor of bond.

//    static std::vector<Index> get_neighbor_sites(size_t length, Index site);   // 4 site neighbor of site
//    static std::vector<BondIndex> get_neighbor_bond_indices(size_t length, BondIndex site); // 6 bond neighbor of bond
//    static std::vector<Index> get_neighbor_indices(size_t length, BondIndex bond);   // 2 site neighbor of bond.

//    std::vector<Index> getNeighborSite(Index site, bool periodicity=false);
//    std::vector<Index> getNeighborSite(BondIndex bond, bool periodicity=false);
//
//    std::vector<BondIndex> getNeighborBond(Index site, bool periodicity=false);
//    std::vector<BondIndex> getNeighborBond(BondIndex site, bool periodicity=false);
//
//    void connection_v1(Index site, std::vector<Index>& neighbors, std::vector<BondIndex>& bonds, bool periodicity);
//    void connection_v1(BondIndex bond, std::vector<Index>& neighbors, std::vector<BondIndex>& bonds, bool periodicity);


};


/**
 * Weighted Planner Stochastic Lattice
 */
class WPSLattice{


public:
    WPSLattice() = default;
    ~WPSLattice() = default;

};



#endif //SITEPERCOLATION_LATTICE_H
