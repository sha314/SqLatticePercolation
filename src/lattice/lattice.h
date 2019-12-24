//
// Created by shahnoor on 10/2/2017.
//

#ifndef SITEPERCOLATION_LATTICE_H
#define SITEPERCOLATION_LATTICE_H

#include <vector>
#include <cmath>

#include "../percolation/cluster/cluster.h"
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

    void setRelativeIndex(Index index, IndexRelative ir);
    IndexRelative getRelativeIndex(Index index);

    std::vector<Index> getSites();
    std::vector<BondIndex> getBonds();

/******************************************************************************
 * Get Neighbor from given index
 ******************************************************************************/
    std::vector<Index> get_neighbor_sites(Index site);   // site neighbor of site
    std::vector<BondIndex> get_neighbor_bonds(Index site);   // bond neighbor of site
    std::vector<BondIndex> get_neighbor_bonds(BondIndex site); // bond neighbor of bond
    std::vector<Index> get_neighbor_sites(BondIndex bond);   // sites neighbor of bond.

    void get_neighbors(Index site, std::vector<Index> &site_neighbor, std::vector<BondIndex> &bond_neighbor);
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
 * The square Lattice. In 1D form.
 * There is a translator that translate
 *      (a) 1D index to 2D and vice versa for Index.
 *          (r,c) = r*L + c
 *      (b) 1D index to 2D BondIndex and vice versa.
 *          (t,r,c) = t*L^2 + r*L + c
 *
 *      here, r -> row index
 *            c -> column index
 *            t -> bond type. 0 if horizontal and 1 if vertical
 *
 * Site and Bonds are always present But they will not be counted unless they are activated
 * always return by references, so that values in the class itself is modified
 */
class SqLattice_v2 {
//    std::vector<std::vector<Index>> _clusters;  // only store index in the cluster
    std::vector<Site> _sites;  // holds all the site_index_sequence
    std::vector<Bond_v2> _bonds;  // holds all horizontal bonds and vertical bonds

    bool _bond_resetting_flag=true; // so that we can reset all bonds
    bool _site_resetting_flag=true; // and all site_index_sequence

    value_type _length{};
    value_type _length_squared{};

private:
    void reset_sites();
    void reset_bonds();
public:
    ~SqLattice_v2() = default;
    SqLattice_v2() = default;
    SqLattice_v2(SqLattice_v2&) = default;
    SqLattice_v2(SqLattice_v2&&) = default;
    SqLattice_v2& operator=(const SqLattice_v2&) = default;
    SqLattice_v2& operator=(SqLattice_v2&&) = default;

    explicit SqLattice_v2(value_type length);
    void init(bool activate_bonds, bool activate_sites, bool bond_reset, bool site_reset);

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

    void view_bonds();

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
    void deactivate_bond(BondIndex bond);


    value_type length() const { return  _length;}

//    Site getSite(Index index);
//    Bond get_h_bond(Index set_ID);
//    Bond get_v_bond(Index set_ID);

    Site& getSite(Index index);
//    Site&& getSiteR(Index index);

    Bond_v2& getBond(BondIndex);

    const Site& getSite(Index index) const ;

    void setGroupID(Index index, int group_id);
    void setGroupID(BondIndex index, int group_id);

    void setGroupIDSite(size_t index, int group_id);
    void setGroupIDBond(size_t index, int group_id);

    int getGroupID(Index index);
    int getGroupID(BondIndex index);

    int getGroupIDSite(size_t index);
    int getGroupIDBond(size_t index);

    // relative index is only for sites
    void setRelativeIndex(Index index, IndexRelative ir);
    IndexRelative getRelativeIndex(Index index);

    void setRelativeIndex(size_t index, IndexRelative ir);
    IndexRelative getRelativeIndex(size_t index);

    std::vector<Index> getSites();
    std::vector<BondIndex> getBonds();

    std::vector<index_type> getSiteIndices();
    std::vector<index_type> getBondIndices();

/******************************************************************************
 * Get Neighbor from given index
 ******************************************************************************/
    std::vector<Index> get_neighbor_sites(Index site);   // site neighbor of site
    std::vector<BondIndex> get_neighbor_bonds(Index site);   // bond neighbor of site
    std::vector<BondIndex> get_neighbor_bonds(BondIndex site); // bond neighbor of bond
    std::vector<Index> get_neighbor_sites(BondIndex bond);   // sites neighbor of bond.

    void get_neighbors(Index site, std::vector<Index> &site_neighbor, std::vector<BondIndex> &bond_neighbor);


    /**************************
     * Translator
     */
    index_type translateBondTo1D(BondIndex bondIndex) const ;
    index_type translateSiteTo1D(Index index) const;
    BondIndex translate1DToBond(index_type bi) const;
    Index translate1DToSite(index_type i) const;

};




/*************************
 * Date : 2019.12.23
 * Author : Shahnoor
 *
 * Specially build for v12
 *
 * nested vector have two indices. first is for row and second is for column.
 * and row corresponds to y value and column to x value when the lattice is placed on the
 * upper right part of a 2D grid
 */
class Lattice_v12{
    std::vector<std::vector<Site_v12>> _sites_2d{};  // holds all the site_index_sequence
    std::vector<std::vector<Bond_v12>> _bonds_2d{};  // holds all horizontal bonds

/**
 * Finding site or bond Index by their id.
 * id is used index for the array the the value of the element is the Index
 */
    std::vector<Index> _sites{};  // holds all sites
    std::vector<Index> _bonds{};  // holds all bonds

    bool _bond_resetting_flag=true; // so that we can reset all bonds
    bool _site_resetting_flag=true; // and all site_index_sequence

    int _length{};

public:

    ~Lattice_v12() = default;
    Lattice_v12() = default;
    explicit Lattice_v12(int length);

    void reset(bool reset_all=false);

    Index getSiteLocation(int id){return _sites[id];}
    Index getBondLocation(int id){return _bonds[id];}

    Bond_v12 getBond(int id);
    Site_v12 getSite(int id);

    Bond_v12 getBond(Index index){return  _bonds_2d[index.row_][index.column_];}
    Site_v12 getSite(Index index){return  _sites_2d[index.row_][index.column_];}

    /***************************************
     * I/O functions
     **************************************/
    void view_as_assigned();
    void view_all(); // view lattice bonds and site_index_sequence together
    void view_sites();
    void view_sites_by_id();
    void view_sites_by_gid();
    void view_sites_by_id_gid();
    void view_by_id();
    void view_by_gid();
    void view_sites_by_relative_index();
    void view_by_relative_index();
    void view_sites_list();
    void view_bonds_list();
    void print_box(int i, int j);




    /************************************
     * Activation functions
     ***********************************/
    void activateAllSite();
    void activateAllBond();

    void activateSite(Index index){_sites_2d[index.row_][index.column_].activate();}
    void activateBond(Index index){_bonds_2d[index.row_][index.column_].activate();}

    void activateSite(int id){auto index = _sites[id]; _sites_2d[index.row_][index.column_].activate();}
    void activateBond(int id){auto index = _bonds[id]; _bonds_2d[index.row_][index.column_].activate();}

    void deactivate_site(Index index);
    void deactivate_bond(Bond bond);


    int length() const { return  _length;}

    Link getBondLink(Index index){return _bonds_2d[index.row_][index.column_].getIndex();}

    int getBondID(Index index){return _bonds_2d[index.row_][index.column_].get_id();}
    int getSiteID(Index index){return _sites_2d[index.row_][index.column_].get_id();}

    void setGroupIDBond(Index index, int group_id);
    void setGroupIDSite(Index index, int group_id);

    void setGroupIDBond(int id, int group_id);
    void setGroupIDSite(int id, int group_id);

    int getGroupIDBond(int id);
    int getGroupIDSite(int id);


    void setRelativeIndex(Index index, IndexRelative ir);
    IndexRelative getRelativeIndex(Index index);

    std::vector<Index> getSites() const {return _sites;};
    std::vector<Index> getBonds() const {return _bonds;};

/******************************************************************************
 * Get Neighbor from given index
 ******************************************************************************/
    std::vector<Index> get_neighbor_sites(Index site);   // site neighbor of site
    std::vector<BondIndex> get_neighbor_bonds(Index site);   // bond neighbor of site
    std::vector<BondIndex> get_neighbor_bonds(BondIndex site); // bond neighbor of bond
    std::vector<Index> get_neighbor_sites(BondIndex bond);   // sites neighbor of bond.

    void get_neighbors(Index site, std::vector<Index> &site_neighbor, std::vector<Index> &bond_neighbor);

    std::string get_string(const IndexRelative& index)const;
    std::string get_string(const Link &lnk) const;
    std::string get_string(const Index& index) const;
};




#endif //SITEPERCOLATION_LATTICE_H
