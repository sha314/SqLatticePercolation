//
// Created by shahnoor on 6/28/19.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_SITE_V10_H
#define SQLATTICEPERCOLATION_PERCOLATION_SITE_V10_H

#include "../percolation.h"

/**
 * class    : SitePercolation_ps_v10
 * superclass : SqLatticePercolation
 * used classes :
 *
 * Site Percolation by Placing Sites
 *
 * version 10
 *
 * First it randomizes the site index list then use it.
 *
 * Paradigm Shift:
 * 1. Does not delete cluster only makes it empty so that index and id remains the same.
 *    This way Searching for index of the cluster using id can be omitted.
 * 2. In site percolation, initially all bonds are present and they all form a cluster of unit size.
 *    This is maintained here unlike Version 9
 *
 * Feature :
 * 1. Can turn on and off both horizontal and boundary condition
 *
 * 2. Uses class Cluster for storing clusters
 *
 * 3. Uses Group_ID for Bonds and Sites to identify that they are in the same cluster
 *
 * 4. Occupation probability is calculated by sites,
 *      i.e., number of active sites divided by total number of sites
 *
 * 5. Spanning is calculated by number of bonds in a spanning clusters with periodicity turned off,
 *      i.e., number of bonds in the spanning clusters divided by total number of bonds
 *
 * 6. Unweighted relabeling is ommited in this version
 *
 * 7. Runtime is significantly improved. For example, if L=200 program will take ~1 min to place all sites.
 *
 * 8. Unnecessary methods of previous version is eliminated
 *
 * 9. Checking spanning by keeping track of boundary sites is implemented
 *
 * 10. last modified cluster id can be obtained from @var _last_placed_site
 *
 *
 */
class SitePercolation_ps_v10 : public SqLatticePercolation{
protected:

    // flags to manipulate method
    bool _periodicity{false};

    value_type min_index; // minimum index = 0
    value_type max_index; // maximum index = length - 1

    // index sequence
    std::vector<Index> site_indices;  // initialized once
    std::vector<BondIndex> bond_indices;  // initialized once
//    std::vector<Index> randomized_index_sequence;
    std::vector<value_type> randomized_index;

    // every birthTime we create a cluster we assign an set_ID for them
    int _cluster_id{};

    Index _last_placed_site;    // keeps track of last placed site

    /**************
     * Spanning variables
     ************/

    /*Holds indices on the edges*/
    std::vector<Index> _top_edge, _bottom_edge, _left_edge, _right_edge;

    /*
     * Contains one site of a spanning cluster
     * so the using the site indices saved here we can find spanning cluster ids and indices
     * from lattice.
     * each time one site is saved here the indices of the cluster that it belongs to must be erased from
     * BoundaryLR and BoundaryTB containers. It would make the spanning detection process faster. todo
     */

    bool _spanning_occured = false;

    std::vector<Index> _spanning_sites; //todo : remove one of these
    std::vector<Index> _wrapping_sites;

    std::vector<value_type> number_of_sites_to_span;
    std::vector<value_type> number_of_bonds_to_span;

    /*************************************
     * Logging variable
     ************************************/
    bool _logging_flag{false};
    value_type _total_relabeling{};
    double time_relabel{};

    void relabel_sites(const std::vector<Index> &sites, int id_a, int delta_x_ab, int delta_y_ab) ;

public:
    static constexpr const char* signature = "SitePercolation_ps_v8";

    ~SitePercolation_ps_v10() = default;
    SitePercolation_ps_v10() = default;
    SitePercolation_ps_v10(SitePercolation_ps_v10 & ) = default;
    SitePercolation_ps_v10(SitePercolation_ps_v10 && ) = default;
    SitePercolation_ps_v10(value_type length, bool periodicity=true);

    SitePercolation_ps_v10& operator=(SitePercolation_ps_v10 & ) = default;
//    SitePercolation_ps_v8&& operator=(SitePercolation_ps_v8 && ) = default;

    void init(bool gen_random=false); // collections of statements that must be called before using

    double get_relabeling_time() {return time_relabel;}
    value_type relabeling_count() const {return _total_relabeling;}

    virtual void reset();


    /*************************************************
     * Flags
     ************************************************/
    bool periodicity() const {return _periodicity;}


    virtual /***********************************************
     * Properties of Percolation class
     ***********************************************/
    std::string getSignature();

    void numberOfActiveSites() const {std::cout << "Number of active site_index_sequence " << _number_of_occupied_sites << std::endl;}
    double activeSites() const { return _number_of_occupied_sites;}

    value_type count_number_of_active_site();

    int birthTimeOfSpanningCluster() const;
    int birthTimeOfACluster(int id) const;


    /****************************************************************
     * Calculations
     ***************************************************************/

    value_type box_counting(value_type delta); // todo both at once
    value_type box_counting_spanning(value_type delta) ;
    std::array<value_type, 2> box_counting_v2(value_type delta);

    void add_entropy_for_bond(value_type index);
    void subtract_entropy_for_bond(const std::vector<BondIndex> bonds);

    /*************************************************
     * Site placing methods
     *
     ************************************************/


    virtual bool occupy();

    value_type placeSite(Index site);
    value_type placeSite(Index site, const std::vector<Index>& sites, const std::vector<BondIndex> &bonds);


    value_type manageClusters(const std::vector<Index>& sites, const std::vector<BondIndex> &bonds);

    Index selectSite(); // selecting site


    void connection_v1(Index site, std::vector<Index> &neighbors, std::vector<BondIndex> &bonds);
    void connection_v2(Index site, std::vector<Index> &site_neighbor, std::vector<BondIndex> &bond_neighbor);
    void connection_periodic(Index site,
                             std::vector<Index> &site_neighbor,
                             std::vector<BondIndex> &bond_neighbor); // 2019.01.04

    /*************************************************
     * Relabeling methods
     *************************************************/

    // applicable to weighted relabeling

    value_type relabel(value_type index_1, value_type index_2);
    void relabel_sites(const Cluster_v3&  clstr, int id);  // todo pass cluster as reference
    void relabel_sites_v4(Index root_a, const Cluster_v3& clstr_b); // relative index is set accordingly
    void relabel_sites_v5(Index root_a, const Cluster_v3& clstr_b); // relative index is set accordingly
    void relabel_sites_v6(Index root_a, const Cluster_v3& clstr_b, int id); // relative index is set accordingly
    void relabel_bonds(const Cluster_v3&  clstr, int id);  // todo


    // relabel site_index_sequence and bonds in the cluster cluster
    void relabel_cluster(Index site_a, const Cluster_v3& clstr_b, size_t bond_pos=0, size_t site_pos=0);

    /**********************************************
     * Information about current state of Class
     **********************************************/
    double numberOfOccupiedSite() const { return _number_of_occupied_sites;}
    double occupationProbability() const { return double(_number_of_occupied_sites)/maxSites();}
    double spanningProbability() const; // number of bonds in spanning cluster / total number of bonds (2*l*l - 2*l)
    long double entropy(); // the shannon entropy
    long double entropy_v1(); // the shannon entropy
    long double entropy_v2(); // the shannon entropy


    double orderParameter() const;  // number of bonds in the largest cluster / total number of bonds
    double orderParameter_v2() const;  // number of bonds in the largest cluster / total number of bonds

    value_type numberOfBondsInTheLargestCluster();
    value_type numberOfSitesInTheLargestCluster();

    value_type numberOfSitesInTheSpanningClusters()  ;
    value_type numberOfBondsInTheSpanningClusters()  ;

    value_type numberOfSitesInTheWrappingClusters()  ;
    value_type numberOfBondsInTheWrappingClusters()  ;

    value_type wrappingClusterSize() {return numberOfBondsInTheWrappingClusters();}


    /***********************************
     * Spanning and Wrapping
     **********************************/
    bool detectSpanning_v5();
    bool detectSpanning_v6();

    void save_index_if_in_boundary_v2(const Index& site);
    bool check_if_id_matches(Index site, const std::vector<Index> &edge);
    bool check_if_id_matches_and_erase(Index site, std::vector<Index> &edge);


//    bool isSpanned() const { return !_spanning_sites.empty();}
    bool isSpanned() const { return !_spanning_sites.empty();}
    void scanEdges();

    bool detectWrapping();
    bool detect_wrapping_v1();


    /************************************
     *  Tracker
     *  Must be called each time a site is placed
     ***********************************/
    void track_numberOfBondsInLargestCluster() override;
    void track_numberOfSitesInLargestCluster() override;

    /*********************************
     * I/O functions
     * Printing Status
     ********************************/
    Index lastPlacedSite() const { return _last_placed_site;}

    void spanningIndices() const;
    void wrappingIndices() const;

    /***********************************************
     * Visual data for plotting
     *********************************************/
    // lattice visual data for python
    void writeVisualLatticeData(const std::string& filename, bool only_spanning=true);

protected:
    void initialize(); // collections of statements that must be called once while constructing
    void initialize_index_sequence();
    void initialize_cluster();
    void randomize_v2(); // better random number generator

    value_type placeSite_weighted_v2(Index site);
    value_type placeSite_weighted_v3(Index site, const std::vector<Index>& sites, const std::vector<BondIndex> &bonds);
    value_type merge_cluster_v2(const std::vector<Index>& sites, const std::vector<BondIndex> &bonds);
    value_type merge_cluster_v3(const std::vector<Index>& sites, const std::vector<BondIndex> &bonds);

    bool anyActiveSite(value_type r, value_type c, value_type delta);
    bool anyActiveSpanningSite(value_type row, value_type col, value_type delta);

public:
    // on test
//    IndexRelative getRelativeIndex(Index root, Index site_new);

    const std::vector<double> clusterSizeDistribution() const override; // 2019.06.17

    void viewRandomizedIndices();
    std::string getClassName(){ return "SitePercolation_ps_v10";}
};




#endif //SQLATTICEPERCOLATION_PERCOLATION_SITE_V10_H
