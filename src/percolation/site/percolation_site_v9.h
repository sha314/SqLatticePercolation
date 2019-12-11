//
// Created by shahnoor on 6/28/19.
//



#ifndef SQLATTICEPERCOLATION_PERCOLATION_SITE_V9_H
#define SQLATTICEPERCOLATION_PERCOLATION_SITE_V9_H

#include "../percolation.h"


/**
 * Site Percolation by Placing Sites
 *
 * version 9
 *
 * First it randomizes the site index list then use it.
 * Paradigm Shift:
 * Does not delete cluster only makes it empty so that index and id remains the same.
 * This way Searching for index of the cluster using id can be omitted.
 *
 * Feature :
 * 1. Can turn on and off both horizontal and boundary condition
 *
 * 2. Uses class Cluster_v2 for storing clusters
 *
 * 3. Uses Group_ID for Bonds and Sites to identify that they are in the same cluster
 *
 * 4. Occupation probability is calculated by sites,
 *      i.e., number of active sites divided by total number of sites
 *
 * 5. Spanning is calculated by number of bonds in a spanning clusters with periodicity turned off,
 *      i.e., number of bonds in the spanning clusters divided by total number of bonds
 *
 * 6. Unweighted relabeling is ommited in this version ??
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
class SitePercolation_ps_v9 : public SqLatticePercolation{
private:


protected:
    /*
     * invariant property of lattice
     */

    // flags to manipulate method
    bool _periodicity{false};

    value_type min_index; // minimum index = 0
    value_type max_index; // maximum index = length - 1


    // index sequence
    std::vector<Index> index_sequence;  // initialized once
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

    std::vector<Index> _spanning_sites;
    std::vector<Index> _wrapping_sites;
//    std::set<Site, ClusterIDComparator> _spanning_sites; // so that we don't have repeated values

    std::vector<value_type> number_of_sites_to_span;
    std::vector<value_type> number_of_bonds_to_span;
//    std::unordered_set<int> spanning_cluster_ids;
//    vector<int> wrapping_cluster_ids; // todo make global to the class


    /*************************************
     * Logging variable
     ************************************/
    bool _logging_flag{false};
    value_type _total_relabeling{};

    /*****************************************
     * Private Methods
     ******************************************/

    void relabel_sites(const std::vector<Index> &sites, int id_a, int delta_x_ab, int delta_y_ab) ;


    double time_relabel{};
public:
    static constexpr const char* signature = "SitePercolation_ps_v8";

    ~SitePercolation_ps_v9() = default;
    SitePercolation_ps_v9() = default;
    SitePercolation_ps_v9(SitePercolation_ps_v9 & ) = default;
    SitePercolation_ps_v9(SitePercolation_ps_v9 && ) = default;
    SitePercolation_ps_v9(value_type length, bool periodicity=true);

    SitePercolation_ps_v9& operator=(SitePercolation_ps_v9 & ) = default;
//    SitePercolation_ps_v8&& operator=(SitePercolation_ps_v8 && ) = default;
    double get_relabeling_time() {return time_relabel;}
    value_type relabeling_count() const {return _total_relabeling;}
    void init(){randomize_v2();}
    virtual void reset();



    /*************************************************
     * Flags
     ************************************************/
    bool periodicity() const {return _periodicity;}


    /***********************************************
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
    void subtract_entropy_for_bond(const std::set<value_type> &found_index_set, int base=-1);

    /*************************************************
     * Site placing methods
     *
     ************************************************/


    virtual bool occupy();

    value_type placeSite(Index site,
                         std::vector<Index>& neighbor_sites,
                         std::vector<BondIndex>& neighbor_bonds);
    value_type placeSite(Index site);
    value_type placeSite_weighted(Index site); // uses weighted relabeling by first identifying the largest cluster
    value_type placeSite_weighted(Index site,
                                  std::vector<Index>& neighbor_sites,
                                  std::vector<BondIndex>& neighbor_bonds);


    Index selectSite(); // selecting site


    void connection_v1(Index site, std::vector<Index> &neighbors, std::vector<BondIndex> &bonds);
    void connection_v2(Index site, std::vector<Index> &site_neighbor, std::vector<BondIndex> &bond_neighbor);


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


    /**********************************************
     * Information about current state of Class
     **********************************************/
    double numberOfOccupiedSite() const { return _number_of_occupied_sites;}
    double occupationProbability() const { return double(_number_of_occupied_sites)/maxSites();}
    double spanningProbability() const; // number of bonds in spanning cluster / total number of bonds (2*l*l - 2*l)
    long double entropy(); // the shannon entropy



    double orderParameter() const;  // number of bonds in the largest cluster / total number of bonds
    double orderParameter_v2() const;  // number of bonds in the largest cluster / total number of bonds

    value_type numberOfBondsInTheLargestCluster();
    value_type numberOfBondsInTheLargestCluster_v2();
    value_type numberOfBondsInTheSpanningCluster();
    value_type numberOfSitesInTheLargestCluster();

    double numberOfSitesInTheSpanningClusters()  ;
    double numberOfBondsInTheSpanningClusters()  ;
    value_type numberOfSitesInTheSpanningClusters_v2()  ;
    value_type numberOfBondsInTheSpanningClusters_v2()  ;

    value_type numberOfSitesInTheWrappingClusters()  ;
    value_type numberOfBondsInTheWrappingClusters()  ;

    value_type wrappingClusterSize() {return numberOfBondsInTheWrappingClusters();}


    /***********************************
     * Spanning Detection
     **********************************/

    bool detectSpanning_v5(const Index& site);
    bool detectSpanning_v6(const Index& site);

    void save_index_if_in_boundary_v2(const Index& site);
    bool check_if_id_matches(Index site, const std::vector<Index> &edge);
    bool check_if_id_matches_and_erase(Index site, std::vector<Index> &edge);


//    bool isSpanned() const { return !_spanning_sites.empty();}
    bool isSpanned() const { return !_spanning_sites.empty();}
    void scanEdges();


    /***********************************
     * Wrapping Detection
     **********************************/
    bool detectWrapping();
    bool detect_wrapping_v1();


    /************************************
     *  Tracker
     *  Must be called each time a site is placed
     ***********************************/
    void track_numberOfBondsInLargestCluster();
    void track_numberOfSitesInLargestCluster();

    /*********************************
     * I/O functions
     * Printing Status
     * // todo declare these as constants
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
    void initialize();
    void initialize_index_sequence();
    void randomize_v2(); // better random number generator

    std::set<value_type> find_index_for_placing_new_bonds(const std::vector<Index> &neighbors);
    int find_cluster_index_for_placing_new_bonds(const std::vector<Index> &neighbors, std::set<value_type> &found_indices);

    value_type manage_clusters(
            const std::set<value_type> &found_index_set,
            std::vector<BondIndex> &hv_bonds,
            Index &site);

    value_type manage_clusters(
            const std::set<value_type> &found_index_set,
            std::vector<BondIndex> &hv_bonds,
            Index &site,
            int base_id // since id and index is same
    );

    bool anyActiveSite(value_type r, value_type c, value_type delta);
    bool anyActiveSpanningSite(value_type row, value_type col, value_type delta);

public:
    // on test
//    IndexRelative getRelativeIndex(Index root, Index site_new);

    const std::vector<double> clusterSizeDistribution() const override; // 2019.06.17


};



#endif //SQLATTICEPERCOLATION_PERCOLATION_SITE_V9_H
