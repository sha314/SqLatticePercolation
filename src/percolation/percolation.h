//
// Created by shahnoor on 10/5/2017.
//

#ifndef SITEPERCOLATION_PERCOLATION_H
#define SITEPERCOLATION_PERCOLATION_H


#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <climits>

#include "../flags.h"
#include "../types.h"
#include "../lattice/lattice.h"
#include "../index/index.h"
#include "calculation.h"
#include "site_position.h"
#include "id_index.h"
#include "boundary.h"
#include "../ext_libs/inverse_array.h"


template <typename T>
T sum(const std::vector<T>& vec){
    T a{};
    for(auto b: vec){
        a += b;
    }
    return a;
}


/**
* to sort sites with respect to their group id
*/
struct ClusterIDComparator{
    bool operator()(const Site& s1, const Site& s2){
        return s1.groupID() < s2.groupID();
    }
};


/**
 * The Square Lattice Percolation class
 */
class SqLatticePercolation{
protected:
    value_type  _length;
    value_type  _length_squared;
    value_type _max_bonds;

    // structural variables of lattice
    SqLattice _lattice;
    InverseArray<int> _cluster_index_from_id;   // used to track id and index of cluster

    // cluster
    // todo. on the next version only count number of sites
    // todo. and bonds and keep track of root sites and bonds only
    // todo. root site or bond means the first one that is stored

    std::vector<Cluster_v2> _clusters;   // check and remove reapeted index manually
    // every birthTime we create a cluster we assign an ID for them
    int _cluster_id{};

    value_type min_index{};
    value_type max_index{};

public:
    static constexpr const char* signature = "SqLatticePercolation";

    virtual ~SqLatticePercolation() = default;
    SqLatticePercolation(value_type length);

    bool occupy();


    /*********
     *  I/O functions
     */
    virtual void viewCluster();
    virtual void viewClusterExtended();
    virtual void view_bonds(){
        _lattice.view_bonds();
    }
    virtual void viewLattice(){
        _lattice.view_sites();

    }

    /**
 * Also shows the cluster index of the sites
 */
    virtual void viewLatticeExtended(){
        _lattice.view_sites_extended();

    }

    /**
    * Displays group ids of sites in a matrix form
    */
    virtual void viewLatticeByID(){
        _lattice.view_sites_by_id();
        _lattice.view_bonds_by_id();
    }

    virtual void viewSiteByID(){
        _lattice.view_sites_by_id();
    }

    virtual void viewBondByID(){
        _lattice.view_bonds_by_id();
    }

    virtual void viewSiteByRelativeIndex();

    void occupationProbability();
    void entropy();
    void orderParameter();

};


/**
 * Site Percolation by Placing Sites
 *
 * version 8
 *
 * First it randomizes the site index list then use it
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
class SitePercolation_ps_v8 : public SqLatticePercolation{
private:
    value_type _impure_sites;

protected:
    /*
     * invariant property of lattice
     */

    // flags to manipulate method
    bool _periodicity{false};
//    bool _measure_spanning_cluster_by_bond{false};
//    std::vector<CalculationFlags> _calculation_flags;

    value_type min_index; // minimum index = 0
    value_type max_index; // maximum index = length - 1


    // index sequence
    std::vector<Index> index_sequence;  // initialized once
    std::vector<Index> randomized_index_sequence;
    value_type _number_of_occupied_sites{};
    value_type _index_sequence_position{}; // usually starts from 0. but when there is impurity it starts from non zer0


    // every birthTime we create a cluster we assign an ID for them
    int _cluster_id{};
//    int _impurity_id{-2};   // id of the impure sites

    /// Cluster tracker
    // key      -> id of cluster
    // value    -> index of cluster
    InverseArray<int> _cluster_index_from_id;   // used to track id and index of cluster
    value_type _index_last_modified_cluster{};  // id of the last modified cluster


    //// quantity to calculate
//    std::vector<double> _number_of_occupied_sites;
//    std::vector<double> _entropy;
    // entropy
    double _entropy{}, _entropy_to_add{}, _entropy_to_subtract{};

    // entropy calculation ingrediants
    //
    /**
     * <id of cluster, mu_i*log(mu_i)> pair
     * no need to erase any keys because _cluster_id_set ensures that
     *      only existing id is being referred to.
     */
    std::map<int, double> _cluster_entropy; // entropy of individual cluster. used to calculate entropy
    value_type _bonds_in_cluster_with_size_two_or_more{0};   // total number of bonds in the clusters. all cluster has bonds > 1

    // order parameter calculation ingradiants
    // id of the cluster which has maximum number of bonds. used to calculate order parameter
    value_type _number_of_bonds_in_the_largest_cluster{};
    value_type _number_of_sites_in_the_largest_cluster{};   // might be useful later


    Index _last_placed_site;    // keeps track of last placed site

    /**************
     * Spanning variables
     ************/

//    int _first_spanning_cluster_id{-2};


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
    value_type _sites_required_to_min_span;


    /*************************************
     * Logging variable
     ************************************/
    bool _logging_flag{false};
    value_type _total_relabeling{};


public:
    static constexpr const char* signature = "SitePercolation_ps_v8";

    ~SitePercolation_ps_v8() = default;
    SitePercolation_ps_v8(value_type length, bool periodicity=true);

    virtual void reset();
    void configure(std::vector<Index> site_indices);
//    void markImpureSites();

    /*************************************************
     * Flags
     ************************************************/
    void logging_flag(bool f){_logging_flag = f;};
    void periodicityFlag(bool p){ cout << "deprecated : line " << __LINE__ << endl; _periodicity = p;};

//    void measureSpanningClusterByBond(bool flag){_measure_spanning_cluster_by_bond = flag;};


    /***********************************************
     * Properties of Percolation class
     ***********************************************/
    std::string getSignature();
    value_type length() const { return _length;}
    value_type unoccupiedSite() const { return randomized_index_sequence.size();}

    std::vector<value_type> number_of_site_in_clusters();
    std::vector<value_type> number_of_bonds_in_clusters();

    value_type numberOfBondsInCluster_by_id(value_type id); // todo
    value_type numberOfSitesInCluster_by_id(value_type id); // todo

    void numberOfActiveSites() const {std::cout << "Number of active sites " << _number_of_occupied_sites << std::endl;}
    double activeSites() const { return _number_of_occupied_sites;}

    value_type count_number_of_active_site();

    value_type impuritySite() const { return _impure_sites;}
    double impurityPercentage() const { return _impure_sites/double(_length_squared);}

    int birthTimeOfSpanningCluster() const;
    int birthTimeOfACluster(int id) const;



    /****************************************************************
     * Calculations
     ***************************************************************/
    void calculate_entropy();
    void calculate_occupation_probability();
    void calculate_spanning_probability_by_largest_cluster();
    void calculation_short_cut();    // calculate entropy and order parameter

    value_type box_counting(value_type delta); // todo both at once
    value_type box_counting_spanning(value_type delta) ;
    array<value_type, 2> box_counting_v2(value_type delta);
    void add_entropy_for(value_type index);
    void subtract_entropy_for(const set<value_type>& found_index_set);
    void subtract_entropy_for(const vector<value_type>& found_index);

    /*************************************************
     * Site placing methods
     *
     ************************************************/
    void placeAllSites(value_type step=1){
        placeSites(_length_squared, step);
    }
    void placeSites(value_type n, value_type step=1);

    virtual bool occupy();
    value_type placeSite_v6();    // does not erase any element of untouched_site_index
    value_type placeSite_v7();    // does not erase any element of untouched_site_index
    value_type placeSite_weighted_v8();    // does not erase any element of untouched_site_index
    value_type placeSiteWeightedRelabeling_v9();
    value_type placeSite_v10();    // implementing for wrapping

    bool placeSiteForSpanning();

    value_type placeSiteUntilFirstSpanning_v2();

    void placeSite_sequentially(double p);  // used to calculate quantities using convolution process

    void connection(Index site, std::vector<Index>& neighbors, std::vector<BondIndex>& bonds);
    void connection2(Index site, std::vector<Index>& site_neighbor, std::vector<BondIndex>& bond_neighbor);

    /**
     * For explosive percolation
     */
    value_type placeSite_explosive_UntilFirstSpanning(int rule);
    bool placeSite_explosive(int rule);
    bool placeSite_explosive_product_rule_ForSpanning();
    bool placeSite_explosive_sum_rule_ForSpanning();
    value_type placeSite_explosive_prodduct_rule(); // for explosive percolation
    value_type placeSite_explosive_sum_rule();
    value_type placeSite_explosive_test();
    value_type cluster_length_for_placing_site_product_rule(Index site);
    value_type cluster_length_for_placing_site_product_rule_v2(Index site);
    value_type cluster_length_for_placing_site_sum_rule(Index site);
    value_type cluster_length_for_placing_site_sum_rule_v2(Index site);

    /***
     * todo todo
     * For percolation with weight containing sites:
     *  1. occupy sites randomly
     *  2. if a site is selected twice, increase its weight by 1
     *  3. repeate this until there is no unoccupied sites left, i.e.,
     *      number of occupied sites = maximum number of sites = Length * length
     *  4. For spanning purposes do the same with periodicity turned off
     */
    value_type placeSite_with_weight();


    /*************************************************
     * Relabeling methods
     *************************************************/

    // applicable to weighted relabeling

    value_type relabel(value_type index_1, value_type index_2);
    void relabel_sites(const Cluster_v2&  clstr, int id);  // todo pass cluster as reference
    void relabel_sites_v4(Index root_a, const Cluster_v2& clstr_b);
    void relabel_bonds(const Cluster_v2&  clstr, int id);  // todo


    /**********************************************
     * Information about current state of Class
     **********************************************/
    double numberOfOccupiedSite() const { return _number_of_occupied_sites;}
    double occupationProbability() const { return double(_number_of_occupied_sites)/_length_squared;}
    double spanningProbability() const; // number of bonds in spanning cluster / total number of bonds (2*l*l - 2*l)
    double entropy() const; // the shannon entropy
    double entropy_v2(); // the shannon entropy
    double entropy_v3() const;   // the shannon entropy
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



    double numberOfcluster() const { return _clusters.size();}
    int firstSpanningClusterID() const {return _lattice.getSite(_spanning_sites.front()).groupID();}
    int firstSpanningClusterID_v2() const {if(!_spanning_sites.empty()){
            return _lattice.getSite(_spanning_sites[0]).groupID();
        }
        return -1;
    }


    /***********************************
     * Spanning Detection
     **********************************/

    bool detectSpanning_v5(const Index& site);
    bool detectSpanning_v6(const Index& site);

    void save_index_if_in_boundary_v2(const Index& site);
    bool check_if_id_matches(Index site, const vector<Index> &edge);
    bool check_if_id_matches_and_erase(Index site, vector<Index> &edge);


//    bool isSpanned() const { return !_spanning_sites.empty();}
    bool isSpanned() const { return !_spanning_sites.empty();}
    void scanEdges();


    /***********************************
     * Wrapping Detection
     **********************************/
    bool detectWrapping_v1(Index site);


    /************************************
     *  Tracker
     *  Must be called each time a site is placed
     ***********************************/
    void track_numberOfBondsInLargestCluster();
    void track_numberOfSitesInLargestCluster();
    void track_entropy();

    /*********************************
     * I/O functions
     * Printing Status
     * // todo declare these as constants
     ********************************/
    Index lastPlacedSite() const { return _last_placed_site;}
    void viewCluster_id_index();
    void periodicity_status();

    void spanningIndices() const;


    /*****************************
     * debuggind functions
     */
    void view_cluster_from_ground_up();

    void uncooupied_site_count()const {std::cout << "unactivated site number " << randomized_index_sequence.size() << std::endl;}

    std::string getLatticeIDs() {return _lattice.getLatticeIDs();}

    void view_index_order() {for(auto a: randomized_index_sequence){std::cout << a << ',' << std::endl;}}

    void logs();

    /*************************************
     * Checking metods
     *
     */
    void check(Index current_site);

    // temporary public
protected:
    void initialize();
    void initialize_index_sequence();
    void randomize();

    std::set<value_type> find_index_for_placing_new_bonds_v3(const vector<Index>& neighbors);
    std::vector<value_type> find_index_for_placing_new_bonds_v4(std::vector<Index> neighbors);

    value_type manage_clusters_weighted_v5(   // todo
            const std::set<value_type> &found_index_set,
            std::vector<BondIndex> &hv_bonds,
            Index &site);

    value_type manage_clusters_v4(
            const std::set<value_type> &found_index_set,
            std::vector<BondIndex> &hv_bonds,
            Index &site);

    value_type manage_clusters_v6(
            const std::set<value_type> &found_index_set,
            std::vector<BondIndex> &hv_bonds,
            Index &site);

    value_type manage_clusters_v7(
            const std::set<value_type> &found_index_set,
            std::vector<BondIndex> &hv_bonds,
            Index &site);

    value_type manage_clusters_weighted_v8(    // suitable for spanning detection
            const std::set<value_type> &found_index_set,
            std::vector<BondIndex> &hv_bonds,
            Index &site);

    value_type manage_clusters_v9(
            const std::vector<value_type> &found_index_set,
            std::vector<BondIndex> &hv_bonds,
            Index &site);

    value_type manage_clusters_v10(
            const set<value_type> &found_index_set,
            vector<BondIndex> &hv_bonds,
            Index &site);

    bool anyActiveSite(value_type r, value_type c, value_type delta);
    bool anyActiveSpanningSite(value_type row, value_type col, value_type delta);

public:
     // on test
    IndexRelative getRelativeIndex(Index root, Index site_new);

};


/******************************************************************************
 * Explosive site percolation in square lattice with sum rule and product rule
 */
class SitePercolationExplosive: public SitePercolation_ps_v8{

public:
    ~SitePercolationExplosive() = default;
    SitePercolationExplosive(value_type length);
    std::string getSignature() {
        string s = "sq_lattice_site_percolation_explosive_";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        s += to_string(_length);
        return s;
    }

};

/*******************************************************************************
 * Site Percolation Ballistic Deposition
 *
 * *************************************************************/

/*******************************************************************************
 * definition:
 *
 *
 * specification:
 *
 *
 * performence:
 *      1. takes power law iteration with length with exponent ~2.2 in debug mode
 *          iteration_needed = constant*Length ^ (2.2)
 *      2. time ?? todo
 */
class SitePercolationBallisticDeposition: public SitePercolation_ps_v8{
protected:
    // elements of @indices_tmp will be erased if needed but not of @indices
    std::vector<value_type> indices, indices_tmp;
public:
    static constexpr const char* signature = "SitePercolation_BallisticDeposition_v1";
    ~SitePercolationBallisticDeposition() = default;
    SitePercolationBallisticDeposition(value_type length, bool periodicity);

    virtual bool occupy();

    /************************************
     * Site selection methods
     */
    Index select_site(vector<Index> &sites, vector<BondIndex> &bonds);
    Index select_site_upto_1nn(vector<Index> &sites, vector<BondIndex> &bonds);
    Index select_site_upto_2nn(vector<Index> &sites, vector<BondIndex> &bonds);

    /**
    * General site placing method
    * n = 0 : occpy only central site
    * n = 1 : occpy 1st nearest neighbor if central site is occupied
    * n = 2 : occpy 2st nearest neighbor in the direction of motion if 1st nearest neighbor is occuoied
    */
    value_type placeSite_nn_v0(int n=0); // debugging version

    void reset(); // todo
    void initialize_indices();
//    void randomize_index();

    virtual std::string getSignature() {
        string s = "sq_lattice_site_percolation_ballistic_deposition_";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        return s;
    }


    /***********************************
     * occupy upto 1st nearset neighbor.
     * If the randomly selected site is occupied then select one of the nearest neighor randomly
     * If it is also occupied skip the rest setps and start next iteration Else occupy it
     *
     *
     */
    value_type placeSite_1nn_v0(); // debugging version
    value_type placeSite_1nn_v1();


    /*********************************
     * occupy upto 2nd nearest neighbor.
     * If the randomly selected site is occupied then select one of the nearest neighor randomly
     * If it is also occupied, select the next neighbor in the direction of motion Else occupy it.
     * If the 2nd nearest neighbor in the direction of motion is also occupied then skip the rest of the steps
     *      and start the next iteration
     *
     */

    value_type placeSite_2nn_v0();

};

/***********
 * Only L1
 */
class SitePercolationBallisticDeposition_L1: public SitePercolationBallisticDeposition{
public:
    ~SitePercolationBallisticDeposition_L1() = default;
    SitePercolationBallisticDeposition_L1(value_type length, bool periodicity)
            :SitePercolationBallisticDeposition(length, periodicity){}

    bool occupy() {
        // if no site is available then return false

        if(_number_of_occupied_sites == _length_squared){
            return false;
        }

        try {
//        value_type v = placeSite_1nn_v0(); // debugging version
            value_type v = placeSite_1nn_v1();

            return v != ULLONG_MAX;
        }catch (OccupiedNeighbor& on){
//        on.what();
//        cout << "line : " << __LINE__ << endl;
            return false;
        }

    }

    std::string getSignature() {
        string s = "sq_lattice_site_percolation_ballistic_deposition_L1";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        return s;
    }

};

/*********************
 *
 */
class SitePercolationBallisticDeposition_L2: public SitePercolationBallisticDeposition{
public:
    ~SitePercolationBallisticDeposition_L2() = default;
    SitePercolationBallisticDeposition_L2(value_type length, bool periodicity)
            :SitePercolationBallisticDeposition(length, periodicity){}

    bool occupy() {
        // if no site is available then return false

        if(_number_of_occupied_sites == _length_squared){
            return false;
        }

        try {

            value_type v = placeSite_2nn_v0();

            return v != ULLONG_MAX;
        }catch (OccupiedNeighbor& on){
//        on.what();
//        cout << "line : " << __LINE__ << endl;
            return false;
        }

    }

    std::string getSignature() {
        string s = "sq_lattice_site_percolation_ballistic_deposition_L2";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        return s;
    }

};



/**
 * Bond Percolation by Placing Bonds
 *
 * version 0
 *
 * Feature :
 * 1. Can turn on and off both horizontal and boundary condition
 *
 * 2. Uses class Cluster_v2 for storing clusters
 *
 * 3. Uses Group_ID for Bonds and Sites to identify that they are in the same cluster
 *
 * 4. Occupation probability is calculated by bonds,
 *      i.e., number of active bonds divided by total number of bonds
 *
 * 5. Spanning is calculated by number of sites in a spanning clusters with periodicity turned off,
 *      i.e., number of sites in the spanning clusters divided by total number of sites
 */
class BondPercolation_pb_v0 : public SqLatticePercolation{

    // flags to manipulate method
    bool _periodicity{false};
    bool _measure_spanning_cluster_by_bond{false};


    // index sequence
    std::vector<value_type> indices;
    std::vector<BondIndex> index_sequence; // initialized only once
    std::vector<BondIndex> randomized_index_sequence;
    value_type _index_sequence_position{};
    BondIndex _last_placed_bond;

    double _total_number_of_active_bonds{};
    double max_cluster_length;


    /// Cluster tracker
    // key      -> id of cluster
    // value    -> index of cluster

    value_type _index_last_modified_cluster{};  // id of the last modified cluster
    value_type _bonds_in_cluster_with_size_two_or_more{};

    //// quantity to calculate
    value_type _number_of_occupied_bonds;
    std::vector<double> _entropy;
    std::vector<value_type> number_of_bonds_to_span;
    std::vector<value_type> number_of_sites_to_span;
    value_type _sites_required_to_min_span;
    std::vector<CalculationFlags> _calculation_flags;


    // first value the cluster id and second value is the cluster index
    std::vector<std::pair<value_type, value_type>> _spanning_cluster_indices;

public:
    static constexpr const char* signature = "BondPercolation_v1";
    ~BondPercolation_pb_v0() = default;
    BondPercolation_pb_v0(value_type length, bool periodicity=true);

    std::string getSignature() {
        string s = "sq_lattice_bond_percolation_";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        return s;
    }

    void periodicityFlag(bool p){_periodicity=p;};

    value_type length() const { return _length;}
    void reset();
    void configure(std::vector<BondIndex> bond_indices);


    void calculationFlags(std::vector<CalculationFlags> cf);

    std::vector<double> entropy() const { return _entropy;}
    double occupationProbability() const ;
    std::vector<double> spanningProbability() const ;

    void calculate_entropy();
    void calculate_spanning_probability();
    void calculate_spanning_probability_by_largest_cluster();


    void placeAllBonds(){
        while(occupy());
    }

    value_type placeBond_v0();
    bool occupy();

    void relabel_sites(const Cluster_v2& clstr, int id);
    void relabel_bonds(const Cluster_v2& clstr, int id);


    void numberOfActiveSites() const {std::cout << "Number of active sites " << _total_number_of_active_bonds << std::endl;}
    double activeSites() const { return _total_number_of_active_bonds;}

    value_type count_number_of_active_site();

    /***********************************
     * Spanning Detection
     **********************************/
//    value_type numberOfSiteInSpanningClusters(std::vector<int> g_ids);
    value_type number_of_site_in_spanning_clusters(std::unordered_set<int> g_ids);

    bool detectSpanning();

    /*********************************
     * Printing Status
     * // todo declare these as constants
     ********************************/
    void periodicity_status();


    BondIndex lastPlacedBond() {
        if (_index_sequence_position == 0) {
            cout << "Empty lattice : line " << __LINE__ << endl;
            return randomized_index_sequence[_index_sequence_position];
        }
        return randomized_index_sequence[_index_sequence_position-1];
    }

    // debuggind functions


private:
    void initialize();
    void initialize_index_sequence();
    /**
    * Called only once when the object is constructed for the first time
    */
    void initialize_indices();

    void randomize();

//    std::vector<Index> get_Sites_for_bonds(BondIndex);


    value_type manage_clusters_v1(
            const set<value_type> &found_index_set,
            vector<Index> &sites,
            BondIndex &bond
    );

    void connection_v1(BondIndex bond, std::vector<Index> &site_neighbor, std::vector<BondIndex> &bond_neighbor);
    void connection_v2(BondIndex bond, std::vector<Index> &site_neighbor, std::vector<BondIndex> &bond_neighbor);
    void mark_sites(vector<Index> &sites);
    std::set<value_type> find_index_for_placing_new_bonds_v1(const std::vector<BondIndex>& neighbors);


    void connection_2_horizontal_no_periodicity(const BondIndex &bond, vector<Index> &site_neighbor,
                                                vector<BondIndex> &bond_neighbor,
                                                value_type next_column, value_type prev_column,
                                                value_type prev_row);

    void connection_2_vertical_no_peridicity(const BondIndex &bond, vector<Index> &site_neighbor,
                                             vector<BondIndex> &bond_neighbor, value_type prev_column, value_type prev_row,
                                             value_type next_row);

    void connection_2_periodic(const BondIndex &bond, vector<Index> &site_neighbor, vector<BondIndex> &bond_neighbor,
                               value_type next_column, value_type prev_column, value_type prev_row, value_type next_row);
};

#endif //SITEPERCOLATION_PERCOLATION_H

