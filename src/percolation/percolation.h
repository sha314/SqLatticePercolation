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
#include <fstream>

#include "../flags.h"
#include "../types.h"
#include "../lattice/lattice.h"
#include "../index/index.h"
#include "calculation.h"
#include "site_position.h"
#include "id_index.h"
#include "boundary.h"
#include <random>


/**
 * The Square Lattice Percolation class
 */
class SqLatticePercolation{
    // constants
    value_type  _length;
    value_type _max_number_of_bonds;
    value_type _max_number_of_sites;
    char type{'0'}; // percolation type. 's' -> site percolation. 'b' -> bond percolation
protected:

    // structural variables of lattice
    SqLattice _lattice;

    value_type _index_sequence_position{};
    // cluster
    // todo. on the next version only count number of sites
    // todo. and bonds and keep track of root sites and bonds only
    // todo. root site or bond means the first one that is stored

    std::vector<Cluster> _clusters;   // check and remove reapeted index manually
    // every birthTime we create a cluster we assign an set_ID for them
    int _cluster_id{};

    value_type min_index{};
    value_type max_index{};

    double _occuption_probability {};
    // entropy
    double _entropy{};
    double _entropy_current{};
    double _entropy_previous{};
    double _largest_jump_entropy{}; // lrgest jump in entropy
    double _entropy_jump_pc{}; // at what pc there have been the largest jump
    size_t _cluster_count{};

    value_type _total_relabeling{};
    double time_relabel{};
    value_type _number_of_occupied_sites{};
    value_type _number_of_occupied_bonds{};
    value_type _max_iteration_limit{};
    std::vector<value_type> randomized_index;
    std::random_device _random_device;
    std::mt19937 _random_generator;

    void set_type(char t){type = t;} // setting percolation type
public:
    static constexpr const char* signature = "SqLatticePercolation";

    virtual ~SqLatticePercolation() = default;
    SqLatticePercolation(value_type length);
    void reset();
    void set_cluster_measuring_unit(int i){
        std::cout << "Cluster measuring unit = " << ((i==0) ? "bond" : "site")
                  << " : line " << __LINE__ << std::endl;
    }

    bool occupy();
    value_type length() const { return _length;}
    value_type maxSites() const {return _max_number_of_sites;}
    value_type maxBonds() const { return _max_number_of_bonds;}

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

    virtual void viewSiteByRelativeIndex(){
        _lattice.view_sites_by_relative_index();
    }
    virtual void viewBondByRelativeIndex(){
        _lattice.view_bonds_by_relative_index_v4();
    }

    virtual void viewByRelativeIndex(){
        _lattice.view_by_relative_index();
    }

    virtual void view(){
        _lattice.view();
    }

    virtual double occupationProbability() const { return _occuption_probability;}
    virtual double entropy() { return _entropy_current;}
    double entropy_by_site(); // for future convenience
    double entropy_by_bond(); // for future convenience
    double orderParameter();
    size_t numberOfcluster() const {return _cluster_count;}
    void jump();
    double largestEntropyJump()const { return _largest_jump_entropy;}
    double largestEntropyJump_pc()const { return _entropy_jump_pc;}


    void get_cluster_info(
            std::vector<value_type> &site,
            std::vector<value_type> &bond
    );

    char get_type() const {return type;} // get percolation type
    virtual value_type maxIterationLimit() {return _max_iteration_limit;};

    double get_relabeling_time() const {return time_relabel;}
    value_type relabeling_count() const {return _total_relabeling;}
};


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
//    bool _measure_spanning_cluster_by_bond{false};
//    std::vector<CalculationFlags> _calculation_flags;

    value_type min_index; // minimum index = 0
    value_type max_index; // maximum index = length - 1


    // index sequence
    std::vector<Index> index_sequence;  // initialized once
//    std::vector<Index> randomized_index_sequence;
    std::vector<value_type> randomized_index;



    // every birthTime we create a cluster we assign an set_ID for them
    int _cluster_id{};
//    int _impurity_id{-2};   // id of the impure sites

    /// Cluster tracker
    // key      -> id of cluster
    // value    -> index of cluster

    value_type _index_last_modified_cluster{};  // id of the last modified cluster

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

    virtual void reset();



    /*************************************************
     * Flags
     ************************************************/
    bool periodicity() const {return _periodicity;}


    /***********************************************
     * Properties of Percolation class
     ***********************************************/
    std::string getSignature();


    void numberOfActiveSites() const {std::cout << "Number of active sites " << _number_of_occupied_sites << std::endl;}
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
    void subtract_entropy_for_bond(const std::set<value_type> &found_index_set);

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
    void relabel_sites(const Cluster&  clstr, int id);  // todo pass cluster as reference
    void relabel_sites_v4(Index root_a, const Cluster& clstr_b); // relative index is set accordingly
    void relabel_sites_v5(Index root_a, const Cluster& clstr_b); // relative index is set accordingly
    void relabel_sites_v6(Index root_a, const Cluster& clstr_b, int id); // relative index is set accordingly
    void relabel_bonds(const Cluster&  clstr, int id);  // todo


    /**********************************************
     * Information about current state of Class
     **********************************************/
    double numberOfOccupiedSite() const { return _number_of_occupied_sites;}
    double occupationProbability() const { return double(_number_of_occupied_sites)/maxSites();}
    double spanningProbability() const; // number of bonds in spanning cluster / total number of bonds (2*l*l - 2*l)
    double entropy(); // the shannon entropy
    double entropy_by_site(); // the shannon entropy. the full calculations. time consuming
    double entropy_by_bond(); // the shannon entropy. the full calculations. time consuming


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
    IndexRelative getRelativeIndex(Index root, Index site_new);

};


/******************************************************************************
 * Explosive site percolation in square lattice with sum rule and product rule
 */
class SitePercolationExplosive: public SitePercolation_ps_v9{

public:
    ~SitePercolationExplosive() = default;
    SitePercolationExplosive(value_type length);
    std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_explosive_";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        s += std::to_string(length());
        return s;
    }

};

/*******************************************************************************
 * Site Percolation Ballistic Deposition
 *
 * *************************************************************/
/************************************
 * Extended from SitePercolation_ps_v9
 */

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
class SitePercolationBallisticDeposition_v2: public SitePercolation_ps_v9{
protected:
    // elements of @indices_tmp will be erased if needed but not of @indices
    std::vector<value_type> indices;
    std::vector<value_type> indices_tmp;
public:
    static constexpr const char* signature = "SitePercolation_BallisticDeposition_v2";
    virtual ~SitePercolationBallisticDeposition_v2(){
        indices.clear();
        indices_tmp.clear();
    };
    SitePercolationBallisticDeposition_v2(value_type length, bool periodicity);

    virtual bool occupy();

    /************************************
     * Site selection methods
     */
    Index select_site(std::vector<Index> &sites, std::vector<BondIndex> &bonds);
    Index select_site_upto_1nn(std::vector<Index> &sites, std::vector<BondIndex> &bonds);
    Index select_site_upto_2nn(std::vector<Index> &sites, std::vector<BondIndex> &bonds);


    void reset(); // todo
    void initialize_indices();
//    void randomize_index();

    virtual std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_ballistic_deposition_";
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

    value_type placeSite_1nn_v2();


    /*********************************
     * occupy upto 2nd nearest neighbor.
     * If the randomly selected site is occupied then select one of the nearest neighor randomly
     * If it is also occupied, select the next neighbor in the direction of motion Else occupy it.
     * If the 2nd nearest neighbor in the direction of motion is also occupied then skip the rest of the steps
     *      and start the next iteration
     *
     */

    value_type placeSite_2nn_v1();

};

/***********
 * Only L1
 */
class SitePercolationBallisticDeposition_L1_v2: public SitePercolationBallisticDeposition_v2{
public:
    ~SitePercolationBallisticDeposition_L1_v2() = default;
    SitePercolationBallisticDeposition_L1_v2(value_type length, bool periodicity)
            :SitePercolationBallisticDeposition_v2(length, periodicity){}

    bool occupy() {
        // if no site is available then return false

        if(_number_of_occupied_sites == maxSites()){
            return false;
        }

        try {
//        value_type v = placeSite_1nn_v0(); // debugging version
            value_type v = placeSite_1nn_v2();
            _occuption_probability = occupationProbability(); // for super class
            return v != ULLONG_MAX;
        }catch (OccupiedNeighbor& on){
//        on.what();
//        cout << "line : " << __LINE__ << endl;
            return false;
        }

    }

    std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_ballistic_deposition_L1";
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
class SitePercolationBallisticDeposition_L2_v2: public SitePercolationBallisticDeposition_v2{
public:
    ~SitePercolationBallisticDeposition_L2_v2() = default;
    SitePercolationBallisticDeposition_L2_v2(value_type length, bool periodicity)
            :SitePercolationBallisticDeposition_v2(length, periodicity){}

    bool occupy() {
        // if no site is available then return false

        if(_number_of_occupied_sites == maxSites()){
            return false;
        }

        try {

//            value_type v = placeSite_2nn_v0();
            value_type v = placeSite_2nn_v1();
            _occuption_probability = occupationProbability(); // for super class

            return v != ULLONG_MAX;
        }catch (OccupiedNeighbor& on){
//        on.what();
//        cout << "line : " << __LINE__ << endl;
            return false;
        }

    }

    std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_ballistic_deposition_L2";
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
 * version 1
 *
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
 * 4. Occupation probability is calculated by bonds,
 *      i.e., number of active bonds divided by total number of bonds
 *
 * 5. Spanning is calculated by number of sites in a spanning clusters with periodicity turned off,
 *      i.e., number of sites in the spanning clusters divided by total number of sites
 */
class BondPercolation_pb_v1 : public SqLatticePercolation{

    // flags to manipulate method
    bool _periodicity{false};
    bool _measure_spanning_cluster_by_bond{false};


    // index sequence
    std::vector<value_type> indices;
    std::vector<BondIndex> index_sequence; // initialized only once
    std::vector<BondIndex> randomized_index_sequence;
    BondIndex _last_placed_bond;

    double _total_number_of_active_bonds{};


    /// Cluster tracker
    // key      -> id of cluster
    // value    -> index of cluster

    value_type _index_last_modified_cluster{};  // id of the last modified cluster


    std::vector<value_type> number_of_bonds_to_span;
    std::vector<value_type> number_of_sites_to_span;
    value_type _sites_required_to_min_span;


    // first value the cluster id and second value is the cluster index
    std::vector<std::pair<value_type, value_type>> _spanning_cluster_indices;


    value_type sites_in_cluster_with_size_greater_than_one{};

    std::vector<Index> _wrapping_indices;

public:
    static constexpr const char* signature = "BondPercolation_v1";
    ~BondPercolation_pb_v1() = default;

    explicit BondPercolation_pb_v1(value_type length, bool periodicity=true);

    std::string getSignature() {
        std::string s = "sq_lattice_bond_percolation_";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        return s;
    }

    void periodicityFlag(bool p){_periodicity=p;};

    void reset();


    double occupationProbability() const override;
    std::vector<double> spanningProbability() const ;

    double entropy();
    double entropy_slow();

    void calculate_spanning_probability();
    void calculate_spanning_probability_by_largest_cluster();

    void placeAllBonds(){
        while(occupy());
    }

    value_type placeBond_v0();
    value_type placeBond_v1();
    bool occupy();

    /**********************
     * Relabeling
     */
    void relabel_sites(const Cluster& clstr, int id);
    void relabel_bonds(const Cluster& clstr, int id);
    void relabel_bonds_v1(BondIndex site_a, const Cluster &clstr_b); // implemented on 17 Aug 2018
    void relabel_bonds(const std::vector<BondIndex> &sites, int id_a, int delta_x_ab, int delta_y_ab); // implemented on 17 Aug 2018

    void relabel_v1(BondIndex last_bond, const Cluster& clstr_b); // relative index is set accordingly. implemented on 17 Aug 2018
    void relabel_sites(const std::vector<Index> &sites, int id_a, int delta_x_ab, int delta_y_ab) ;

    void relabel_sites_relative(BondIndex bond, const std::vector<Index> &sites);
    void relabel_new_sites_relative(const std::vector<Index> &sites, int id);

    // relabel sites and bonds in the cluster cluster
    void relabel_cluster(BondIndex bond, const std::vector<Index>& sites);
    void relabel_cluster(BondIndex bond, const Cluster& clstr_b, size_t bond_pos, size_t site_pos);

    void numberOfActiveSites() const {std::cout << "Number of active sites " << _total_number_of_active_bonds << std::endl;}
    double activeSites() const { return _total_number_of_active_bonds;}

    value_type count_number_of_active_site();

    void subtract_entropy_for_site(const std::set<value_type> &found_index_set);
    void add_entropy_for_site(value_type found_index);
    /***********************************
     * Spanning and Wrapping
     **********************************/
//    value_type numberOfSiteInSpanningClusters(std::vector<int> g_ids);
    value_type number_of_site_in_spanning_clusters(std::unordered_set<int> g_ids);

    bool detectSpanning();
    bool detectWrapping();

    IndexRelative getRelativeIndex(BondIndex root, BondIndex bond_new); // implemented on 16 Aug 2018
    IndexRelative getRelativeIndex(Index root, Index site_new);
    IndexRelative getRelativeIndex_v2(BondIndex root, BondIndex bond_new); // implemented on 17 Aug 2018

//    const std::vector<BondIndex>& wrapping_bonds() const { return  _wrapping_indices;}
    const std::vector<Index>& wrapping_indices() const { return  _wrapping_indices;}

    /*********************************
     * Printing Status
     * // todo declare these as constants
     ********************************/
    void periodicity_status();


    BondIndex lastPlacedBond() {
        if (_index_sequence_position == 0) {
            std::cerr << "Empty lattice : line " << __LINE__ << std::endl;
            return _last_placed_bond;
        }
        return _last_placed_bond;
    }


private:
    void initialize();
    void initialize_index_sequence();
    /**
    * Called only once when the object is constructed for the first time
    */
    void initialize_indices();

    void randomize();
    void randomize_v2();

//    std::vector<Index> get_Sites_for_bonds(BondIndex);


    value_type manage_clusters(
            const std::set<value_type> &found_index_set,
            std::vector<Index> &sites,
            BondIndex &bond
    );

    value_type manage_clusters(
            const std::set<value_type> &found_index_set,
            std::vector<Index> &sites,
            BondIndex &bond,
            int base_id
    );

    void connection_v1(BondIndex bond, std::vector<Index> &site_neighbor, std::vector<BondIndex> &bond_neighbor);
    void connection_v2(BondIndex bond, std::vector<Index> &site_neighbor, std::vector<BondIndex> &bond_neighbor);

    std::set<value_type> find_index_for_placing_new_bonds(const std::vector<BondIndex> &neighbors);
    int find_cluster_index_for_placing_new_bonds(
            const std::vector<BondIndex> &neighbors,
            std::set<value_type> &found_indices
    );

    int find_cluster_index_for_placing_new_bonds_v2(
            const std::vector<BondIndex> &neighbors,
            std::set<value_type> &found_indices
    );


    void connection_2_horizontal_no_periodicity(const BondIndex &bond, std::vector<Index> &site_neighbor,
                                                std::vector<BondIndex> &bond_neighbor,
                                                value_type next_column, value_type prev_column,
                                                value_type prev_row);

    void connection_2_vertical_no_peridicity(const BondIndex &bond, std::vector<Index> &site_neighbor,
                                             std::vector<BondIndex> &bond_neighbor, value_type prev_column, value_type prev_row,
                                             value_type next_row);

    void connection_2_periodic(const BondIndex &bond, std::vector<Index> &site_neighbor, std::vector<BondIndex> &bond_neighbor,
                               value_type next_column, value_type prev_column, value_type prev_row, value_type next_row);


};

/**
 * Universal class to access all types of percolation in square lattice.
 * More userfriendly.
 */
template <class T>
class Percolation{
public:
    enum class Types{
       SitePercolation,
       BondPercolation,
       SitePercolationBallisticDepositionL1,
       SitePercolationBallisticDepositionL2
    };

    ~Percolation() = default;
    T getClass(Types t, value_type length, bool periodicity=true);
};

#endif //SITEPERCOLATION_PERCOLATION_H

