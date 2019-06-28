//
// Created by shahnoor on 6/28/19.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_BOND_V2_H
#define SQLATTICEPERCOLATION_PERCOLATION_BOND_V2_H


#include "../percolation.h"


/**
 * Bond Percolation by Placing Bonds
 *
 * version 2
 *
 * First it randomizes the bond index list then use it.
 * Paradigm Shift:
 * 1. Does not delete cluster only makes it empty so that index and id remains the same.
 *    This way Searching for index of the cluster using id can be omitted.
 * 2. In bond percolation, initially all sites are present and they all form a cluster of unit size.
 *    This is maintained here unlike Version 1
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
 * 6. Unweighted relabeling is ommited in this version ??
 *
 * 7. Runtime is significantly improved. For example, if L=200 program will take ??? min to place all sites.
 *
 * 8. Unnecessary methods of previous version is eliminated
 *
 * 9. Checking spanning by keeping track of boundary sites is implemented
 *
 * 10. last modified cluster id can be obtained from @var _last_placed_site
 *
 *
 */
class BondPercolation_pb_v2 : public SqLatticePercolation{

    // flags to manipulate method
    bool _periodicity{false};

    // index sequence
    std::vector<BondIndex> bond_index_sequence; // initialized only once
    std::vector<Index> site_index_sequence;
    BondIndex _last_placed_bond;

    double _total_number_of_active_bonds{};

    std::vector<value_type> number_of_bonds_to_span;
    std::vector<value_type> number_of_sites_to_span;

    value_type sites_in_cluster_with_size_greater_than_one{};
    std::vector<Index> _wrapping_indices;

public:
    static constexpr const char* signature = "BondPercolation_v1";
    ~BondPercolation_pb_v2() = default;

    explicit BondPercolation_pb_v2(value_type length, bool periodicity=true);
    void init(bool random_seed=false); // some methods that are directly affected by the seed of random number generator
    // + methods that requires a bit computaion (so not in constructor)
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


    value_type placeBond_v1();
    bool occupy();

    /**********************
     * Relabeling
     */
    void relabel_sites(const Cluster_v3& clstr, int id);
    void relabel_bonds(const Cluster_v3& clstr, int id);
    void relabel_bonds_v1(BondIndex site_a, const Cluster_v3 &clstr_b); // implemented on 17 Aug 2018
    void relabel_bonds(const std::vector<BondIndex> &sites, int id_a, int delta_x_ab, int delta_y_ab); // implemented on 17 Aug 2018

    void relabel_v1(BondIndex last_bond, const Cluster_v3& clstr_b); // relative index is set accordingly. implemented on 17 Aug 2018
    void relabel_sites(const std::vector<Index> &sites, int id_a, int delta_x_ab, int delta_y_ab) ;

    void relabel_sites_relative(BondIndex bond, const std::vector<Index> &sites);
    void relabel_new_sites_relative(const std::vector<Index> &sites, int id);

    // relabel site_index_sequence and bonds in the cluster cluster
    void relabel_cluster(BondIndex bond, const std::vector<Index>& sites);
    void relabel_cluster(BondIndex bond, const Cluster_v3& clstr_b, size_t bond_pos, size_t site_pos);

    void numberOfActiveSites() const {std::cout << "Number of active site_index_sequence " << _total_number_of_active_bonds << std::endl;}
    double activeSites() const { return _total_number_of_active_bonds;}

    value_type count_number_of_active_site();

    void subtract_entropy_for_site(const std::set<value_type> &found_index_set, int base=-1);
    void add_entropy_for_site(value_type found_index);
    void track_numberOfBondsInLargestCluster();
    void track_numberOfSitesInLargestCluster();
    /***********************************
     * Spanning and Wrapping
     **********************************/
//    value_type numberOfSiteInSpanningClusters(std::vector<int> g_ids);
    value_type number_of_site_in_spanning_clusters(std::unordered_set<int> g_ids);

    bool detectSpanning();
    bool detectWrapping();
    bool detect_wrapping_v1();
//    bool detect_wrapping_v2();

    IndexRelative getRelativeIndex(BondIndex root, BondIndex bond_new); // implemented on 16 Aug 2018
    IndexRelative getRelativeIndex(Index root, Index site_new);
    IndexRelative getRelativeIndex_v2(BondIndex root, BondIndex bond_new); // implemented on 17 Aug 2018

//    const std::vector<BondIndex>& wrapping_bonds() const { return  _wrapping_indices;}
    const std::vector<Index> wrapping_indices() const { return  _wrapping_indices;}

    value_type numberOfBondsInTheWrappingClusters();
    value_type numberOfBondsInTheLargestCluster();
    value_type numberOfSitesInTheLargestCluster();
    value_type numberOfSitesInTheWrappingClusters();

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

    const std::vector<double> clusterSizeDistribution() const ;

private:
    void initialize();
    void initialize_index_sequence();
    /**
    * Called only once when the object is constructed for the first time
    */
    void initialize_cluster();

    void randomize_v2();


    value_type manage_clusters_v3(
            const std::set<value_type> &found_index_set,
            std::vector<Index> &sites,
            BondIndex &bond,
            int base_id
    );

};



#endif //SQLATTICEPERCOLATION_PERCOLATION_BOND_V2_H
