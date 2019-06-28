//
// Created by shahnoor on 6/28/19.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_BOND_V1_H
#define SQLATTICEPERCOLATION_PERCOLATION_BOND_V1_H

#include "../percolation.h"


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
    void init() {randomize_v2();};
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

//    IndexRelative getRelativeIndex(BondIndex root, BondIndex bond_new); // implemented on 16 Aug 2018
//    IndexRelative getRelativeIndex(Index root, Index site_new);
//    IndexRelative getRelativeIndex_v2(BondIndex root, BondIndex bond_new); // implemented on 17 Aug 2018

//    const std::vector<BondIndex>& wrapping_bonds() const { return  _wrapping_indices;}
    const std::vector<Index>& wrapping_indices() const { return  _wrapping_indices;}

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

#endif //SQLATTICEPERCOLATION_PERCOLATION_BOND_V1_H
