//
// Created by shahnoor on 12/30/19.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_BOND_V12_H
#define SQLATTICEPERCOLATION_PERCOLATION_BOND_V12_H


#include "../percolation_v12.h"

/**
 *
 *
 *
 */
class BondPercolation_pb_v12 : public Percolation_v12{
    size_t index_counter{};
    int id_last_bond{};
    std::vector<int> _wrapping_site_ids;
public:
    BondPercolation_pb_v12() = default;
    explicit BondPercolation_pb_v12(int length);
    void reset();
    std::string getSignature();
    std::string getClassName(){return "BondPercolation_pb_v12";};
    void init();

    bool occupy();

    IndexRelative  relabel_new_site(int id_current);
    void relabel_v3(int id_current_bond,  std::vector<int>& id_sites_connected, Cluster_v12 &clstr_b);

    IndexRelative getRelativeIndexDX(Index root, Index site_new);
    IndexRelative getRelativeIndexDX_v2(Index root, Index site_new);

    int lastBond() const {return id_last_bond;}
    Link lastBondIndex()  {return _lattice.getBond(id_last_bond).getIndex();}
    int lastBondGroupID() { return  _lattice.getBond(id_last_bond).get_groupID();}

    int sign(int a);

    bool detectWrapping();
    int wrappingSite_id();
    Index wrappingSite();

    // numerical properties
    long double entropy_v1();
    long double entropy_v2();

    long double entropy() override;
    double occupationProbability() override;
    size_t wrappingClusterSize() override;
    size_t maxIterationLimit()override {return _max_iteration_limit;}
    size_t numberOfSitesInTheWrappingClusters() override ;
    size_t numberOfBondsInTheWrappingClusters() override ;

    void manageCluster();

    void track_clusters(int root);
    void add_entropy(int root);
    void subtract_entropy(const std::set<int>& gids);
};


#endif //SQLATTICEPERCOLATION_PERCOLATION_BOND_V12_H
