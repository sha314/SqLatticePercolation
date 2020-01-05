//
// Created by shahnoor on 12/30/19.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_SITE_V12_H
#define SQLATTICEPERCOLATION_PERCOLATION_SITE_V12_H


#include "../percolation_v12.h"

/**
 *
 *
 *
 *
 */
class SitePercolation_ps_v12 :  public Percolation_v12{
    size_t index_counter{};
    int id_last_site{};
    std::vector<int> _wrapping_site_ids;
public:
    SitePercolation_ps_v12() = default;
    explicit SitePercolation_ps_v12(int length);
    void reset();
    std::string getSignature();
    std::string getClassName(){return "SitePercolation_ps_v12";};
    void init();

    bool occupy();

    IndexRelative  relabel_new_site(int id_current);
    void relabel(Cluster_v12& clstr, int id_current);
    void relabel_v2(Cluster_v12& clstr, int id_current, IndexRelative dx_dy);
    void relabel_v3(int id_current_a, std::vector<Index>& neighbors_a, Cluster_v12 &clstr_b);

    IndexRelative getRelativeIndexDX(Index root, Index site_new);
    IndexRelative getRelativeIndexDX_v2(Index root, Index site_new);

    int lastSite() const {return id_last_site;}
    Index lastSiteIndex()  {return _lattice.getSite(id_last_site).get_index();}
    int lastSiteGroupID() { return  _lattice.getSite(id_last_site).get_groupID();}

    int sign(int a);

    bool detectWrapping();
    int wrappingSite_id();
    Index wrappingSite();

    // numerical properties
    long double entropy_v1_site();
    long double entropy_v1_bond();
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


#endif //SQLATTICEPERCOLATION_PERCOLATION_SITE_V12_H
