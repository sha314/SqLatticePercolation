//
// Created by shahnoor on 11/10/19.
//
#include "test_v10.h"
#include "../percolation/site/percolation_site_v10.h"
#include <iostream>

using namespace std;


void test_v10(int argc, char **argv) {
//    test_v9_lattice(argc, argv);
    test_v10_number_of_clusters_per_site(argc, argv);
}

void test_v10_lattice(int argc, char **argv) {

    int length = atoi(argv[1]);
    SitePercolation_ps_v10 sp(length, true);
    sp.setRandomState(0, false);
    sp.init();

//    sp.reset();

    sp.view();
    sp.viewSiteByID();

    while(sp.occupy()) {

        sp.viewSiteByID();
        sp.viewCluster();
//        sp.viewClusterExtended();
    }
}


void test_v10_number_of_clusters_per_site(int argc, char **argv) {

    int length = atoi(argv[1]);
    int ensemble_size = atoi(argv[2]);

    SitePercolation_ps_v10 sp(length, true);
    sp.setRandomState(0, false);
    sp.init();

//    sp.reset();

//    sp.view();
//    sp.viewSiteByID();
    std::vector<size_t> number_of_clusters(sp.maxSites());
    size_t k{};
    for(int e{1}; e <= ensemble_size; ++e) {
        cout << "Iteration " << e << endl;
        sp.reset();
        k=0;
        number_of_clusters[k++] += sp.numberOfClusters();
        while (sp.occupy()) {
//            cout << "***********" << endl;
//        sp.viewSiteByID();
//            sp.viewCluster();
//        sp.viewClusterExtended();
//            cout << "number of clusters = " << sp.numberOfClusters() << endl;
            number_of_clusters[k++] += sp.numberOfClusters();
        }

//        sp.viewSiteByID();
//        sp.viewCluster();
//        sp.viewClusterExtended();
    }

    double size = sp.maxSites();
    ofstream fout("cluster_per_site");
    fout << "#p=occupation probability" << endl;
    fout << "#n(p)=number of clusters per site at p" << endl;
    fout << "#n(p)=(number of clusters at p) / (L*L)" << endl;
    fout << "#<p><n(p)>" << endl;
    for(size_t m{}; m < number_of_clusters.size(); ++m){
        fout << m/size << "\t" << number_of_clusters[m]/(size * ensemble_size) << endl;
    }
    fout.close();
}
