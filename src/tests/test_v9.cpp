//
// Created by shahnoor on 11/10/19.
//

#include "test_v9.h"
#include "../percolation/site/percolation_site_v9.h"

using namespace std;


void test_v9(int argc, char **argv) {
    test_v9_lattice(argc, argv);
    test_v9_number_of_clusters_per_site(argc, argv);
}

void test_v9_lattice(int argc, char **argv) {

    int length = atoi(argv[1]);
    SitePercolation_ps_v9 sp(length, true);
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


void test_v9_number_of_clusters_per_site(int argc, char **argv) {

    int length = atoi(argv[1]);
    SitePercolation_ps_v9 sp(length, true);
    sp.setRandomState(0, false);
    sp.init();

//    sp.reset();

    sp.view();
    sp.viewSiteByID();

    while(sp.occupy()) {

        sp.viewSiteByID();
        sp.viewCluster();
//        sp.viewClusterExtended();
        cout << "number of clusters = " << sp.numberOfClusters() << endl;
    }
}
