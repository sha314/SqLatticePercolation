//
// Created by shahnoor on 2/15/21.
//

#include "test_v13.h"
#include "../v13/lattice_v13.h"
#include "../v13/percolation_v13.h"
#include "../v13/run_v13.h"
#include "../v13/percolation_rsbd_v13.h"
#include <iostream>

void test_reset();

using namespace std;

void percolation_seed_length_pair(int length, value_type seed);

void test_lattice(int argc, char** argv) {
    cout << "length        = argv[1]" << endl;
    cout << "ensemble size = argv[2]" << endl;
//    cout << "rsbd l        = argv[3]" << endl;
    if (argc <=2 ) cerr << "Not enough arguments " << endl;
    int length = stoi(argv[1]);
    int ensemble_size = stoi(argv[2]);

    auto lattice = Lattice_v13(length);
//    lattice.view(0);
//    lattice.view(1);
//    lattice.view(2);
//    lattice.view_relative_index();
//    lattice.view_site_gids();
//    lattice.test_id_index_conversion();
//    lattice.print_site_bond_list();
}

void test_detect_wrapping(int argc, char **argv) {
    cout << "length        = argv[1]" << endl;
    cout << "ensemble size = argv[2]" << endl;
//    cout << "rsbd l        = argv[3]" << endl;
    if (argc <=2 ) cerr << "Not enough arguments " << endl;
    int length = stoi(argv[1]);
    int ensemble_size = stoi(argv[2]);
//# take arguments from commandline
//    cout << "Got here " << __LINE__ << endl;
    auto sq_lattice_p = SitePercolationL0_v13(length, 18, false);
//    cout << "Got here " << __LINE__ << endl;
//    sq_lattice_p.setRandomState(0, true);

//# sq_lattice_p.viewLattice(3)
//# sq_lattice_p.viewCluster()
    int i = 0;
//    cout << "Got here " << __LINE__ << endl;
    while (sq_lattice_p.place_one_site()) {
        double H1 = sq_lattice_p.entropy_v1();
        double H2 = sq_lattice_p.entropy_v2();
//        cout << "p= " << sq_lattice_p.occupation_prob() <<
//             " entropy_v1 " << H1 <<
//             " entropy_v2 " << H2 <<
//             " order " << sq_lattice_p.order_param_wrapping() << endl;
//        sq_lattice_p.viewLattice(3);
//        sq_lattice_p.viewLattice(4);
//        sq_lattice_p.lattice_ref.print_bonds();
#ifdef UNIT_TEST
    if (abs(H1-H2) > 1e-6){
        cout << "Error : entropy not equal " << __LINE__ << endl;
        exit(-1);
    }

#endif
//        sq_lattice_p.viewCluster(1);
        i += 1;
        if (sq_lattice_p.detect_wrapping()) {
            cout << "p= " << sq_lattice_p.occupation_prob() <<
                 " entropy_v1 " << H1 <<
                 " entropy_v2 " << H2 <<
                 " order " << sq_lattice_p.order_param_wrapping() << endl;
            cout << "Wrapping detected ***************** <<<" << endl;
            break;
        }

//        if (i > 8) break;

    }
//    sq_lattice_p.viewLattice(3);
//    sq_lattice_p.viewLattice(4);
//    sq_lattice_p.viewLattice(1);
//    sq_lattice_p.viewCluster();
}
void percolation_seed_length_pairL1(int length, value_type seed) {

    auto sq_lattice_p = SitePercolationL1_v13(length, seed, false);
//    sq_lattice_p.setRandomState(0, true);

//# sq_lattice_p.viewLattice(3)
//# sq_lattice_p.viewCluster()
    int iteration_count = 0;
    double pp = sq_lattice_p.occupation_prob();
//    cout << "p= " << pp << endl;
    while (sq_lattice_p.place_one_site()) {
        double H1 = sq_lattice_p.entropy_v1();
        double H2 = sq_lattice_p.entropy_v2();
//        cout << "p= " << sq_lattice_p.occupation_prob() <<
//             " entropy_v1 " << H1 <<
//             " entropy_v2 " << H2 <<
//             " order " << sq_lattice_p.order_param_wrapping() << endl;
//        sq_lattice_p.viewLattice(3);
//        sq_lattice_p.viewLattice(4);
//        sq_lattice_p.lattice_ref.print_bonds();
#ifdef UNIT_TEST
        if (abs(H1-H2) > 1e-6){
            cout << "Error : entropy not equal " << __LINE__ << endl;
            exit(-1);
        }

#endif
//        sq_lattice_p.viewCluster(1);
        iteration_count += 1;
        sq_lattice_p.detect_wrapping();
//        if (i > 8) break;

    }

#ifdef UNIT_TEST
    if(iteration_count <= length*length){
        cout << "Not full iteration limit reached : " << __LINE__ << endl;
        cout << "Number of iteration " << iteration_count << endl;
        exit(-1);
    }
//    sq_lattice_p.viewCluster(0);
    double P2 = sq_lattice_p.order_param_wrapping();
    double P1 = sq_lattice_p.order_param_largest_clstr();

    if (abs(P1-1.0) > 1e-6){
        cout << "Error : order parameter P1 not equal to 1.0. line " << __LINE__ << endl;
        cout << "P1 = " << P1 << endl;
        exit(-1);
    }

    if (abs(P2-1.0) > 1e-6){
        cout << "Error : order parameter P2 not equal to 1.0. line " << __LINE__ << endl;
        cout << "P2 = " << P2 << endl;
        exit(-1);
    }
    double p = sq_lattice_p.occupation_prob();

    if (abs(p-1.0) > 1e-6){
        cout << "Error : occupation_prob p not equal to 1.0. line " << __LINE__ << endl;
        cout << "p = " << p << endl;
        exit(-1);
    }

#endif
//    sq_lattice_p.viewLattice(3);
//    sq_lattice_p.viewLattice(4);
//    sq_lattice_p.viewLattice(1);
//    sq_lattice_p.viewCluster();
}
void test_percolation_L1() {
    std::random_device _rd;
    value_type _random_state = _rd();
    std::mt19937 _random_engine;

    _random_engine.seed(_random_state); // seeding


    int max_length = 100;
    for (int i=0; i < 100; ++ i) {

        int length = 5 +_random_engine() % (max_length - 5);
        value_type seed = _random_engine();
        cout << "Test run " << i <<  " => Length, seed " << length << ", " << seed << " ****************" << endl;
        percolation_seed_length_pairL1(length, seed);
        cout << "Test run " << i << " status : success" << endl;
    }

}

void test_percolation_L0() {
    std::random_device _rd;
    value_type _random_state = _rd();
    std::mt19937 _random_engine;

    _random_engine.seed(_random_state); // seeding


    int max_length = 6;
    for (int i=0; i < 100; ++ i) {

        int length = 5 + _random_engine() % max_length;
        value_type seed = _random_engine();
        cout << "Test run " << i <<  " => Length, seed " << length << ", " << seed << " ****************" << endl;
        percolation_seed_length_pair(length, seed);
        cout << "Test run " << i << " status : success" << endl;
    }

}

void percolation_seed_length_pair(int length, value_type seed) {

    auto sq_lattice_p = SitePercolationL0_v13(length, seed, false);
//    sq_lattice_p.setRandomState(0, true);

//# sq_lattice_p.viewLattice(3)
//# sq_lattice_p.viewCluster()
    int i = 0;
    while (sq_lattice_p.place_one_site()) {
        double H1 = sq_lattice_p.entropy_v1();
        double H2 = sq_lattice_p.entropy_v2();
//        cout << "p= " << sq_lattice_p.occupation_prob() <<
//             " entropy_v1 " << H1 <<
//             " entropy_v2 " << H2 <<
//             " order " << sq_lattice_p.order_param_wrapping() << endl;
//        sq_lattice_p.viewLattice(3);
//        sq_lattice_p.viewLattice(4);
//        sq_lattice_p.lattice_ref.print_bonds();
#ifdef UNIT_TEST
        if (abs(H1-H2) > 1e-6){
            cout << "Error : entropy not equal " << __LINE__ << endl;
            exit(-1);
        }

#endif
//        sq_lattice_p.viewCluster(1);
        i += 1;
        sq_lattice_p.detect_wrapping();

        sq_lattice_p.test_cluster();
        if (i > 16) break;

    }
    sq_lattice_p.test_lattice();
#ifdef UNIT_TEST
    double P2 = sq_lattice_p.order_param_wrapping();
    double P1 = sq_lattice_p.order_param_largest_clstr();

    if (abs(P1-1.0) > 1e-6){
        cout << "Error : order parameter P1 not equal to 1.0. line " << __LINE__ << endl;
        exit(-1);
    }

    if (abs(P2-1.0) > 1e-6){
        cout << "Error : order parameter P2 not equal to 1.0. line " << __LINE__ << endl;
        cout << "P2 = " << P2 << endl;
        exit(-1);
    }
    double p = sq_lattice_p.occupation_prob();

    if (abs(p-1.0) > 1e-6){
        cout << "Error : occupation_prob p not equal to 1.0. line " << __LINE__ << endl;
        exit(-1);
    }

#endif
//    sq_lattice_p.viewLattice(3);
//    sq_lattice_p.viewLattice(4);
//    sq_lattice_p.viewLattice(1);
//    sq_lattice_p.viewCluster();
}

void test_v13(int argc, char **argv) {
//    test_lattice(argc, argv);
//    test_detect_wrapping(argc, argv);
//    test_percolation_L0();
//    test_percolation_L1();
//    percolation_seed_length_pairL1(6, 455251785);
//    test_reset();
    run_ensemble_v13(argc, argv);
}

void test_reset() {
//# take arguments from commandline
    auto sq_lattice_p = SitePercolationL0_v13(6, 18, false);
//    sq_lattice_p.setRandomState(0, true);

//    sq_lattice_p.viewLattice(3);
//    sq_lattice_p.viewCluster(1);
    sq_lattice_p.reset();
//    sq_lattice_p.viewLattice(3);
//    sq_lattice_p.viewCluster(1);

    int i = 0;
    while (sq_lattice_p.place_one_site()) {
        cout << "p= " << sq_lattice_p.occupation_prob() <<
             " entropy_v1 " << sq_lattice_p.entropy_v1() <<
             " entropy_v2 " << sq_lattice_p.entropy_v2() <<
             " order_largest " << sq_lattice_p.order_param_largest_clstr() <<
             " order_wrapping " << sq_lattice_p.order_param_wrapping() << endl;
//        sq_lattice_p.viewLattice(3);
//        sq_lattice_p.viewLattice(4);
//        sq_lattice_p.lattice_ref.print_bonds();
//        sq_lattice_p.viewCluster(1);
        i += 1;
        if (sq_lattice_p.detect_wrapping()) {
//            cout << "p= " << sq_lattice_p.occupation_prob() <<
//                 " entropy_v1 " << sq_lattice_p.entropy_v1() <<
//                 " entropy_v2 " << sq_lattice_p.entropy_v2() <<
//                 " order_largest " << sq_lattice_p.order_param_largest_clstr() <<
//                 " order_wrapping " << sq_lattice_p.order_param_wrapping() << endl;
//            cout << "Wrapping detected ***************** <<<" << endl;
//            break;
        }

//        if (i > 8) break;

    }
    sq_lattice_p.viewLattice(3);
//    sq_lattice_p.viewLattice(4);
//    sq_lattice_p.viewLattice(1);
//    sq_lattice_p.viewCluster();
}

void run_ensemble_v13(int argc, char **argv){
    cout << "length        = argv[1]" << endl;
    cout << "ensemble size = argv[2]" << endl;
//    cout << "rsbd l        = argv[3]" << endl;
    if (argc <=2 ) cerr << "Not enough arguments " << endl;
    int length = stoi(argv[1]);
    int ensemble_size = stoi(argv[2]);


//    test_v9_lattice(argc, argv);
//    test_v10_number_of_clusters_per_site(argc, argv);
//    test_v10_lattice(argc, argv);

//    run_v10_regular(length, ensemble_size);

//    int rsbd_l = stoi(argv[3]);
//    run_v13_rsbd_L0(length, ensemble_size);
    run_v13_rsbd_L1(length, ensemble_size);
//    run_v13_rsbd<SitePercolationL0_v13>(length, ensemble_size);
//    run_v13_rsbd<SitePercolationL1_v13>(length, ensemble_size);

}
