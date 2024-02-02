//
// Created by shahnoor on 2/15/21.
//

#include "test_v13.h"
#include "../v13/lattice_v13.h"
#include "../v13/percolation_v13.h"
#include "../v13/percolation_v13_bond.h"
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

void test_detect_wrapping() {
//# take arguments from commandline
    cout << "Got here " << __LINE__ << endl;
    auto sq_lattice_p = SitePercolationL0_v13(6, 18, false);
    cout << "Got here " << __LINE__ << endl;
//    sq_lattice_p.setRandomState(0, true);

//# sq_lattice_p.viewLattice(3)
//# sq_lattice_p.viewCluster()
    int i = 0;
    cout << "Got here " << __LINE__ << endl;
    while (sq_lattice_p.place_one_site()) {
        double H1 = sq_lattice_p.entropy_v1();
        double H2 = sq_lattice_p.entropy_v2();
        cout << "p= " << sq_lattice_p.occupation_prob() <<
             " entropy_v1 " << H1 <<
             " entropy_v2 " << H2 <<
             " order " << sq_lattice_p.order_param_wrapping() << endl;
//        sq_lattice_p.viewLattice(3);
//        sq_lattice_p.viewLattice(4);
//        sq_lattice_p.lattice_ref.print_bonds();
#ifdef UNIT_TEST
    if (abs(H1-H2) > 1e-6){
        cout << "Error : entropy not equal " << __LINE__ << endl;
        exit(-1);
    }

#endif
        sq_lattice_p.viewCluster(1);
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
    sq_lattice_p.viewLattice(3);
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
    //    sq_lattice_p.viewLattice(3);
       sq_lattice_p.viewLattice(1);
//        sq_lattice_p.lattice_ref.print_bonds();
#ifdef UNIT_TEST
        if (abs(H1-H2) > 1e-6){
            cout << "Error : entropy not equal " << __LINE__ << endl;
            exit(-1);
        }

#endif
    //    sq_lattice_p.viewCluster(0);
    //    sq_lattice_p.viewCluster(2, true);
        iteration_count += 1;
        sq_lattice_p.detect_wrapping();
    //    if (iteration_count >40) break;
        // cout << "Press Enter to continue..." << endl;
        // cin.get();
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
    percolation_seed_length_pairL1(5, 1);

    // std::random_device _rd;
    // value_type _random_state = _rd();
    // std::mt19937 _random_engine;

    // _random_engine.seed(_random_state); // seeding


    // int max_length = 100;
    // for (int i=0; i < 100; ++ i) {

    //     int length = 5 +_random_engine() % (max_length - 5);
    //     value_type seed = _random_engine();
    //     cout << "Test run " << i <<  " => Length, seed " << length << ", " << seed << " ****************" << endl;
    //     percolation_seed_length_pairL1(length, seed);
    //     cout << "Test run " << i << " status : success" << endl;
    // }

}

void test_percolation_L0() {
    std::random_device _rd;
    value_type _random_state = _rd();
    std::mt19937 _random_engine;

    _random_engine.seed(_random_state); // seeding


    int max_length = 6;
    for (int i=0; i < 1; ++ i) {

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
        // if (i > 16) break;

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

void test_percolation_bond(){
    int length=5;
    value_type seed=0;
    auto sq_lattice_p = BondPercolation_v13(length, seed, false);
//    sq_lattice_p.setRandomState(0, true);
    // sq_lattice_p.viewLattice(3);
    // sq_lattice_p.viewCluster();
    int i = 0;
    while (sq_lattice_p.place_one_bond()) {
        // double H1 = sq_lattice_p.entropy_v1();
        // double H2 = sq_lattice_p.entropy_v2();
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
        // sq_lattice_p.detect_wrapping();

        // sq_lattice_p.test_cluster();
        // if (i > 16) break;

    }
    // sq_lattice_p.test_lattice();
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
//    test_detect_wrapping();
//    test_percolation_L0();
//    test_percolation_L1();
//    percolation_seed_length_pairL1(6, 455251785);
//    test_reset();
    // run_ensemble_v13(argc, argv);
    // percolation_fractalDimension_by_spanning_site_v13(1000, 100, 500, 100);
    // quantities_at_pc_v13(argc, argv);

    // 2024.02.01
    test_percolation_bond();
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
    run_v13_rsbd_L0(length, ensemble_size);
//    run_v13_rsbd_L1(length, ensemble_size);
//    run_v13_rsbd<SitePercolationL0_v13>(length, ensemble_size);
//    run_v13_rsbd<SitePercolationL1_v13>(length, ensemble_size);

}


/**
 * data for log(Mass) vs log(Length) curve
 * for spanning sites only.
 * To measure the fractal dimension df or d_f from log(M_pc) vs log(L) graph.
 * M_pc = number of sites in wrapping cluster at pc.
 * 
 */
void percolation_fractalDimension_by_spanning_site_v13(
        value_type ensemble_size, value_type L_start, value_type l_end, value_type delta_L
)
{
    double M{}, tc{};
    
    ofstream fout("fractal_dimension_-"+currentTime()+".txt");
    fout << "#Measuring fractal dimension by only spanning site_index_sequence" << endl;
    fout << "#MSitePercolationL1_v13" << endl;
    fout << "#<Length>\t<Mass>\t<tc>" << endl;
    clock_t t, t_outer;
    for(value_type len{L_start}; len <= l_end ; len += delta_L){
        cout << "Length " << len << " : " << endl;
        SitePercolationL1_v13 sp(len);
        M = 0;
        tc=0;
        t_outer = clock();
        for(value_type i{} ; i != ensemble_size ; ++i){
            t = clock();
            sp.reset();
            while(!sp.detect_wrapping()){
                sp.place_one_site();
            }
            tc += sp.occupation_prob();
            // M += sp.numberOfSitesInTheSpanningClusters();
            M += sp.get_wrapping_cluster_site_count_at_pc();
            cout << "\t\tIteration " << i << " . birthTime " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
        }
        fout << len << '\t' << M / double(ensemble_size) << '\t' << tc/ensemble_size << endl; // writing data to file
        cout << "Length " << len << " . Time taken "
             << getFormattedTime((clock() - t_outer) / double(CLOCKS_PER_SEC)) << endl;
    }
    fout.close();
}


/**
 * data for log(Mass) vs log(Length) curve
 * for spanning sites only.
 * To measure the fractal dimension df or d_f from log(M_pc) vs log(L) graph.
 * M_pc = number of sites in wrapping cluster at pc.
 * 
 * Cluster size distribution for exponent tau?
 * at pc
 * 
 * s = 1,2,3,.... cluster sizes
 * n(s) = number of cluster of size `s`
 * s vs n(s) gives tau. Ofcourse you need to normalize n(s) by dividing it by the total number of clusters
 */
void quantities_at_pc_v13(
        int argc, char **argv
)
{
    cout << "length        = argv[1]" << endl;
    cout << "ensemble size = argv[2]" << endl;
//    cout << "rsbd l        = argv[3]" << endl;
    if (argc <=2 ) cerr << "Not enough arguments " << endl;
    int length = stoi(argv[1]);
    int ensemble_size = stoi(argv[2]);

    double M{}, tc{};
    std::vector<double> tmp, cluster_sz_dist;
    std::string tm = getCurrentTime();
    

    
    clock_t t, t_outer;

    SitePercolationL1_v13 sp(length);
    auto signature = sp.get_signature();
    std::string extension = signature + "_L" + std::to_string(length) + '_' + tm + ".txt";

    M = 0;
    tc=0;
    t_outer = clock();
    for(value_type i{} ; i != ensemble_size ; ++i){
        t = clock();
        // SitePercolationL1_v13 sp(length);
        sp.reset();
        auto flag = true;
        while(flag){
            sp.place_one_site();
            // cout << "Got to Line " << __LINE__ << endl;
            flag = !sp.detect_wrapping();
            // cout << "Got to Line " << __LINE__ << endl;

        }
        tmp = sp.clusterSizeDistribution(); // TODO
        if (tmp.size() > cluster_sz_dist.size()){
            cluster_sz_dist.resize(tmp.size());
        }
        for(int i=0; i < tmp.size(); ++i){
            cluster_sz_dist[i] += tmp[i];
        }
        tc += sp.occupation_prob();
        // M += sp.numberOfSitesInTheSpanningClusters();
        M += sp.get_wrapping_cluster_site_count_at_pc();
        cout << "\t\tIteration " << i << " . time taken " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }
    std::ostringstream header_info;
    header_info << "{"
                << R"("length":)" << length
                << R"(,"ensemble_size":)" << ensemble_size
                << R"(,"En":)" << ensemble_size
                << R"(,"random_seed":)" << sp.getRandomState()
                << R"(,"signature":")" << sp.get_signature() << "\""
                << R"(,"datetime":")" << tm << "\""
                //                << R"(,"cols":)" << R"(["S","n_S"])"
                << "}";

    ofstream fout("fractal_dimension_"+extension);
    fout << "#" << header_info.str() << endl;
    fout << "#Measuring fractal dimension by only spanning site_index_sequence" << endl;
    fout << "#" << signature << endl;
    fout << "#<Length>\t<Mass>\t<tc>" << endl;
    fout << length << '\t' << M / double(ensemble_size) << '\t' << tc/ensemble_size << endl; // writing data to file
    fout.close();


    ofstream fout2("cluster_size_dist_pc_"+extension);
    fout2 << "#" << header_info.str() << endl;
    fout2 << "#Measuring fractal dimension by only spanning site_index_sequence" << endl;
    fout2 << "#" << signature << endl;
    fout2 << "#<s>,<n(s)>" << endl;
    for(int i=1; i < cluster_sz_dist.size(); ++i){
        if (cluster_sz_dist[i] == 0) continue;
        fout2 << i << "\t" << cluster_sz_dist[i]/(ensemble_size) << endl;
    }
    fout2.close();
    
    cout << "Length " << length << " . Time taken "
            << getFormattedTime((clock() - t_outer) / double(CLOCKS_PER_SEC)) << endl;
    
    
}
/**
 * @brief Plotting the state of lattice at pc
 * showing sites of same cluster in same color for visual appeal
 * 
 * @param length 
 */

void visualize_lattice(int length){
    cout << "length        = argv[1]   : " << length << endl;


    SitePercolationL0_v13 sp(length);
    auto signature = sp.get_signature();
    std::string tm = getCurrentTime();
    std::string extension = signature + "_L" + std::to_string(length) + '_' + tm + ".txt";

    while(!sp.detect_wrapping()){
        sp.place_one_site();
    }
    
    std::ostringstream header_info;
    header_info << "{"
                << R"("length":)" << length
                << R"(,"random_seed":)" << sp.getRandomState()
                << R"(,"signature":")" << sp.get_signature() << "\""
                << R"(,"datetime":")" << tm << "\""
                //                << R"(,"cols":)" << R"(["S","n_S"])"
                << "}";

    ofstream fout("sq_lattice_site_percolation_visual_"+extension);
    fout << "#" << header_info.str() << endl;
    fout << "#Only non empty site info are recorded" << endl;
    fout << "#x and y coordinate and the color of the site" << endl;
    fout << "#" << signature << endl;
    fout << "#<x>\t<y>\t<color>" << endl;

    int max_bonds = length*length;
    cout << "Writing data to file" << endl;
    for (int i=0; i < max_bonds; ++i){
        auto site = sp.lattice_ref.get_site_by_id(i);
        if(site.is_occupied()){
            fout << site.get_index().row() << "\t" << site.get_index().col() << "\t" << site.get_gid() << endl;
        }
    }

    fout.close();

    
}