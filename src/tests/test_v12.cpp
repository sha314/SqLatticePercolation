//
// Created by shahnoor on 12/23/19.
//

#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include "test_v12.h"
#include "../lattice/lattice.h"
#include "../percolation/percolation_v12.h"
#include "../util/time_tracking.h"
#include "../percolation/bond/percolation_bond_v12.h"
#include "../percolation/site/percolation_site_v12.h"

using namespace std;

void test_v12(int argc, char **argv) {
    cout << "length        = argv[1]" << endl;
    cout << "ensemble size = argv[2]" << endl;
//    cout << "rsbd l        = argv[3]" << endl;

    int length = stoi(argv[1]);
    int ensemble_size = stoi(argv[2]);


//    test_v12_lattice(length);
//    test_v12_percolation_site(length);
//    test_v12_percolation_bond(length);

//    run_v12_regular_site(length, ensemble_size);
//    run_v12_regular_bond(length, ensemble_size);
    demarcationLine_v12(length);
}

void test_v12_lattice(int length) {
    Lattice_v12 lattice(length);

//    lattice.view_as_assigned();
    lattice.view_all_v2();
//    lattice.view_by_id_v2();
//    lattice.view_by_gid_v2();
//    lattice.view_by_relative_index_v2();
//    lattice.view_sites_list();
    lattice.view_bonds_list();

//    lattice.view_all_v1();
//    lattice.view_from_str_v1();

//    lattice.view_sites_v1();
//    lattice.view_sites_by_id_v1();
//    lattice.view_sites_by_gid_v1();
//    lattice.view_by_id_v1();
//    lattice.view_by_gid_v1();
//    lattice.view_sites_by_relative_index_v1();
//    lattice.view_by_relative_index_v1();

//    lattice.view_sites_list();
//    lattice.view_bonds_list();

//    cout << lattice.getSiteLocation(0) << endl;
//    cout << lattice.getSiteLocation(5) << endl;
//    cout << lattice.getSiteLocation(12) << endl;
//    cout << lattice.getSiteLocation(22) << endl;


//    cout << lattice.getBondLocation(11) << endl;
//    cout << lattice.getBondLocation(17) << endl;
//    cout << lattice.getBondLocation(22) << endl;
//    cout << lattice.getBondLocation(35) << endl;
//    cout << lattice.getBondLocation(51) << endl;
//
//
//    cout << lattice.getBond(0)<< endl;
//    cout << lattice.getBond(11)<< endl;
//    cout << lattice.getBond(17)<< endl;
//    cout << lattice.getBond(22)<< endl;
//    cout << lattice.getBond(35)<< endl;
}


void test_v12_percolation_site(int length) {
    SitePercolation_ps_v12 percolation(length);
    percolation.setRandomState(1549510744, false);
//    percolation.init({6,14,9,11,10});
    percolation.init();

//    percolation.viewCluster();
    percolation.viewLattice_by_id();
//    percolation.viewLattice_by_gid();

//    percolation.occupy();
    int i=0;
    cout << "entropy redefined   " << percolation.entropy_v1_bond() << endl;
    while(percolation.occupy()) {
        cout << "********************** last site id " << percolation.lastSite()
             << " index " << percolation.lastSiteIndex()
             << " group id " << percolation.lastSiteGroupID() << endl;
        cout << "entropy traditional " << percolation.entropy_v1_site() << endl;
        cout << "entropy traditional " << percolation.entropy_v2_site() << endl;

        cout << "entropy redefined   " << percolation.entropy_v1_bond() << endl;
        cout << "entropy redefined   " << percolation.entropy_v2_bond() << endl;
        cout << "entropy list        " << percolation.entropy_v3_list() << endl;
//    percolation.viewCluster();
//    percolation.viewLattice_by_id();
        percolation.viewLattice_by_gid();
//        percolation.viewLattice_sites_by_gid();
//    percolation.viewLattice();
//        percolation.viewLattice_by_relative_index();

//        if(percolation.detectWrapping()){
////            percolation.viewLattice_by_relative_index();
//            cout << "wrapping site " << percolation.wrappingSite() << endl;
//            cout << "wrapping site id " << percolation.wrappingSite_id() << endl;
//            cout << "p = " << percolation.occupationProbability() << endl;
//            break;
//        }

        ++i;
//        if(i >= 10) break;
    }


}


void test_v12_percolation_bond(int length) {
    BondPercolation_pb_v12 percolation(length);
    percolation.setRandomState(1549510744, false);
    percolation.init();

//    percolation.viewCluster();
//    percolation.viewLattice_by_id();
//    percolation.viewLattice_by_gid();

//    percolation.occupy();
    int i=0;
    while(percolation.occupy()) {
        cout << "*********" << i <<" th ************* last bond id " << percolation.lastBond()
             << " index " << percolation.lastBondIndex().get_string()
             << " group id " << percolation.lastBondGroupID() << endl;
//        percolation.viewCluster();
//    percolation.viewLattice_by_id();
//        percolation.viewLattice_by_gid();
//        percolation.viewLattice_sites_by_gid();
//    percolation.viewLattice();
//        percolation.viewLattice_by_relative_index();

        if(percolation.detectWrapping()){
            percolation.viewLattice_by_relative_index();
            cout << "wrapping site " << percolation.wrappingSite() << endl;
            cout << "wrapping site id " << percolation.wrappingSite_id() << endl;
            cout << "p = " << percolation.occupationProbability() << endl;
            break;
        }

        ++i;
//        if(i >= 8) break;
    }


}

/**
 * TODO : ensemble size 6, length 10, at critical point. seed 2562449578
 *
 * @param length
 * @param ensemble_size
 */
void run_v12_regular_site(int length, int ensemble_size) {

    std::cout << "length=" << length << " and ensemble_size=" << ensemble_size << std::endl;

    long length_squared = length*length;
    long twice_length_squared = 2 * length_squared;

    SitePercolation_ps_v12 lp(length);
//    SitePercolationBallisticDeposition_L2_v3 lp(length, true);
//    SitePercolationRSBD_L1_v10 lp(length, true);
//    SitePercolationRSBD_L2_v10 lp(length, true);

    lp.setRandomState(2562449578, true);
    lp.init();

//    lp.viewLattice();
//    lp.viewLatticeByID();
//    lp.viewLatticeExtended();

//    lp.viewCluster();
//    lp.viewRandomizedIndices();
//    lp.viewRemainingIndices();
//    lp.viewBondByID();
    // simulation starts here
    auto a = size_t(ensemble_size);
    std::vector<double> pcs(a), sites_pc(a), bonds_pc(a);
    value_type counter{};
    value_type limit = lp.maxIterationLimit();
    std::vector<double> entropy(limit), dHs(limit);
    std::vector<double> nob_wraping(limit),
            nob_largest(limit), dPs(limit);

    for(int i{} ; i != ensemble_size ; ++i){

        lp.reset();

        bool successful = false;
        auto t_start = std::chrono::system_clock::now();
        counter = 0;
        bool detect_wrapping{true};
        while (true){
//            lp.viewRemainingSites();
//            lp.viewIndices();
            successful = lp.occupy();

            if(successful) {
//                cout << counter << " th site ************########********** last " << lp.lastPlacedSite() << endl;
//                lp.viewRandomizedIndices();
//                lp.viewRemainingIndices();
//                lp.viewSiteByID();
//                lp.viewBondByID();
//                lp.viewCluster();
//                lp.ckeckCluster();
//                auto H1 = lp.entropy_v1_bond();
//                auto H2 = lp.entropy_v2_bond();
//                cout
//                     << H1
//                     << " vs " << H2
//                     << endl;
//                if(abs(H1 - H2) > 1e-12){
//                    cout.precision(numeric_limits<long double>::digits10);
//                    cout    << "iteration " << counter << " "
//                            << H1
//                            << " vs " << H2
//                            << endl;
//                    cout << "Entropy mismatched" << endl;
//                    exit(0);
//                }
                entropy[counter] += lp.entropy_v2_bond();
//                entropy[counter] += lp.entropy_v3_list();
                nob_wraping[counter] += lp.numberOfBondsInTheWrappingClusters();
                nob_largest[counter] += lp.numberOfBondsInTheLargestCluster();
//                lp.jump();
//                dHs[counter] += lp.jump_entropy();
//                dPs[counter] += lp.jump_wrapping_cluster();
                if(detect_wrapping && lp.detectWrapping()){

                    pcs[i] = lp.occupationProbability();
                    sites_pc[i] = lp.numberOfSitesInTheWrappingClusters();
                    bonds_pc[i] = lp.numberOfBondsInTheWrappingClusters();
//                    cout << "wrapping at "   << pcs[i] << ","
//                              << sites_pc[i] << ","
//                              << bonds_pc[i] << std::endl;


//                    std::vector<value_type> site, bond;
//                    lattice_percolation.get_cluster_info(site, bond);
//                    for(value_type j{}; j != site.size(); ++j){
//                        fout_s << site[j] << ',';
//                    }
//                    for(value_type j{}; j != bond.size(); ++j){
//                        fout_b << bond[j] <<',';
//                    }
//                    fout_s << std::endl;
//                    fout_b << std::endl;

                    detect_wrapping = false;
//                    break; // only for debugging
                }


                ++counter;
//                if(counter >= 3) break;
            }
            if(counter >= limit){ // twice_length_squared is the number of bonds
                break;
            }
        }
//        fout_jump << lattice_percolation.largestEntropyJump() << "," << lattice_percolation.largestEntropyJump_pc() << std::endl;

        auto t_end = std::chrono::system_clock::now();
        std::cout << "Iteration " << (i+1)
                  //                 << " . Thread " << std::this_thread::get_id()
                  << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << std::endl;

//        cout << "Relabeling time " << lattice_percolation.get_relabeling_time() << endl;
//        lp.viewClusterExtended();
    }
//    cout << "wrapping site id " << lp.wrappingSite_id();
//    cout << "  wrapping site index " << lp.wrappingSite() << endl;
//    lp.viewLattice_by_id();
//    lp.viewLattice_by_gid();
//    lp.viewLattice_sites_by_gid();
//    lp.viewLattice_by_relative_index();
//    lp.viewCluster();

//    fout_b.close();
//    fout_s.close();
//    fout_jump.close();

    std::string tm = getCurrentDateTime();
    const char delimiter = '\t';
    std::ostringstream header_info;
    header_info << "{"
                << R"("length":)" << length
                << R"(,"ensemble_size":)" << ensemble_size
                << R"(,"random_seed":)" << lp.getRandomState()
                << R"(,"signature":")" << lp.getSignature() << "\""
                << R"(,"datetime":")" << tm << "\""
                << R"(,"date":")" << getCurrentDate() << "\""
                << R"(,"time":")" << getCurrentTime() << "\""
                << R"(,"classname":")" << lp.getClassName() << "\""
                //                << R"(,"cols":)" << R"(["S","n_S"])"
                << "}";


    std::string extension = "_L" + std::to_string(length) + '_' + tm + ".txt";
//    std::string filename_s = lattice_percolation.getSignature() + "_cluster_by_site" + extension;
//    std::string filename_b = lattice_percolation.getSignature() + "_cluster_by_bond" + extension;
    std::string filename_critical = lp.getClassName() + "_critical" + extension;
//    std::string filename = lattice_percolation.getSignature() + "_entropy-jump" + extension;
    std::string filename_entropy_order_parameter = lp.getClassName() + "_entropy-order" + extension;

    std::string filename_specific_heat_susceptibility = lp.getClassName() + "_specific_heat-susceptibility" + extension;


//    std::ofstream fout_jump(filename);
//    // JSON formated header
//    fout_jump << '#' << header_info.str() << std::endl;
//    fout_jump << "#each line is an independent realization" << std::endl;
//    fout_jump << "#each line contains information about all clusters at critical point" << std::endl;
//    fout_jump << "#cluster size is measured by number of bonds in it" << std::endl;
//    fout_jump << "#<entropy jump>,<occupation probability at which the jump occurs>" << std::endl;
//
//    std::ofstream fout_s(filename_s);
//    // JSON formated header
//    fout_s << '#' << header_info.str() << std::endl;
//    fout_s << "#each line is an independent realization" << std::endl;
//    fout_s << "#each line contains information about all clusters at critical point" << std::endl;
//    fout_s << "#cluster size is measured by number of site_index_sequence in it" << std::endl;
//
//    std::ofstream fout_b(filename_b);
//    // JSON formated header
//    fout_b << '#' << header_info.str() << std::endl;
//    fout_b << "#each line is an independent realization" << std::endl;
//    fout_b << "#each line contains information about all clusters at critical point" << std::endl;
//    fout_b << "#cluster size is measured by number of bonds in it" << std::endl;

    std::ofstream fout_critical(filename_critical);
    fout_critical << '#' << header_info.str() << std::endl;
    fout_critical << "#data at critical occupation probability or pc" << std::endl;
    fout_critical << "#<pc>,<sites in wrapping cluster>,<bonds in wrapping cluster>" << std::endl;
    for(size_t i{}; i < ensemble_size; ++i){
        if(pcs[i] == 0) continue;
        fout_critical << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
                      << pcs[i] << delimiter
                      << sites_pc[i] << delimiter
                      << bonds_pc[i] << std::endl;
//        cout << pcs[i] << endl;
    }
    fout_critical.close();



    std::ofstream fout(filename_entropy_order_parameter);
    fout << '#' << header_info.str() << std::endl;
    fout << "#<p>,<H(p,L)>,<P1(p,L)>,<P2(p,L)>" << std::endl;
    fout << "#p = occupation probability" << std::endl;
    fout << "#H(p,L) = Entropy = sum( - u_i * log(u_i))" << std::endl;
    fout << "#P1(p,L) = Order parameter = (number of bonds in largest cluster) / (total number of bonds)" << std::endl;
    fout << "#P2(p,L) = Order parameter = (number of bonds in spanning or wrapping cluster) / (total number of bonds)" << std::endl;
    fout << "#C(p,L) = Specific heat = -T dH/dT" << std::endl;
    fout << "#X(p,L) = Susceptibility = dP/dp" << std::endl;
    fout << "#u_i = (number of bonds in the i-th cluster) / (total number of bonds)" << std::endl;
//    std::cout.precision(12);
    fout.precision(numeric_limits<double>::digits10);
    for(size_t i{}; i < limit; ++i){

        fout << (i + 1) / double(limit) << delimiter;
        fout << entropy[i] / double(ensemble_size) << delimiter;
        fout << nob_largest[i] / double(ensemble_size * lp.maxBonds()) << delimiter;
        fout << nob_wraping[i] / double(ensemble_size * lp.maxBonds()) ;
        fout << std::endl;
    }
    fout.close();
//
//
//
//    fout.open(filename_specific_heat_susceptibility);
//    fout << '#' << header_info.str() << std::endl;
//    fout << "#p = occupation probability" << std::endl;
//    fout << "#H(p,L) = Entropy = sum( - u_i * log(u_i))" << std::endl;
//    fout << "#P(p,L) = Order parameter = S_spanning / (2*L*L)" << std::endl;
//    fout << "#C(p,L) = Specific heat = -T dH/dT" << std::endl;
//    fout << "#X(p,L) = Susceptibility = dP/dp, using wrapping cluster" << std::endl;
//    fout << "#u_i = (number of bonds in the i-th cluster) / (total number of bonds)" << std::endl;
//    fout << "#<p>,<C(p,L)>,<X(p,L)>" << std::endl;
////    std::cout.precision(12);
//    fout.precision(numeric_limits<double>::digits10);
//    double p;
//    /*
//     * dp = 1/L^2
//     * dP/dp = d(S_max/(2 L^2))/dp
//     *       = d(S_max)/(2 L^2) * dp
//     *       = d(S_max)/2
//     */
//    for(size_t i{}; i < limit; ++i){
//        p = (i + 1) / double(limit);
//        fout << p << delimiter;
//        fout << -1*(1.0-p)*dHs[i]*length_squared / double(ensemble_size) << delimiter;
//        fout << dPs[i] / double(2*ensemble_size) ;
//        fout << std::endl;
//    }
//    fout.close();
}

void run_v12_regular_bond(int length, int ensemble_size) {
        std::cout << "length=" << length << " and ensemble_size=" << ensemble_size << std::endl;

    long length_squared = length*length;
    long twice_length_squared = 2 * length_squared;

    BondPercolation_pb_v12 bp(length);


    bp.setRandomState(2562449578, true);
    bp.init();

//    lp.viewLattice();
//    lp.viewLatticeByID();
//    lp.viewLatticeExtended();

//    lp.viewCluster();
//    lp.viewRandomizedIndices();
//    lp.viewRemainingIndices();
//    lp.viewBondByID();
    // simulation starts here
    auto a = size_t(ensemble_size);
    std::vector<double> pcs(a), sites_pc(a), bonds_pc(a);
    value_type counter{};
    value_type limit = bp.maxIterationLimit();
    std::vector<double> entropy(limit), dHs(limit);
    std::vector<double> nob_wraping(limit),
            nob_largest(limit), dPs(limit);

    for(int i{} ; i != ensemble_size ; ++i){

        bp.reset();

        bool successful = false;
        auto t_start = std::chrono::system_clock::now();
        counter = 0;
        bool detect_wrapping{true};
        while (true){
//            lp.viewRemainingSites();
//            lp.viewIndices();
            successful = bp.occupy();

            if(successful) {
//                cout << counter << " th site ************########********** last " << lp.lastPlacedSite() << endl;
//                lp.viewRandomizedIndices();
//                lp.viewRemainingIndices();
//                lp.viewSiteByID();
//                lp.viewBondByID();
//                lp.viewCluster();
//                lp.ckeckCluster();
//                auto H1 = lp.entropy_v1();
//                auto H2 = lp.entropy_v2();
//                cout
//                     << H1
//                     << " vs " << H2
//                     << endl;
//                if(abs(H1 - H2) > 1e-10){
//                    cout << "Entropy mismatched" << endl;
//                    exit(0);
//                }
                entropy[counter] += bp.entropy();
                nob_wraping[counter] += bp.numberOfSitesInTheWrappingClusters();
                nob_largest[counter] += bp.numberOfSitesInTheLargestCluster();
//                lp.jump();
//                dHs[counter] += lp.jump_entropy();
//                dPs[counter] += lp.jump_wrapping_cluster();
                if(detect_wrapping && bp.detectWrapping()){

                    pcs[i] = bp.occupationProbability();
                    sites_pc[i] = bp.numberOfSitesInTheWrappingClusters();
                    bonds_pc[i] = bp.numberOfBondsInTheWrappingClusters();
                    cout << "wrapping at "   << pcs[i] << ","
                              << sites_pc[i] << ","
                              << bonds_pc[i] << std::endl;


//                    std::vector<value_type> site, bond;
//                    lattice_percolation.get_cluster_info(site, bond);
//                    for(value_type j{}; j != site.size(); ++j){
//                        fout_s << site[j] << ',';
//                    }
//                    for(value_type j{}; j != bond.size(); ++j){
//                        fout_b << bond[j] <<',';
//                    }
//                    fout_s << std::endl;
//                    fout_b << std::endl;

                    detect_wrapping = false;
//                    break; // only for debugging
                }


                ++counter;
//                if(counter >= 3) break;
            }
            if(counter >= limit){ // twice_length_squared is the number of bonds
                break;
            }
        }
//        fout_jump << lattice_percolation.largestEntropyJump() << "," << lattice_percolation.largestEntropyJump_pc() << std::endl;

        auto t_end = std::chrono::system_clock::now();
        std::cout << "Iteration " << (i+1)
                  //                 << " . Thread " << std::this_thread::get_id()
                  << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << std::endl;

//        cout << "Relabeling time " << lattice_percolation.get_relabeling_time() << endl;
//        lp.viewClusterExtended();
    }
//    cout << "wrapping site id " << lp.wrappingSite_id();
//    cout << "  wrapping site index " << lp.wrappingSite() << endl;
//    lp.viewLattice_by_id();
//    lp.viewLattice_by_gid();
//    lp.viewLattice_sites_by_gid();
//    lp.viewLattice_by_relative_index();
//    lp.viewCluster();

//    fout_b.close();
//    fout_s.close();
//    fout_jump.close();

    std::string tm = getCurrentDateTime();
    const char delimiter = '\t';
    std::ostringstream header_info;
    header_info << "{"
                << R"("length":)" << length
                << R"(,"ensemble_size":)" << ensemble_size
                << R"(,"random_seed":)" << bp.getRandomState()
                << R"(,"signature":")" << bp.getSignature() << "\""
                << R"(,"datetime":")" << tm << "\""
                << R"(,"date":")" << getCurrentDate() << "\""
                << R"(,"time":")" << getCurrentTime() << "\""
                << R"(,"classname":")" << bp.getClassName() << "\""
                //                << R"(,"cols":)" << R"(["S","n_S"])"
                << "}";


    std::string extension = "_L" + std::to_string(length) + '_' + tm + ".txt";
//    std::string filename_s = lattice_percolation.getSignature() + "_cluster_by_site" + extension;
//    std::string filename_b = lattice_percolation.getSignature() + "_cluster_by_bond" + extension;
    std::string filename_critical = bp.getClassName() + "_critical" + extension;
//    std::string filename = lattice_percolation.getSignature() + "_entropy-jump" + extension;
    std::string filename_entropy_order_parameter = bp.getClassName() + "_entropy-order" + extension;

    std::string filename_specific_heat_susceptibility = bp.getClassName() + "_specific_heat-susceptibility" + extension;


//    std::ofstream fout_jump(filename);
//    // JSON formated header
//    fout_jump << '#' << header_info.str() << std::endl;
//    fout_jump << "#each line is an independent realization" << std::endl;
//    fout_jump << "#each line contains information about all clusters at critical point" << std::endl;
//    fout_jump << "#cluster size is measured by number of bonds in it" << std::endl;
//    fout_jump << "#<entropy jump>,<occupation probability at which the jump occurs>" << std::endl;
//
//    std::ofstream fout_s(filename_s);
//    // JSON formated header
//    fout_s << '#' << header_info.str() << std::endl;
//    fout_s << "#each line is an independent realization" << std::endl;
//    fout_s << "#each line contains information about all clusters at critical point" << std::endl;
//    fout_s << "#cluster size is measured by number of site_index_sequence in it" << std::endl;
//
//    std::ofstream fout_b(filename_b);
//    // JSON formated header
//    fout_b << '#' << header_info.str() << std::endl;
//    fout_b << "#each line is an independent realization" << std::endl;
//    fout_b << "#each line contains information about all clusters at critical point" << std::endl;
//    fout_b << "#cluster size is measured by number of bonds in it" << std::endl;

    std::ofstream fout_critical(filename_critical);
    fout_critical << '#' << header_info.str() << std::endl;
    fout_critical << "#data at critical occupation probability or pc" << std::endl;
    fout_critical << "#<pc>,<sites in wrapping cluster>,<bonds in wrapping cluster>" << std::endl;
    for(size_t i{}; i < ensemble_size; ++i){
        if(pcs[i] == 0) continue;
        fout_critical << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
                      << pcs[i] << delimiter
                      << sites_pc[i] << delimiter
                      << bonds_pc[i] << std::endl;
//        cout << pcs[i] << endl;
    }
    fout_critical.close();



    std::ofstream fout(filename_entropy_order_parameter);
    fout << '#' << header_info.str() << std::endl;
    fout << "#<p>,<H(p,L)>,<P1(p,L)>,<P2(p,L)>" << std::endl;
    fout << "#p = occupation probability" << std::endl;
    fout << "#H(p,L) = Entropy = sum( - u_i * log(u_i))" << std::endl;
    fout << "#P1(p,L) = Order parameter = (number of sites in largest cluster) / (total number of sites)" << std::endl;
    fout << "#P2(p,L) = Order parameter = (number of sites in spanning or wrapping cluster) / (total number of sites)" << std::endl;
    fout << "#C(p,L) = Specific heat = -T dH/dT" << std::endl;
    fout << "#X(p,L) = Susceptibility = dP/dp" << std::endl;
    fout << "#u_i = (number of sites in the i-th cluster) / (total number of sites)" << std::endl;
//    std::cout.precision(12);
    fout.precision(numeric_limits<double>::digits10);
    for(size_t i{}; i < limit; ++i){

        fout << (i + 1) / double(limit) << delimiter;
        fout << entropy[i] / double(ensemble_size) << delimiter;
        fout << nob_largest[i] / double(ensemble_size * bp.maxBonds()) << delimiter;
        fout << nob_wraping[i] / double(ensemble_size * bp.maxBonds()) ;
        fout << std::endl;
    }
    fout.close();
//
//
//
//    fout.open(filename_specific_heat_susceptibility);
//    fout << '#' << header_info.str() << std::endl;
//    fout << "#p = occupation probability" << std::endl;
//    fout << "#H(p,L) = Entropy = sum( - u_i * log(u_i))" << std::endl;
//    fout << "#P(p,L) = Order parameter = S_spanning / (2*L*L)" << std::endl;
//    fout << "#C(p,L) = Specific heat = -T dH/dT" << std::endl;
//    fout << "#X(p,L) = Susceptibility = dP/dp, using wrapping cluster" << std::endl;
//    fout << "#u_i = (number of bonds in the i-th cluster) / (total number of bonds)" << std::endl;
//    fout << "#<p>,<C(p,L)>,<X(p,L)>" << std::endl;
////    std::cout.precision(12);
//    fout.precision(numeric_limits<double>::digits10);
//    double p;
//    /*
//     * dp = 1/L^2
//     * dP/dp = d(S_max/(2 L^2))/dp
//     *       = d(S_max)/(2 L^2) * dp
//     *       = d(S_max)/2
//     */
//    for(size_t i{}; i < limit; ++i){
//        p = (i + 1) / double(limit);
//        fout << p << delimiter;
//        fout << -1*(1.0-p)*dHs[i]*length_squared / double(ensemble_size) << delimiter;
//        fout << dPs[i] / double(2*ensemble_size) ;
//        fout << std::endl;
//    }
//    fout.close();
}

void demarcationLine_v12(int length){
    cout << "length " << length << endl;

//    SitePercolation_ps_v12 percolation(length);
//    size_t sz = percolation.maxSites();

    BondPercolation_pb_v12 percolation(length);
    size_t sz = percolation.maxBonds();

    percolation.setRandomState(0, true);
    percolation.init();



//    double t;
    double P{}, P_old{}, dP;
    bool c{false};




    vector<double> dP_list(sz);
    vector<int> color_list(sz);

    size_t i{0};

//    sp.occupy();
//    P_old = sp.numberOfBondsInTheLargestCluster();
//    t = net.occupationProbability();
//    fout << t << '\t' << 0 << '\t' << c << endl;

//    sp.viewLattice_sites_by_gid();
    P_old = 0;
    while (percolation.occupy()){
//        cout << " ****************************** gid "
//             << percolation.lastBondGroupID() << " id " << percolation.lastBond() << endl;
//            << percolation.lastSiteGroupID() << " id " << percolation.lastSite() << endl;

        P = percolation.numberOfBondsInTheLargestCluster();

        dP = abs(P- P_old);

        P_old = P;
        c = percolation.isSelfClusterJump();

        dP_list[i] = dP;
        color_list[i] = c;
        ++i;


//        percolation.viewClusterExtended();
//        percolation.viewLattice_sites_by_gid();
//        percolation.viewLattice();
//        percolation.viewCluster();
//        cout << "****log   " << percolation.occupationProbability() << '\t' << P << '\t' << dP << " self jump "<< c << endl; // takes ~80% of the total time
//        if(i >= 6) break;
//        if(percolation.detectWrapping()){break;}
    }

    auto tm = getCurrentDateTime();
    string signature = percolation.getSignature();
    string filename = signature + "_L" + to_string(length) + "_demarcationLine_" + tm + ".txt";



    ofstream fout(filename);
    stringstream ss;
    ss   << "{\"signature\":" << "\"" << signature << "\""
         << R"*(,"class":")*" << percolation.getClassName() << "\""
         << ",\"L\":" << length  // in case of explosive percolation
         << R"*(,"date":")*" << getCurrentDate() << "\""
         << R"*(,"time":")*" << getCurrentTime() << "\""
         << R"*(,"datetime":")*" << tm << "\""
         << R"*(,"L":")*" << length << "\""
         << R"*(,"length":")*" << length << "\""
         << ",\"cols\":" << "[\"t\", \"dP\", \"color\"]" << "}";


    fout << "#" << ss.str() << endl;
    fout << "#p = occupation probability" << endl;
    fout << "#cluster size = number of bonds in the cluster" << endl;
    fout << "#P = order parameter = largest cluster size / N" << endl;
    fout << "#N = L*L" << endl;
    fout << "#S_max = largest cluster size" << endl;
    fout << "#dS = largest clusteter jump" << endl;
    fout << "#c = 0 or 1 (different cluster jump or self cluster jump)"<< endl;
    fout << "#<t>\t<dP>\t<c>" << endl;


    for(size_t i{}; i < sz; ++i){
        if (dP_list[i] != 0) {
            fout << (i+1)/double(sz) << '\t' << dP_list[i] << '\t' << color_list[i] << endl;

//            auto ch = (color_list[i] == 0)? "b" : "r";
//            fout << (i+1)/double(sz) << '\t' << dP_list[i] << '\t'
//                 << ch << endl;
            if(color_list[i] == 0){
                cout << "got one" << endl;
            }
        }
    }

    fout.close();

    auto sss = set<int>(color_list.begin(), color_list.end());
    cout << sss.size() << endl;
}