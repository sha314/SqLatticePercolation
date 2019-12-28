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

using namespace std;

void test_v12(int argc, char **argv) {
    cout << "length        = argv[1]" << endl;
    cout << "ensemble size = argv[2]" << endl;
//    cout << "rsbd l        = argv[3]" << endl;

    int length = stoi(argv[1]);
    int ensemble_size = stoi(argv[2]);


//    test_v12_lattice(length);
//    test_v12_percolation(length);

    run_v12_regular(length, ensemble_size);
}

void test_v12_lattice(int length) {
    Lattice_v12 lattice(length);

//    lattice.view_as_assigned();
//    lattice.view_all();
    lattice.view_from_str();
//    lattice.view_all_v2();
//    lattice.view_sites();
//    lattice.view_sites_by_id();
//    lattice.view_sites_by_gid();
//    lattice.view_by_id();
//    lattice.view_by_gid();
//    lattice.view_sites_by_relative_index();
//    lattice.view_by_relative_index();

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


void test_v12_percolation(int length) {
    SitePercolation_ps_v12 percolation(length);
    percolation.setRandomState(0, false);
    percolation.init();

//    percolation.viewCluster();
//    percolation.viewLattice_by_id();
//    percolation.viewLattice_by_gid();

//    percolation.occupy();
    int i=0;
    while(percolation.occupy()) {
        cout << "********************** last site " << percolation.lastSite() << endl;
//    percolation.viewCluster();
//    percolation.viewLattice_by_id();
//        percolation.viewLattice_by_gid();
//    percolation.viewLattice();
        percolation.viewLattice_by_relative_index();
        if(percolation.detectWrapping()){
//            percolation.viewLattice_by_relative_index();
            cout << "wrapping site " << percolation.wrappingSite() << endl;
            cout << "wrapping site id " << percolation.wrappingSite_id() << endl;
            cout << "p = " << percolation.occupationProbability() << endl;
            break;
        }

        ++i;
//        if(i >= 11) break;
    }


}


void run_v12_regular(int length, int ensemble_size) {

    std::cout << "length=" << length << " and ensemble_size=" << ensemble_size << std::endl;

    long length_squared = length*length;
    long twice_length_squared = 2 * length_squared;

    SitePercolation_ps_v12 lp(length);
//    SitePercolationBallisticDeposition_L2_v3 lp(length, true);
//    SitePercolationRSBD_L1_v10 lp(length, true);
//    SitePercolationRSBD_L2_v10 lp(length, true);

    lp.setRandomState(0, true);
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
//                entropy[counter] += lp.entropy();
//                nob_wraping[counter] += lp.numberOfBondsInTheWrappingClusters();
//                nob_largest[counter] += lp.numberOfBondsInTheLargestCluster();
//                lp.jump();
//                dHs[counter] += lp.jump_entropy();
//                dPs[counter] += lp.jump_wrapping_cluster();
                if(detect_wrapping && lp.detectWrapping()){
                    pcs[i] = lp.occupationProbability();
                    sites_pc[i] = lp.numberOfSitesInTheWrappingClusters();
                    bonds_pc[i] = lp.numberOfBondsInTheWrappingClusters();


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

//    fout_b.close();
//    fout_s.close();
//    fout_jump.close();

    std::string tm = getCurrentTime();
    const char delimiter = '\t';
    std::ostringstream header_info;
    header_info << "{"
                << R"("length":)" << length
                << R"(,"ensemble_size":)" << ensemble_size
                << R"(,"random_seed":)" << lp.getRandomState()
                << R"(,"signature":")" << lp.getSignature() << "\""
                << R"(,"datetime":")" << tm << "\""
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



//    std::ofstream fout(filename_entropy_order_parameter);
//    fout << '#' << header_info.str() << std::endl;
//    fout << "#<p>,<H(p,L)>,<P1(p,L)>,<P2(p,L)>" << std::endl;
//    fout << "#p = occupation probability" << std::endl;
//    fout << "#H(p,L) = Entropy = sum( - u_i * log(u_i))" << std::endl;
//    fout << "#P1(p,L) = Order parameter = (number of bonds in largest cluster) / (total number of bonds)" << std::endl;
//    fout << "#P2(p,L) = Order parameter = (number of bonds in spanning or wrapping cluster) / (total number of bonds)" << std::endl;
//    fout << "#C(p,L) = Specific heat = -T dH/dT" << std::endl;
//    fout << "#X(p,L) = Susceptibility = dP/dp" << std::endl;
//    fout << "#u_i = (number of bonds in the i-th cluster) / (total number of bonds)" << std::endl;
////    std::cout.precision(12);
//    fout.precision(numeric_limits<double>::digits10);
//    for(size_t i{}; i < limit; ++i){
//
//        fout << (i + 1) / double(limit) << delimiter;
//        fout << entropy[i] / double(ensemble_size) << delimiter;
//        fout << nob_largest[i] / double(ensemble_size * lp.maxBonds()) << delimiter;
//        fout << nob_wraping[i] / double(ensemble_size * lp.maxBonds()) ;
//        fout << std::endl;
//    }
//    fout.close();
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
