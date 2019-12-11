//
// Created by shahnoor on 9/20/18.
//

#include "test_percolation.h"
#include "../types.h"
#include "../percolation/percolation.h"
#include "../util/time_tracking.h"
#include "../percolation/site/percolation_site_v9.h"
#include "../percolation/bond/percolation_bond_v1.h"
#include "../percolation/site/percolation_site_v10.h"
#include "../percolation/bond/percolation_bond_v2.h"
#include "../percolation/site/percolation_site_explosive.h"
#include "../percolation/site/percolation_site_ballistic.h"
#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>

using namespace std;

/**
 * writes x and y coordinate for visual data of the lattice
 * @param argc
 * @param argv
 */
void visual(value_type length){

    cout << "length " << length << endl;

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    SitePercolation_ps_v9 lattice_percolation(length, true);

    string tm = getCurrentTime();

    // simulation starts here
    bool successful{false};
    size_t counter{};
    while (true){
        successful = lattice_percolation.occupy();
        if(successful) {
            if(lattice_percolation.detectWrapping()){
                break;
            }
            ++counter;
        }
        if(counter >= lattice_percolation.maxIterationLimit()){ // twice_length_squared is the number of bonds
            break;
        }
    }

    lattice_percolation.writeVisualLatticeData(
            lattice_percolation.getSignature()+"L"+to_string(length)+"-visual-"+tm+".csv",
            false);

}
/**
 * This function simulates for all fundamental data for site percolation on square lattice.
 * 1. cluster size distribution by sites and bonds
 * 2. critical point data
 * 3. entropy jump
 * 4. entropy and order parameter
 *
 * length of the lattice and the ensemble size must be provided via command line arguments.
 * @param argc
 * @param argv
 */
void simulate_site_percolation(value_type length, value_type ensemble_size) {

    std::cout << "length " << length << " ensemble_size " << ensemble_size << std::endl;

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    SitePercolation_ps_v10 lattice_percolation(length, true);
    lattice_percolation.setRandomState(0, true);
    lattice_percolation.init();

    std::string tm = getCurrentTime();
    const char delimiter = '\t';
    std::ostringstream header_info;
    header_info << "{"
                << R"("length":)" << length
                << R"(,"ensemble_size":)" << ensemble_size
                << R"(,"random_seed":)" << lattice_percolation.getRandomState()
                << R"(,"signature":")" << lattice_percolation.getSignature() << "\""
                << R"(,"time":")" << tm << "\""
                //                << R"(,"cols":)" << R"(["S","n_S"])"
                << "}";


    std::string extension = "_L" + std::to_string(length) + '_' + tm + ".txt";
//    std::string filename_s = lattice_percolation.getSignature() + "_cluster_by_site" + extension;
//    std::string filename_b = lattice_percolation.getSignature() + "_cluster_by_bond" + extension;
    std::string filename_critical = lattice_percolation.getSignature() + "_critical" + extension;
//    std::string filename = lattice_percolation.getSignature() + "_entropy-jump" + extension;
    std::string filename_entropy_order_parameter = lattice_percolation.getSignature() + extension;


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


    // simulation starts here
    std::vector<double> pcs(ensemble_size), sites_pc(ensemble_size), bonds_pc(ensemble_size);
    value_type counter{};
    std::vector<double> entropy(lattice_percolation.maxIterationLimit());
    std::vector<double> nob_wraping(lattice_percolation.maxIterationLimit()),
            nob_largest(lattice_percolation.maxIterationLimit());

    for(value_type i{} ; i != ensemble_size ; ++i){

        lattice_percolation.reset();

        bool successful = false;
        auto t_start = std::chrono::system_clock::now();
        counter = 0;
        bool detect_wrapping{true};
        while (true){
            successful = lattice_percolation.occupy();
            if(successful) {
                entropy[counter] += lattice_percolation.entropy();
                nob_wraping[counter] += lattice_percolation.numberOfBondsInTheWrappingClusters();
                nob_largest[counter] += lattice_percolation.numberOfBondsInTheLargestCluster();
                lattice_percolation.jump();
                if(detect_wrapping && lattice_percolation.detectWrapping()){
                    pcs[i] = lattice_percolation.occupationProbability();
                    sites_pc[i] = lattice_percolation.numberOfSitesInTheWrappingClusters();
                    bonds_pc[i] = lattice_percolation.numberOfBondsInTheWrappingClusters();


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
            }
            if(counter >= lattice_percolation.maxIterationLimit()){ // twice_length_squared is the number of bonds
                break;
            }
        }
//        fout_jump << lattice_percolation.largestEntropyJump() << "," << lattice_percolation.largestEntropyJump_pc() << std::endl;
        {
            auto t_end = std::chrono::system_clock::now();
            std::cout << "Iteration " << i
                      //                 << " . Thread " << std::this_thread::get_id()
                      << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << std::endl;
        }
//        cout << "Relabeling time " << lattice_percolation.get_relabeling_time() << endl;
    }

//    fout_b.close();
//    fout_s.close();
//    fout_jump.close();

    std::ofstream fout_critical(filename_critical);
    fout_critical << '#' << header_info.str() << std::endl;
    fout_critical << "#data at critical occupation probability or pc" << std::endl;
    fout_critical << "#<pc>,<sites in wrapping cluster>,<bonds in wrapping cluster>" << std::endl;
    for(size_t i{}; i < ensemble_size; ++i){
        if(pcs[i] == 0) continue;
        fout_critical << std::setprecision(numeric_limits<double>::digits10 + 1) <<  pcs[i] << delimiter
                      << sites_pc[i] << delimiter
                      << bonds_pc[i] << std::endl;
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
    fout.precision(numeric_limits<double>::digits10 + 1);
    for(size_t i{}; i < lattice_percolation.maxIterationLimit(); ++i){

        fout << (i + 1) / double(lattice_percolation.maxIterationLimit()) << delimiter;
        fout << entropy[i] / double(ensemble_size) << delimiter;
        fout << nob_largest[i] / double(ensemble_size * lattice_percolation.maxBonds()) << delimiter;
        fout << nob_wraping[i] / double(ensemble_size * lattice_percolation.maxBonds()) ;
        fout << std::endl;
    }
    fout.close();
}

/**
 * simulating site percolation at each step
 * @param length
 * @param ensemble_size
 */
void simulate_site_percolation_detailed(value_type length, value_type ensemble_size) {
    cout << "length " << length << " ensemble_size " << ensemble_size << endl;
    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    SitePercolation_ps_v9 lattice_percolation(length, true);

    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length
                << ",\"ensemble_size\":" << ensemble_size
                << ",\"signature\":\"" << lattice_percolation.getSignature() << "\""
                << "}" ;

    string tm = getCurrentTime();

    string filename_entropy = lattice_percolation.getSignature()  + "entropy_L" + to_string(length) + '_' + tm;
    string filename_order_parameter = lattice_percolation.getSignature() + "order-parameter_L" + to_string(length) + '_' + tm;


    filename_entropy += ".csv";
    filename_order_parameter += ".csv";



    ofstream fout_entropy(filename_entropy);
    fout_entropy << '#' << header_info.str() << endl;
    fout_entropy << "#each line contains entropy values for one independent realization" << endl;

    ofstream fout_P(filename_order_parameter);
    fout_P << '#' << header_info.str() << endl;
    fout_P << "#each line contains order parameter values for one independent realization" << endl;
    fout_P << "#order parameter is measured in terms of wrapping cluster" << endl;

    // simulation starts here
    value_type counter{};

    double H,P;
    for(value_type i{} ; i != ensemble_size ; ++i){

        lattice_percolation.reset();

        bool successful = false;
        auto t_start = std::chrono::system_clock::now();
        counter = 0;
        bool wrapping_written{false};
        while (true){
            successful = lattice_percolation.occupy();
            if(successful) {
                H = lattice_percolation.entropy();
                lattice_percolation.detectWrapping();
                P = lattice_percolation.numberOfBondsInTheWrappingClusters();
                fout_entropy << H << ',';
                fout_P << P << ',';
                ++counter;
            }
            if(counter >= lattice_percolation.maxIterationLimit()){ // twice_length_squared is the number of bonds
                break;
            }
        }
        fout_entropy << endl;
        fout_P << endl;
        {
            auto t_end = std::chrono::system_clock::now();
            cout << "Iteration " << i
                 //                 << " . Thread " << std::this_thread::get_id()
                 << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << endl;
        }
//        cout << "Relabeling time " << lattice_percolation.get_relabeling_time() << endl;
    }

    fout_entropy.close();
    fout_P.close();
}



void simulate_bond_percolation(int argc, char **argv) {
    value_type length = atoi(argv[1]);
    value_type ensemble_size = atoi(argv[2]);

//    value_type length = 5;
    BondPercolation_pb_v1 bp(length);
    value_type total_bonds= 2*length*length;
    double length_squared = length*length;

    value_type j = 0;
    vector<double> entrpy(total_bonds), nos1(total_bonds), nos2(total_bonds);


    bool successful{false};

    for(value_type e{}; e < ensemble_size; ++e) {
        bp.reset();
//        bp.viewBondByID();
//        bp.viewClusterExtended();

        cout << "Iteration " << e << " : ";
        auto t0 = chrono::system_clock::now();
        j = 0;
        while (true) {
            successful = bp.occupy();
            if (successful) {
//                cout << "*******************************************" << endl;
//                entrpy[j] += bp.entropy_slow();
                entrpy[j] += bp.entropy();
                auto index = bp.lastPlacedBond();

//                cout << j << " th bond" << index << endl;
//                bp.viewBondByID();
//                bp.viewClusterExtended();
//                cout << "Entropy " << H_tmp << endl;


                ++j;
            }
            if (j >= total_bonds) { // length_squared is the number of site
                break;
            }
        }
        auto t1 = chrono::system_clock::now();
        cout << "elapsed time " << chrono::duration<double>(t1 - t0).count() << " sec" << endl;
//        bp.viewBondByID();
//        bp.viewClusterExtended();
    }


    for(value_type e{}; e < entrpy.size(); ++e) {
        entrpy[e] /= double(ensemble_size);
    }


    stringstream ss;
    ss << "{\"length\":" << length
       <<",\"ensemble_size\":"<< ensemble_size
       << ",\"signature\":\"" << bp.getSignature() << "\"}" << endl;
    string header  = ss.str();
    header += R"***(#<p>	<H(p,L)>	<P1(p,L)>	<P2(p,L)>
#p = occupation probability
#H(p,L) = Entropy = sum( - u_i * log(u_i))
#P1(p,L) = Order parameter = (number of site_index_sequence in largest cluster) / (total number of site_index_sequence)
#P2(p,L) = Order parameter = (number of site_index_sequence in spanning or wrapping cluster) / (total number of site_index_sequence)
#C(p,L) = Specific heat = -T dH/dT
#X(p,L) = Susceptibility = dP/dp
#u_i = (number of site_index_sequence in the i-th cluster) / (total number of site_index_sequence))***";

    string tm = currentTime();
    string filename = bp.getSignature() + "_" + to_string(length) + "_" + tm + ".txt";


    ofstream fout(filename);
    fout << header << endl;
    double p, H, P1, P2;
    for(size_t i{}; i!= total_bonds; ++i){
        p = (i+1)/ double(total_bonds);
        H = entrpy[i];                      // Entropy
//        P1 = nos1[i] / length_squared;    // Order parameter
//        P2 = nos2[i] / length_squared;    // Order parameter

        fout << p << '\t' << H << '\t' << '-' << '\t' << '-' << endl;
    }
    fout.close();
}

/**
 * This function simulates for all fundamental data for site percolation on square lattice.
 * 1. cluster size distribution by sites and bonds
 * 2. critical point data
 * 3. entropy jump
 * 4. entropy and order parameter
 *
 * length of the lattice and the ensemble size must be provided via command line arguments.
 * @param argc
 * @param argv
 */
void simulate_bond_percolation_v2(value_type length, value_type ensemble_size) {

    cout << "length " << length << " ensemble_size " << ensemble_size << endl;

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    BondPercolation_pb_v1 lattice_percolation(length, true);

    ostringstream header_info;
    header_info << "{"
                << R"("length":)" << length
                << R"(,"ensemble_size":)" << ensemble_size
                << R"(,"signature":")" << lattice_percolation.getSignature() << "\""
                << R"(,"delimiter":)" << "\"\t\""
                << "}" ;

    std::string tm = getCurrentTime();
    std::string extension  = "_L" + std::to_string(length) + '_' + tm + ".txt";
//    std::string filename_s = lattice_percolation.getSignature() + "_cluster_by_site" + extension;
//    std::string filename_b = lattice_percolation.getSignature() + "_cluster_by_bond" + extension;
    std::string filename_critical = lattice_percolation.getSignature() + "_critical" + extension;
//    std::string filename = lattice_percolation.getSignature() + "_entropy-jump" + extension;
    std::string filename_entropy_order_parameter = lattice_percolation.getSignature() + extension;

//    ofstream fout_jump(filename);
//    // JSON formated header
//    fout_jump << '#' << header_info.str() << endl;
//    fout_jump << "#each line is an independent realization" << endl;
//    fout_jump << "#each line contains information about all clusters at critical point" << endl;
//    fout_jump << "#cluster size is measured by number of bonds in it" << endl;
//    fout_jump << "#<entropy jump>,<occupation probability at which the jump occurs>" << endl;

//    ofstream fout_s(filename_s);
//    // JSON formated header
//    fout_s << '#' << header_info.str() << endl;
//    fout_s << "#each line is an independent realization" << endl;
//    fout_s << "#each line contains information about all clusters at critical point" << endl;
//    fout_s << "#cluster size is measured by number of site_index_sequence in it" << endl;
//
//    ofstream fout_b(filename_b);
//    // JSON formated header
//    fout_b << '#' << header_info.str() << endl;
//    fout_b << "#each line is an independent realization" << endl;
//    fout_b << "#each line contains information about all clusters at critical point" << endl;
//    fout_b << "#cluster size is measured by number of bonds in it" << endl;



    // simulation starts here
    value_type counter{};
    vector<double> critical_pc(ensemble_size), critical_sites(ensemble_size), critical_bonds(ensemble_size);
    vector<double> entropy(lattice_percolation.maxIterationLimit());
    vector<double> nob_wraping(lattice_percolation.maxIterationLimit()),
            nob_largest(lattice_percolation.maxIterationLimit());

    for(value_type i{} ; i != ensemble_size ; ++i){

        lattice_percolation.reset();

        bool successful = false;
        auto t_start = std::chrono::system_clock::now();
        counter = 0;
        bool wrapping_written{false};
        while (true){
            successful = lattice_percolation.occupy();
            if(successful) {
                entropy[counter] += lattice_percolation.entropy();
                nob_wraping[counter] += lattice_percolation.numberOfBondsInTheWrappingClusters();
                nob_largest[counter] += lattice_percolation.numberOfBondsInTheLargestCluster();
                lattice_percolation.jump();
                if(!wrapping_written && lattice_percolation.detectWrapping()){
                    critical_pc[i] = lattice_percolation.occupationProbability();
                    critical_sites[i] = lattice_percolation.numberOfSitesInTheWrappingClusters();
                    critical_bonds[i] = lattice_percolation.numberOfBondsInTheWrappingClusters();


//                    vector<value_type> site, bond;
//
//                    lattice_percolation.get_cluster_info(site, bond);
//
//                    for(value_type j{}; j != site.size(); ++j){
//                        fout_s << site[j] << ',';
//                    }
//                    for(value_type j{}; j != bond.size(); ++j){
//                        fout_b << bond[j] <<',';
//                    }


//                    fout_s << endl;
//                    fout_b << endl;
                    wrapping_written = true;
                }


                ++counter;
            }
            if(counter >= lattice_percolation.maxIterationLimit()){ // twice_length_squared is the number of bonds
                break;
            }
        }
//        fout_jump << lattice_percolation.largestEntropyJump() << "," << lattice_percolation.largestEntropyJump_pc() << endl;

        auto t_end = std::chrono::system_clock::now();
        cout << "Iteration " << i
             //                 << " . Thread " << std::this_thread::get_id()
             << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << endl;

//        cout << "Relabeling time " << lattice_percolation.get_relabeling_time() << endl;
    }

//    fout_b.close();
//    fout_s.close();

//    fout_jump.close();

    ofstream fout_critical(filename_critical);
    fout_critical << '#' << header_info.str() << endl;
    fout_critical << "#data at critical occupation probability or pc" << endl;
    fout_critical << "#<pc><site_index_sequence in wrapping cluster><bonds in wrapping cluster>" << endl;

    ofstream fout(filename_entropy_order_parameter);
    fout << '#' << header_info.str() << endl;
    fout << "#<p>,<H(p,L)>,<P1(p,L)>,<P2(p,L)>" << endl;
    fout << "#p = occupation probability" << endl;
    fout << "#H(p,L) = Entropy = sum( - u_i * log(u_i))" << endl;
    fout << "#P1(p,L) = Order parameter = (number of bonds in largest cluster) / (total number of bonds)" << endl;
    fout << "#P2(p,L) = Order parameter = (number of bonds in spanning or wrapping cluster) / (total number of bonds)" << endl;
    fout << "#C(p,L) = Specific heat = -T dH/dT" << endl;
    fout << "#X(p,L) = Susceptibility = dP/dp" << endl;
    fout << "#u_i = (number of bonds in the i-th cluster) / (total number of bonds)" << endl;
    for(size_t i{}; i < lattice_percolation.maxIterationLimit(); ++i){
        fout << (i+1) / double(lattice_percolation.maxIterationLimit()) << "\t";
        fout << entropy[i] / double(ensemble_size) << "\t";
        fout << nob_largest[i] / double(ensemble_size * lattice_percolation.maxBonds()) << "\t";
        fout << nob_wraping[i] / double(ensemble_size * lattice_percolation.maxBonds()) ;
        fout << endl;
    }
    fout.close();

    for(size_t i{}; i < ensemble_size; ++i){
        fout_critical   << critical_pc[i] << '\t'
                        << critical_sites[i] << '\t'
                        << critical_bonds[i] << endl;
    }
    fout_critical.close();
}

/**
 * TODO : SitePercolation_ps_v10 works but SitePercolation_ps_v9 fails because of allocation error
 * @param length
 * @param ensemble_size
 */
void cluster_size_distribution_site(value_type length, value_type ensemble_size) {
    cout << "length " << length << " ensemble_size " << ensemble_size << endl;

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    vector<double> cluster_size_distro;

    SitePercolation_ps_v10 lp(length, true);

    for(value_type i{} ; i != ensemble_size ; ++i){
//        cout << "line " << __LINE__ << endl;

        lp.reset();

        bool successful = false;
        auto t_start = std::chrono::system_clock::now();
        size_t counter = 0;
        bool wrapping_written{false};
//        cout << "line " << __LINE__ << endl;

        while (true){
            successful = lp.occupy();
            if(successful) {
                if(lp.detectWrapping()) {
                    vector<double> tmp = lp.clusterSizeDistribution();
//                    cout << "after returned " << tmp.size() << endl;
                    if (tmp.size() > cluster_size_distro.size()) {
                        cluster_size_distro.resize(tmp.size());
                    }
                    for (size_t k{}; k < tmp.size(); ++k) {
                        cluster_size_distro[k] += tmp[k];
                    }
                    break;
                }
                ++counter;
            }
            if(counter >= lp.maxIterationLimit()){ // twice_length_squared is the number of bonds
                break;
            }
        }

        auto t_end = std::chrono::system_clock::now();
        cout << "Iteration " << i+1
             //                 << " . Thread " << std::this_thread::get_id()
             << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << endl;

//        cout << "Relabeling time " << lp.get_relabeling_time() << endl;
    }


    cout  << "writing to file" << endl;

    ostringstream header_info;
    header_info << "{"
                << R"("length":)" << length
                << R"(,"ensemble_size":)" << ensemble_size
                << R"(,"signature":")" << lp.getSignature() << "\""
                << R"(,"cols":)" << R"(["S","n_S"])"
                << "}";


    std::string tm = getCurrentTime();
    std::string extension  = "_L" + std::to_string(length) + '_' + tm + ".txt";

    std::string filename = lp.getSignature() + "_cluster-size-distribution" + extension;


    ofstream fout(filename);
    fout << '#' << header_info.str() << endl;
    fout << "#<S>,<n_S>" << endl;
    fout << "#S=cluster size" << endl;
    fout << "#n_S=number of cluster of size S" << endl;
    fout << "#cluster size is measured by number of bonds on a cluster in site percolation" << endl;
//    cout << cluster_size_distro.size() << endl;
    double area=0;
    for(size_t S{0}; S < cluster_size_distro.size(); ++S){
        if (cluster_size_distro[S] == 0) continue; // no need to write zero values
        fout << S << "\t" << cluster_size_distro[S]/(ensemble_size *2* length_squared);
        fout << endl;
        area += S * cluster_size_distro[S]/(ensemble_size * 2 * length_squared);
    }
    fout.close();
    cout << "area under curve " << area << endl;
}


/**
 * TODO : SitePercolation_ps_v10 works but SitePercolation_ps_v9 fails because of allocation error
 * @param length
 * @param ensemble_size
 */
void cluster_size_distribution_bond(value_type length, value_type ensemble_size) {
    cout << "length " << length << " ensemble_size " << ensemble_size << endl;

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    vector<double> cluster_size_distro;
    BondPercolation_pb_v2 lp(length, true);
    lp.setRandomState(0, true);
    lp.init();


    for(value_type i{} ; i != ensemble_size ; ++i){
//        cout << "line " << __LINE__ << endl;

        lp.reset();

        bool successful = false;
        auto t_start = std::chrono::system_clock::now();
        size_t counter = 0;
        bool wrapping_written{false};
//        cout << "line " << __LINE__ << endl;

        while (true){
            successful = lp.occupy();
            if(successful) {
                if(lp.detectWrapping()) {
                    vector<double> tmp = lp.clusterSizeDistribution();
//                    cout << "after returned " << tmp.size() << endl;
                    if (tmp.size() > cluster_size_distro.size()) {
                        cluster_size_distro.resize(tmp.size());
                    }
                    for (size_t k{}; k < tmp.size(); ++k) {
                        cluster_size_distro[k] += tmp[k];
                    }
                    break;
                }
                ++counter;
            }
            if(counter >= lp.maxIterationLimit()){ // twice_length_squared is the number of bonds
                break;
            }
        }

        auto t_end = std::chrono::system_clock::now();
        cout << "Iteration " << i+1
             //                 << " . Thread " << std::this_thread::get_id()
             << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << endl;

//        cout << "Relabeling time " << lp.get_relabeling_time() << endl;
    }


    cout  << "writing to file" << endl;

    std::string tm = getCurrentTime();
    ostringstream header_info;
    header_info << "{"
                << R"("length":)" << length
                << R"(,"ensemble_size":)" << ensemble_size
                << R"(,"signature":")" << lp.getSignature() << "\""
                << R"(,"cols":)" << R"(["S","n_S"])"
                << R"(,"random_state":)" << lp.getRandomState()
                << R"(,"date":")" << tm << "\""
                << "}";



    std::string extension  = "_L" + std::to_string(length) + '_' + tm + ".txt";

    std::string filename = lp.getSignature() + "_cluster-size-distribution" + extension;


    ofstream fout(filename);
    fout << '#' << header_info.str() << endl;
    fout << "#<S>,<n_S>" << endl;
    fout << "#S=cluster size" << endl;
    fout << "#n_S=number of cluster of size S" << endl;
    fout << "#cluster size is measured by number of site_index_sequence on a cluster in bond percolation" << endl;
//    cout << cluster_size_distro.size() << endl;
    double area=0;
    for(size_t S{0}; S < cluster_size_distro.size(); ++S){
        if (cluster_size_distro[S] == 0) continue; // no need to write zero values
        fout << S << "\t" << cluster_size_distro[S]/(ensemble_size * length_squared);
        fout << endl;
        area += S * cluster_size_distro[S]/(ensemble_size * length_squared);
    }
    fout.close();
    cout << "area under curve " << area << endl;
}

void entropyJumps(int argc, char** argv){
    value_type length = atoi(argv[1]);
    value_type ensemble_size = atoi(argv[2]);


    cout << "length " << length << " ensemble_size " << ensemble_size << endl;

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    BondPercolation_pb_v1 bp(length, true);

    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length
                << ", \"ensemble_size\":" << ensemble_size
                << ", \"signature\":\"" << bp.getSignature() << "\""
                << "}" ;

    string tm = currentTime();

    string filename = bp.getSignature() + "_entropy-jump_" + to_string(length) + '_' + tm;
    filename += ".csv";

    ofstream fout(filename);
    // JSON formated header
    fout << '#' << header_info.str() << endl;
    fout << "#each line is an independent realization" << endl;
    fout << "#each line contains information about all clusters at critical point" << endl;
    fout << "#cluster size is measured by number of site_index_sequence in it" << endl;

    value_type counter{};
    double H{};
    for(value_type i{} ; i != ensemble_size ; ++i){

        bp.reset();

        bool successful = false;
        auto t_start = std::chrono::system_clock::now();
        counter = 0;
        while (true){
            successful = bp.occupy();
            if(successful) {
                H = bp.entropy();
//                cout << H << endl;
                bp.jump();

                ++counter;
            }
            if(counter >= bp.maxIterationLimit()){ // twice_length_squared is the number of bonds
                break;
            }
        }
        fout << bp.largestEntropyJump() << "," << bp.largestEntropyJump_pc() << endl;
        {
            auto t_end = std::chrono::system_clock::now();
            cout << "Iteration " << i
                 //                 << " . Thread " << std::this_thread::get_id()
                 << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << endl;
        }
//        cout << "Relabeling time " << bp.get_relabeling_time() << endl;
    }

    fout.close();
}



/**
 * data for log(Mass) vs log(Length) curve
 * for spanning sites only
 */
void percolation_to_get_fractal(
        value_type ensemble_size, value_type L_start, value_type l_end, value_type delta_L
)
{
    double M_spanning_sites{}, M_sites{}, M_spanning_bonds{}, M_bonds{}, N_cluster{};
    ofstream fout("spanning_data_diff_length.txt");
    fout << "#All data when spanning occurs for different length" << endl;
    fout << "#<Length>\t<Number of cluster>\t<Number of spanning site_index_sequence>\t<Total number of site_index_sequence> \t<Number of spanning bonds>" << endl;
    clock_t t, t_outer;
    for(value_type len{L_start}; len <= l_end ; len += delta_L){
        cout << "Length " << len << " : " << endl;
        SitePercolation_ps_v9 sp(len);

        M_spanning_sites=0;
        M_sites=0;
        M_spanning_bonds=0;
        M_bonds=0;
        N_cluster=0;

        t_outer = clock();
        for(value_type i{} ; i != ensemble_size ; ++i){
            t = clock();
            sp.reset();
            while(!sp.detectWrapping()){
                sp.occupy();
            }
            M_spanning_sites += sp.numberOfSitesInTheSpanningClusters();
            M_sites += sp.numberOfOccupiedSite();
            M_spanning_bonds += sp.numberOfBondsInTheSpanningClusters();
            N_cluster += sp.numberOfcluster();
            cout << "\t\tIteration " << i << " . birthTime " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
        }
        fout << len << '\t'
             << N_cluster / double(ensemble_size)       << '\t'
             << M_spanning_sites /double(ensemble_size) << '\t'
             << M_sites /double(ensemble_size)          << '\t'
             << M_spanning_bonds /double(ensemble_size) << '\t'
             << M_bonds /double(ensemble_size)          << endl; // writing data to file
        cout << "Length " << len << " . Time taken "
             << getFormattedTime((clock() - t_outer) / double(CLOCKS_PER_SEC)) << endl;
    }
    fout.close();
}


/**
 * data for log(Mass) vs log(Length) curve
 * for spanning sites only
 */
void percolation_fractalDimension_by_spanning_site(
        value_type ensemble_size, value_type L_start, value_type l_end, value_type delta_L
)
{
    double M{};
    ofstream fout("mass_length_data.txt");
    fout << "#Measuring fractal dimension by only spanning site_index_sequence" << endl;
    fout << "#<Length>\t<Mass>" << endl;
    clock_t t, t_outer;
    for(value_type len{L_start}; len <= l_end ; len += delta_L){
        cout << "Length " << len << " : " << endl;
        SitePercolation_ps_v9 sp(len);
        M = 0;
        t_outer = clock();
        for(value_type i{} ; i != ensemble_size ; ++i){
            t = clock();
            sp.reset();
            while(!sp.detectWrapping()){
                sp.occupy();
            }
            M += sp.numberOfSitesInTheSpanningClusters();
            cout << "\t\tIteration " << i << " . birthTime " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
        }
        fout << len << '\t' << M / double(ensemble_size) << endl; // writing data to file
        cout << "Length " << len << " . Time taken "
             << getFormattedTime((clock() - t_outer) / double(CLOCKS_PER_SEC)) << endl;
    }
    fout.close();
}



void measure_entropy_by_site(int argc, char** argv){
    value_type length = atoi(argv[1]);
    value_type ensemble_size = atoi(argv[2]);


    cout << "length " << length << " ensemble_size " << ensemble_size << endl;

    value_type length_squared = length*length;
    vector<double> entropy(length_squared);
//    vector<double> entropy2(length_squared);

    SitePercolation_ps_v9 sp(length, true);
    value_type j{};
    for(value_type i{} ; i != ensemble_size ; ++i){

        sp.reset();

        bool successful = false;
        auto t_start = std::chrono::system_clock::now();
        j = 0;
        while (true){
            successful = sp.occupy();
            if(successful) {
                entropy[j] += sp.entropy(); // slower method
//                entropy[j] += sp.entropy_v5_site_threaded(); // slower method
//                entropy2[j] += sp.entropy_v4(2); // faster method

//                cout << '\t' << entropy[j]  << '\t' << entropy2[j] << endl;

                ++j;
            }
            if(j >= length_squared){ // length_squared is the number of site
                break;
            }
        }
        {
            auto t_end = std::chrono::system_clock::now();
            cout << "Iteration " << i
                 //                 << " . Thread " << std::this_thread::get_id()
                 << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << endl;
        }
//        cout << "Relabeling time " << sp.get_relabeling_time() << endl;
    }

    // Taking Average
    double ensmbl = double(ensemble_size);
    for(size_t i{}; i!= length_squared ; ++i){
        entropy[i] /= ensmbl;
//        entropy2[i] /= ensmbl;
    }


    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length
                << ", \"ensemble_size\":" << ensemble_size
                << ", \"signature\":\"" << sp.getSignature() << "\""
                << "}" << endl;

    string tm = currentTime();

    string filename = sp.getSignature() + "_entropy_by_site_" + to_string(length) + '_' + tm;
    filename += ".txt";

    ofstream fout(filename);
    // JSON formated header
    fout << header_info.str();
    fout << "#p=Occupation probability" << endl;
    fout << "#H=entropy" << endl;
    fout << "#cluster length is measured by site" << endl;
    fout << "#<p>\t<H>\t<H fast>\n";

    for(value_type i{}; i != length_squared;++i){
        fout << (i+1)/double(length_squared)
             << '\t' << entropy[i]
             //             << '\t' << entropy2[i]
             << endl;
    }
    fout.close();
}

void test_bond_percolation(value_type length, value_type ensemble_size) {
    cout << "length " << length << " ensemble_size " << ensemble_size << endl;

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    vector<double> cluster_size_distro;
    double area=0;
    BondPercolation_pb_v2 lp(length, true);
    lp.setRandomState(0, true);
    lp.init();
//    lp.reset();

    bool successful = false;
    auto t_start = std::chrono::system_clock::now();
    size_t counter = 0;
    bool wrapping_written{false};
//        cout << "line " << __LINE__ << endl;
//    lp.ckeckCluster();
//    lp.viewClusterExtended();
//    lp.viewSiteByRelativeIndex();
    while (true){
        successful = lp.occupy();
        if(successful) {
//            cout << "step " << counter << " *************************************" << endl;
//            lp.ckeckCluster();
//            lp.viewClusterExtended();
//            lp.viewSiteByRelativeIndex();
            if(lp.detectWrapping()) {
                cout << "wrapping. pc = " << lp.occupationProbability() << endl;
//                lp.viewSiteByRelativeIndex();
//                cout << lp.wrapping_indices()[0] << " and " << lp.wrapping_indices()[1] << endl;
////                vector<double> tmp = lp.clusterSizeDistribution();
//////                    cout << "after returned " << tmp.size() << endl;
////                if (tmp.size() > cluster_size_distro.size()) {
////                    cluster_size_distro.resize(tmp.size());
////                }
////                for (size_t k{}; k < tmp.size(); ++k) {
////                    cluster_size_distro[k] += tmp[k];
////                }
////                lp.viewClusterExtended();
////                lp.viewLatticeByID();
                break;
            }
//            break;
            ++counter;
        }
//        break;
//        if(counter == 15) break;
        if(counter >= lp.maxIterationLimit()){ // twice_length_squared is the number of bonds
            break;
        }
    }

//    lp.viewClusterExtended();
//    lp.viewSiteByRelativeIndex();
}

/**
 * Site percolation by placing sites
 */
void test_site_percolation(value_type length, value_type ensemble_size) {
    cout << "length " << length << " ensemble_size " << ensemble_size << endl;

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    vector<double> cluster_size_distro;
    double area=0;
    SitePercolationBallisticDeposition_L2_v2 lp(length, true);
    lp.setRandomState(0, true);
    lp.init();
//    lp.reset();

    bool successful = false;
    auto t_start = std::chrono::system_clock::now();
    size_t counter = 0;
    bool detect_wrapping{true};
//        cout << "line " << __LINE__ << endl;
//    lp.ckeckCluster();
//    lp.viewClusterExtended();
//    lp.viewSiteByRelativeIndex();

    while (true){
        successful = lp.occupy();
        if(successful) {
            cout << "step " << counter << " *************************************" << lp.lastPlacedSite() << endl;
//            lp.ckeckCluster();
//            lp.viewClusterExtended();
//            lp.viewSiteByID();
            lp.viewSiteByRelativeIndex();
            if(detect_wrapping && lp.detectWrapping()) {
                cout << "wrapping. pc = " << lp.occupationProbability() << endl;
                detect_wrapping = false;
//                cout << lp.wrapping_indices()[0] << " and " << lp.wrapping_indices()[1] << endl;
//                vector<double> tmp = lp.clusterSizeDistribution();
//                    cout << "after returned " << tmp.size() << endl;
//                if (tmp.size() > cluster_size_distro.size()) {
//                    cluster_size_distro.resize(tmp.size());
//                }
//                for (size_t k{}; k < tmp.size(); ++k) {
//                    cluster_size_distro[k] += tmp[k];
//                }
//                lp.viewClusterExtended();
//                lp.viewLatticeByID();
                break;
            }
            ++counter;
        }
//        break;
//        if(counter == 3) break;
        if(counter >= lp.maxIterationLimit()){ // twice_length_squared is the number of bonds
            break;
        }
    }

    lp.viewClusterExtended();
    lp.viewSiteByRelativeIndex();
}


/**
 * Test Site Percolation Program by generating data using different methods.
 * Current test class : SitePercolation_ps_v10
 */
void test_simulate_site_percolation(value_type length, value_type ensemble_size) {

    std::cout << "length " << length << " ensemble_size " << ensemble_size << std::endl;

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    SitePercolation_ps_v10 lattice_percolation(length, true);
    lattice_percolation.setRandomState(0, true);
    lattice_percolation.init();

//    std::string tm = getCurrentTime();
//    const char delimiter = '\t';
//    std::ostringstream header_info;
//    header_info << "{"
//                << R"("length":)" << length
//                << R"(,"ensemble_size":)" << ensemble_size
//                << R"(,"random_seed":)" << lattice_percolation.getRandomState()
//                << R"(,"signature":")" << lattice_percolation.getSignature() << "\""
//                << R"(,"time":")" << tm << "\""
//                //                << R"(,"cols":)" << R"(["S","n_S"])"
//                << "}";
//
//
//    std::string extension = "_L" + std::to_string(length) + '_' + tm + ".txt";
////    std::string filename_s = lattice_percolation.getSignature() + "_cluster_by_site" + extension;
////    std::string filename_b = lattice_percolation.getSignature() + "_cluster_by_bond" + extension;
//    std::string filename_critical = lattice_percolation.getSignature() + "_critical" + extension;
////    std::string filename = lattice_percolation.getSignature() + "_entropy-jump" + extension;
//    std::string filename_entropy_order_parameter = lattice_percolation.getSignature() + extension;


    // simulation starts here
    std::vector<double> pcs(ensemble_size), sites_pc(ensemble_size), bonds_pc(ensemble_size);
    value_type counter{};
    std::vector<double> entropy(lattice_percolation.maxIterationLimit());
    std::vector<double> nob_wraping(lattice_percolation.maxIterationLimit()),
            nob_largest(lattice_percolation.maxIterationLimit());

    for(value_type i{} ; i != ensemble_size ; ++i){
        lattice_percolation.reset();
        cout << "**************************************************** " << i << endl;
        lattice_percolation.viewClusterExtended();
        lattice_percolation.viewLatticeByID();

        bool successful = false;
        auto t_start = std::chrono::system_clock::now();
        counter = 0;
        bool detect_wrapping{true};
        while (true){
            successful = lattice_percolation.occupy();
            cout << lattice_percolation.lastPlacedSite() << endl;
            lattice_percolation.viewCluster(lattice_percolation.lastPlacedSite());
            lattice_percolation.viewLatticeByID();

            if(successful) {
//                double H1 = lattice_percolation.entropy_v1();
//                double H2 = lattice_percolation.entropy_v2();
//                if(abs(H1 - H2) > 1e-10){
//                    cout << "entropy differs much" << endl;
//                }

//                nob_wraping[counter] += lattice_percolation.numberOfBondsInTheWrappingClusters();
//                nob_largest[counter] += lattice_percolation.numberOfBondsInTheLargestCluster();
//                lattice_percolation.jump();
//                if(detect_wrapping && lattice_percolation.detectWrapping()){
//                    pcs[i] = lattice_percolation.occupationProbability();
//                    sites_pc[i] = lattice_percolation.numberOfSitesInTheWrappingClusters();
//                    bonds_pc[i] = lattice_percolation.numberOfBondsInTheWrappingClusters();
//
//                    detect_wrapping = false;
//                }
//
                ++counter;
            }
            if(counter >= lattice_percolation.maxIterationLimit()){ // twice_length_squared is the number of bonds
                break;
            }
        }
//        fout_jump << lattice_percolation.largestEntropyJump() << "," << lattice_percolation.largestEntropyJump_pc() << std::endl;
        {
            auto t_end = std::chrono::system_clock::now();
            std::cout << "Iteration " << i
                      //                 << " . Thread " << std::this_thread::get_id()
                      << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << std::endl;
        }
//        cout << "Relabeling time " << lattice_percolation.get_relabeling_time() << endl;
    }

//    fout_b.close();
//    fout_s.close();
//    fout_jump.close();
//
//    std::ofstream fout_critical(filename_critical);
//    fout_critical << '#' << header_info.str() << std::endl;
//    fout_critical << "#data at critical occupation probability or pc" << std::endl;
//    fout_critical << "#<pc>,<sites in wrapping cluster>,<bonds in wrapping cluster>" << std::endl;
//    for(size_t i{}; i < ensemble_size; ++i){
//        if(pcs[i] == 0) continue;
//        fout_critical << pcs[i] << delimiter
//                      << sites_pc[i] << delimiter
//                      << bonds_pc[i] << std::endl;
//    }
//    fout_critical.close();
//
//
//
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
//    for(size_t i{}; i < lattice_percolation.maxIterationLimit(); ++i){
//
//        fout << (i + 1) / double(lattice_percolation.maxIterationLimit()) << delimiter;
//        fout << entropy[i] / double(ensemble_size) << delimiter;
//        fout << nob_largest[i] / double(ensemble_size * lattice_percolation.maxBonds()) << delimiter;
//        fout << nob_wraping[i] / double(ensemble_size * lattice_percolation.maxBonds()) ;
//        fout << std::endl;
//    }
//    fout.close();
}
