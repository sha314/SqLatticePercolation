//
// Created by shahnoor on 11/12/19.
//

#ifndef SQLATTICEPERCOLATION_TEST_V10_H
#define SQLATTICEPERCOLATION_TEST_V10_H

#include <iostream>
#include <vector>
#include <chrono>
#include <sstream>
#include <fstream>
#include "../util/time_tracking.h"
#include "../types.h"

void test_v10(int argc, char* argv[]);
void test_v10_lattice(int argc, char* argv[]);
void test_v10_number_of_clusters_per_site(int argc, char **argv);

void run_v10_regular(int length, int ensemble_size);

template<class PType>
void run_v10_rsbd(int length, int ensemble_size){
    std::cout << "length " << length << " ensemble_size " << ensemble_size << std::endl;

    size_t length_squared = length*length;
    size_t twice_length_squared = 2 * length_squared;

    PType lattice_percolation(length, true);
    lattice_percolation.setRandomState(0, true);
    lattice_percolation.init();




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

    std::string tm = getCurrentTime();
    const char delimiter = '\t';
    std::ostringstream header_info;
    header_info << "{"
                << R"("length":)" << length
                << R"(,"ensemble_size":)" << ensemble_size
                << R"(,"random_seed":)" << lattice_percolation.getRandomState()
                << R"(,"signature":")" << lattice_percolation.getSignature() << "\""
                << R"(,"datetime":")" << tm << "\""
                << R"(,"classname":")" << lattice_percolation.getClassName() << "\""
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

    std::ofstream fout_critical(filename_critical);
    fout_critical << '#' << header_info.str() << std::endl;
    fout_critical << "#data at critical occupation probability or pc" << std::endl;
    fout_critical << "#<pc>,<sites in wrapping cluster>,<bonds in wrapping cluster>" << std::endl;
    for(size_t i{}; i < ensemble_size; ++i){
        if(pcs[i] == 0) continue;
        fout_critical << pcs[i] << delimiter
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
    for(size_t i{}; i < lattice_percolation.maxIterationLimit(); ++i){

        fout << (i + 1) / double(lattice_percolation.maxIterationLimit()) << delimiter;
        fout << entropy[i] / double(ensemble_size) << delimiter;
        fout << nob_largest[i] / double(ensemble_size * lattice_percolation.maxBonds()) << delimiter;
        fout << nob_wraping[i] / double(ensemble_size * lattice_percolation.maxBonds()) ;
        fout << std::endl;
    }
    fout.close();
}

#endif //SQLATTICEPERCOLATION_TEST_V10_H
