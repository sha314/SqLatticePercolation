//
// Created by shahnoor on 9/20/18.
//

#ifndef SQLATTICEPERCOLATION_TEST_PERCOLATION_H
#define SQLATTICEPERCOLATION_TEST_PERCOLATION_H

#include <iostream>
#include <string>
#include <chrono>
#include "../types.h"
#include "../percolation/percolation.h"
#include "../util/time_tracking.h"

void visual(value_type );
void cluster_size_distribution_site(value_type length, value_type ensemble_size);
void cluster_size_distribution_bond(value_type length, value_type ensemble_size);
void simulate_site_percolation(value_type, value_type );
void simulate_bond_percolation(value_type, value_type );
void simulate_bond_percolation_v2(value_type length, value_type ensemble_size); // 2019.01.01
void simulate_site_percolation_detailed(value_type length, value_type ensemble_size);

void entropyJumps(int argc, char** argv);

void test_bond_percolation(value_type length, value_type ensemble_size);

/**
 * data for log(Mass) vs log(Length) curve
 * for spanning sites only
 */
void percolation_to_get_fractal(
        value_type ensemble_size, value_type L_start, value_type l_end, value_type delta_L=20
);

/**
 * data for log(Mass) vs log(Length) curve
 * for spanning sites only
 */
void percolation_fractalDimension_by_spanning_site(
        value_type ensemble_size, value_type L_start, value_type l_end, value_type delta_L=20
);
void measure_entropy_by_site(int argc, char** argv);

template<class PType>
void simulate_site_percolation_T(value_type length, value_type ensemble_size) {

    std::cout << "length " << length << " ensemble_size " << ensemble_size << std::endl;

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    PType lattice_percolation(length, true);

    std::ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length
                << ",\"ensemble_size\":" << ensemble_size
                << ",\"signature\":\"" << lattice_percolation.getSignature() << "\""
                << "}" ;

    std::string tm = getCurrentTime();
    std::string extension = "_L" + std::to_string(length) + '_' + tm + ".csv";
    std::string filename_s = lattice_percolation.getSignature() + "_cluster_by_site" + extension;
    std::string filename_b = lattice_percolation.getSignature() + "_cluster_by_bond" + extension;
    std::string filename_critical = lattice_percolation.getSignature() + "_critical" + extension;
    std::string filename = lattice_percolation.getSignature() + "_entropy-jump" + extension;
    std::string filename_entropy_order_parameter = lattice_percolation.getSignature() + extension;


    std::ofstream fout_jump(filename);
    // JSON formated header
    fout_jump << '#' << header_info.str() << std::endl;
    fout_jump << "#each line is an independent realization" << std::endl;
    fout_jump << "#each line contains information about all clusters at critical point" << std::endl;
    fout_jump << "#cluster size is measured by number of bonds in it" << std::endl;
    fout_jump << "#<entropy jump>,<occupation probability at which the jump occurs>" << std::endl;

    std::ofstream fout_s(filename_s);
    // JSON formated header
    fout_s << '#' << header_info.str() << std::endl;
    fout_s << "#each line is an independent realization" << std::endl;
    fout_s << "#each line contains information about all clusters at critical point" << std::endl;
    fout_s << "#cluster size is measured by number of site_index_sequence in it" << std::endl;

    std::ofstream fout_b(filename_b);
    // JSON formated header
    fout_b << '#' << header_info.str() << std::endl;
    fout_b << "#each line is an independent realization" << std::endl;
    fout_b << "#each line contains information about all clusters at critical point" << std::endl;
    fout_b << "#cluster size is measured by number of bonds in it" << std::endl;

    std::ofstream fout_critical(filename_critical);
    fout_critical << '#' << header_info.str() << std::endl;
    fout_critical << "#data at critical occupation probability or pc" << std::endl;
    fout_critical << "#<pc>,<site_index_sequence in wrapping cluster>,<bonds in wrapping cluster>" << std::endl;

    // simulation starts here
    value_type counter{};
    std::vector<double> entropy(lattice_percolation.maxIterationLimit());
    std::vector<double> nob_wraping(lattice_percolation.maxIterationLimit()),
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
                nob_largest[counter] += lattice_percolation.numberOfBondsInTheLargestCluster_v2();
                lattice_percolation.jump();
                if(!wrapping_written && lattice_percolation.detectWrapping()){
                    fout_critical << lattice_percolation.occupationProbability() << ","
                                  << lattice_percolation.numberOfSitesInTheWrappingClusters() << ","
                                  << lattice_percolation.numberOfBondsInTheWrappingClusters()  << std::endl;

                    std::vector<value_type> site, bond;

                    lattice_percolation.get_cluster_info(site, bond);

                    for(value_type j{}; j != site.size(); ++j){
                        fout_s << site[j] << ',';
                    }
                    for(value_type j{}; j != bond.size(); ++j){
                        fout_b << bond[j] <<',';
                    }


                    fout_s << std::endl;
                    fout_b << std::endl;
                    wrapping_written = true;
                }


                ++counter;
            }
            if(counter >= lattice_percolation.maxIterationLimit()){ // twice_length_squared is the number of bonds
                break;
            }
        }
        fout_jump << lattice_percolation.largestEntropyJump() << "," << lattice_percolation.largestEntropyJump_pc() << std::endl;
        {
            auto t_end = std::chrono::system_clock::now();
            std::cout << "Iteration " << i
                 //                 << " . Thread " << std::this_thread::get_id()
                 << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << std::endl;
        }
//        cout << "Relabeling time " << lattice_percolation.get_relabeling_time() << endl;
    }

    fout_b.close();
    fout_s.close();
    fout_critical.close();
    fout_jump.close();



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
        fout << (i+1) / double(lattice_percolation.maxIterationLimit()) << ",";
        fout << entropy[i] / double(ensemble_size) << ",";
        fout << nob_largest[i] / double(ensemble_size /* lattice_percolation.maxBonds()*/) << ",";
        fout << nob_wraping[i] / double(ensemble_size /* lattice_percolation.maxBonds()*/) ;
        fout << std::endl;
    }
    fout.close();
}

#endif //SQLATTICEPERCOLATION_TEST_PERCOLATION_H
