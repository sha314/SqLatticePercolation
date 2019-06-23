//
// Created by shahnoor on 04-May-18.
//

#include "../percolation/percolation.h"
#include "site_percolation.h"
#include "../util/time_tracking.h"
#include "../cmdLineArgs/cmd_args.h"
#include <fstream>

using namespace std;

template <class TPercolation>
void sq_lattice_site_percolation_periodic(size_t length, size_t ensemble_size){
    // todo
    std::cout << "Entry \"sq_lattice_site_percolation_periodic<T>\" : line " << __LINE__ << std::endl;
    std::cout << "On test" << std::endl;
    TPercolation percolation(length, true);
    std::cout << "Signature " << percolation.signature << std::endl;

}



void sq_lattice_site_percolation_periodic(size_t length, size_t ensemble_size){
    SitePercolation_ps_v9 sp(length, true);

    cout << "Signature " << sp.getSignature() << endl;
    clock_t t;

    size_t length_squared = length*length;
    std::vector<double> nos(length_squared), nob_def1(length_squared), nob_def2(length_squared), entrpy(length_squared);
    size_t j{};
    bool flag{};

    string common;
    common += to_string(length);
    common += '_';
    common += currentTime();;
    common += ".txt";

    string filename = "critical_point_data_sq_lattice_sp_periodic_L_" + common;

    ofstream fout_wrapping_cluster(filename);      // contains the time information of the spanning cluster

    /**
     * First uncomment line is the signature
     * Second uncomment line is the parameter value
     * Third and all other uncomment lines are data
     */

    fout_wrapping_cluster << "#<p_c>\t<id_sc>\t<b_t_sc>\t<sites_sc>\t<bonds_sc>" << endl;
    fout_wrapping_cluster << "BEGIN_HEADER" << endl;
    fout_wrapping_cluster << "ensemble_size\t"  << ensemble_size << endl;
    fout_wrapping_cluster << "length\t"         << length << endl;
    fout_wrapping_cluster << "data_line\t"      << 11 << endl;
    fout_wrapping_cluster << "END_HEADER" << endl;
    fout_wrapping_cluster << "#p_c = critical occupation probability" << endl;
    fout_wrapping_cluster << "#id_sc = id of the wrapping cluster" << endl;
    fout_wrapping_cluster << "#b_t_sc = birth time of the wrapping cluster" << endl;
    fout_wrapping_cluster << "#sites_sc = number of site_index_sequence of the wrapping cluster" << endl;
    fout_wrapping_cluster << "#bonds_sc = number of bonds of the wrapping cluster" << endl;


    for(value_type i{} ; i != ensemble_size ; ++i){
        t = clock();
        sp.reset();

        bool wrapping_occured {false};
        while (j < length_squared){
            flag = sp.occupy();
            if(flag) {
                ++j;
                nos[j] += sp.numberOfOccupiedSite();
                nob_def1[j] += sp.numberOfBondsInTheLargestCluster_v2(); // is this or
                nob_def2[j] += sp.numberOfSitesInTheSpanningClusters_v2(); // this
//                entrpy[j] += sp.entropy();  // old method and takes long time
                entrpy[j] += sp.entropy(); // faster method
            }
            if(flag && sp.detectWrapping()){
                if(!wrapping_occured) {
                    wrapping_occured = true;

                    fout_wrapping_cluster << sp.occupationProbability() << '\t'
                                          << sp.birthTimeOfSpanningCluster() << '\t'
                                          << sp.numberOfSitesInTheSpanningClusters_v2() << '\t'
                                          << sp.numberOfBondsInTheSpanningClusters_v2() << endl;
                }
            }
        }


        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    cout << "reached the end " << endl;
    fout_wrapping_cluster.close();
    cout << "file is closed " << endl;

    // Normalizing
    for(size_t i{}; i!= length_squared ; ++i){
        nos[i] /= double(ensemble_size);
        entrpy[i] /= double(ensemble_size);
        nob_def1[i] /= double(ensemble_size);
        nob_def2[i] /= double(ensemble_size);
    }

    cout << "last " << nob_def1.back()  << " : line " << __LINE__ << endl;
    // calculating and writing to a file
    string filename2 = "sq_lattice_sp_periodic_L_" + common;

    calculate_and_write_to_file_v2(length, ensemble_size, 2*length_squared, nob_def1, nob_def2, entrpy, filename2);
}




void calculate_and_write_to_file_v2(value_type length, value_type ensembleSize, value_type total_number_of_bonds,
                                    const std::vector<double> &nob1, const std::vector<double> &nob2,
                                    const vector<double> &entrpy, const string &filename) {
    size_t length_squared = length*length;
    ofstream fout(filename);
    string column_names = "#<p>\t<H(p,L)>\t<P1<p,L>>\t<P2<p,L>>\n";
    fout << column_names;

    string details = "#p = occupation probability\n";
    details += "#H(p,L) = Entropy = sum( - u_i * log(u_i))\n";
    details += "#P1(p,L) = Order parameter = (number of bonds in the largest cluster) / (total number of bonds)\n";
    details += "#u_i = (number of bonds in the i-th cluster) / (total number of bonds)\n";
    fout << details << endl;

    string header = "BEGIN_HEADER\n";
    header += "length\t"+ to_string(length) +"\n";
    header += "ensemble_size\t"+ to_string(ensembleSize) +"\n";
    header += "data_line\t10\n";   // todo change if neccessary;
    header += "END_HEADER\n";
    fout << header << endl;

    double p, H, P1, P2;
    for(size_t i{}; i!= length_squared; ++i){
        p = (i+1)/ double(length_squared); // occupation probability
        H = entrpy[i];                      // Entropy
        P1 = nob1[i] / double(total_number_of_bonds);    // Order parameter
        P2 = nob2[i] / double(total_number_of_bonds);    // Order parameter

        fout << p << '\t' << (1-p) << '\t' << H << '\t' << P1 << '\t' << P2 << endl;
    }
    fout.close();
}