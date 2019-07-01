//
// Created by shahnoor on 9/30/2017.
//

#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>

#include "cmd_args.h"
#include "../percolation/percolation.h"
#include "../types.h"
#include "../util/time_tracking.h"
#include "../percolation/site/percolation_site_v9.h"
#include "../percolation/site/percolation_site_explosive.h"
#include "../percolation/site/percolation_site_ballistic.h"


using namespace std;


/**
 *
 */
void print_help(){
    std::string hlp = R"***(
Percolation on a Square Lattice

Available Argument

-h      Help
-l      Length of the lattice. Default value is 5
-r      Range of length. No default value
-n      Ensemble size. Default value is 1
-i      Number of impure site_index_sequence. Default value is 0
-m      only counted for explosive percolation. default value is 0
        0   ->  sum rule
        1   ->  product rule

-p      place site_index_sequence for
            0   ->  Just place all the site_index_sequence and printout the result in the commandline only once to see
                    what's going on there and how much time it requires.

            1   ->  Place site_index_sequence until first spanning. Without periodicity.

            2   ->  MicroCanonical Ensemble data of Site percolation
                        Number of occupied site_index_sequence,
                        Entropy,
                        number of bonds in the largest cluster
                    Other quantities can be calculated from these data
                    This option writes all data of every realization separately
                    option '4' does the same thing but writes only the ensemble average.

            3   ->  Normal percolation
                    Data of the following quantity at threshold, i.e., when spanning occurs.
                        p_c = critical occupation probability
                        id_sc = id of the spanning cluster
                        b_t_sc = birth time of the spanning cluster
                        sites_sc = number of site_index_sequence of the spanning cluster
                        bonds_sc = number of bonds of the spanning cluster

            4   ->  MicroCanonical Ensemble data of Site percolation
                        Number of occupied site_index_sequence,
                        Entropy,
                        number of bonds in the largest cluster
                    writes only the ensemble average of the above quantity.

            5   ->  explosive percolation
                    Data of the following quantity at threshold, i.e., when spanning occurs.
                        p_c = critical occupation probability
                        id_sc = id of the spanning cluster
                        b_t_sc = birth time of the spanning cluster
                        sites_sc = number of site_index_sequence of the spanning cluster
                        bonds_sc = number of bonds of the spanning cluster

            6   ->  MicroCanonical Ensemble data of explosive Site percolation
                        Number of occupied site_index_sequence,
                        Entropy,
                        number of bonds in the largest cluster
                    writes only the ensemble average of the above quantity.

            7   ->  MicroCanonical Ensemble data of Site percolation
                        Occupation Probability,
                        Temperature
                        Entropy,
                        Order Parameter,
                        Specific heat
                        Susceptibility

                    writes only the ensemble average of the above quantity.

            8   ->  MicroCanonical Ensemble data of explosive Site percolation
                        Occupation Probability,
                        Temperature
                        Entropy,
                        Order Parameter,
                        Specific heat
                        Susceptibility

                    writes only the ensemble average of the above quantity.

            9   ->  Box counting data for global fractal dimension and spanning fractal dimension

            10  ->  Cluster_v3 length data by bonds and site_index_sequence
    )***";
    std::cout << hlp << std::endl;
}




/**
 *
 */
void print_help_ballistic(){
    std::string hlp = R"***(
Percolation on a Square Lattice (Ballistic Deposition)
L0 regular site percolation
L1 upto first nearest neighbor
L2 upto second nearest neighbor

Available Argument

-h      Help
-l      Length of the lattice. Default value is 5
-n      Ensemble size. Default value is 1

-p      place site_index_sequence for
            0   ->  Just place all the site_index_sequence and printout the result in the commandline only once to see
                    what's going on there and how much time it requires.

            1   ->  Critical Point data (L1)
                    Data of the following quantity at threshold, i.e., when spanning occurs.
                        p_c = critical occupation probability
                        id_sc = id of the spanning cluster
                        b_t_sc = birth time of the spanning cluster
                        sites_sc = number of site_index_sequence of the spanning cluster
                        bonds_sc = number of bonds of the spanning cluster

            2   ->  MicroCanonical Ensemble data of Site percolation Ballistic Deposition (Periodic)(L1)
                        Occupation Probability,
                        Temperature
                        Entropy= sum u_i log(u_i),
                        Order Parameter=(bonds in the largest cluster)/(total bonds),
                        Specific heat
                        Susceptibility

            3   ->  MicroCanonical Ensemble data of Site percolation Ballistic Deposition (Non-Periodic)(L1)
                        Occupation Probability,
                        Temperature
                        Entropy= sum u_i log(u_i),
                        Order Parameter=(bonds in the largest cluster)/(total bonds),
                        Specific heat
                        Susceptibility

    )***";
    std::cout << hlp << std::endl;
}



/**
 * Data required to calculate entropy, specific heat
 * @param length
 * @param ensembleSize
 */
void percolation_data_all_site_placed(value_type length, value_type ensembleSize, value_type impure_site){
    SitePercolation_ps_v9 sp(length, true);
    cout << "Signature " << sp.getSignature() << endl;
    // setting up the calculation flag
//    sp.calculationFlags({CalculationFlags::Occupation_Probability, CalculationFlags::Entropy});

    clock_t t;

    string common;
    common += to_string(length);
    common += '_';
    common += currentTime();;
    common += ".txt";

    string bond_filename = "percolation_data_L_" + common;

    ofstream fout(bond_filename);        // the output data file

    fout << sp.getSignature() << endl;
    // here only minimum data is printed so that it would take smaller size on the storage
    fout << "#Percolation data" << endl;
    fout << "#Number of occupied site_index_sequence, Entropy, number of bonds in the largest cluster" << endl;
    fout << "#length = " << endl << length << endl;
    fout << "#Using the information above we can calculate the followings" << endl;
    fout << "#Occupation probability = (Number of occupied site_index_sequence) / (length * length)" << endl;
    fout << "#Order parameter = (number of bonds in the largest cluster) / (2 * length * length)" << endl;
    fout << "#Format of the data" << endl;
    fout << "#\'>\' symbol indicates begining of a completely new system" << endl;
    fout << "#Number of \'>\' is the ensemble size in the file" << endl;
    fout << "#first column is printed as integer so that we can easily take ensemble average form given data" << endl;
    fout << "#<Number of Occupied Sites>\t<Entropy>\t<number of bonds in the largest cluster>" << endl;

    for(value_type i{} ; i != ensembleSize ; ++i){

        t = clock();
        sp.reset();

//        sp.placeAllSites();
        fout << ">" << endl;
        while (sp.occupy()){

            double nos = sp.numberOfOccupiedSite();
//            value_type nob1 = sp.numberOfBondsInTheLargestCluster();
            value_type nob2 = sp.numberOfBondsInTheLargestCluster_v2();
//            double H1 = sp.entropy();
//            double H2 = sp.entropy_v2();
            double H3 = sp.entropy();
//            if(abs(H2 - H3) > 1e-5){
//                cout << "Newly calculated entropy values are not equal" << endl;
//            }
            fout << nos << '\t'
//                 << H1 << '\t'
//                 << H2 << '\t'
                 << H3 << '\t'
//                 << nob1 << '\t'
                 << nob2
                 << endl;

//            if((H1-H2) > 1e-3){
//                cout << "H1 != H2" << endl;
//            }
//            if(nob1 != nob2){
////                cout << "nob1 != nob2" << endl;
//                cout << nob1 << " != " << nob2 << endl;
//            }
//            if(nos == length * length && nob2 != 2 * length * length){
//                cout << "error : line " << endl;
//            }
        }

        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    fout.close();

}



/**
 *
 * @param length
 * @param ensembleSize
 * @param rule
 * @param impure_site
 */
void percolation_data_sq_lattice(value_type length, value_type ensembleSize){
    SitePercolation_ps_v9 sp(length, true);  // periodicity is turned on

    cout << "Signature " << sp.signature << endl;
    // setting up the calculation flag
//    sp.calculationFlags({CalculationFlags::Occupation_Probability, CalculationFlags::Entropy});

    clock_t t;
    size_t length_squared = length*length;
    std::vector<double> nos(length_squared), nob(length_squared), entrpy(length_squared);
    size_t j{};
    for(value_type i{} ; i != ensembleSize ; ++i){

        t = clock();
        sp.reset();
        j = 0;
        while (sp.occupy()){
            nos[j] += sp.numberOfOccupiedSite();
            nob[j] += sp.numberOfBondsInTheLargestCluster_v2();
            entrpy[j] += sp.entropy();  // old method and takes long time
//            entrpy[j] += sp.entropy_v3(); // faster method
            ++j;
        }
//        cout << "back " <<  nos.back() << " " << nob.back() << endl;
        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    // Normalizing
    for(size_t i{}; i!= length_squared ; ++i){
        nos[i] /= double(ensembleSize);
        entrpy[i] /= double(ensembleSize);
        nob[i] /= double(ensembleSize);
    }

    // calculating and writing to a file
    string filename = "sq_lattice_site_percolation_" + to_string(length) + "_calculated_" + currentTime();
    filename += ".txt";
    calculate_and_write_to_file(length, ensembleSize, 2*length_squared, nob, entrpy, filename);

}




/**
 *
 * @param length
 * @param ensembleSize
 * @param rule
 * @param impure_site
 */
void percolation_data_sq_lattice_ballistic_L1(bool periodicity, value_type length, value_type ensembleSize){
    SitePercolationBallisticDeposition_v2 sp(length, periodicity);

    cout << "Signature " << sp.getSignature() << endl;
    // setting up the calculation flag
//    sp.calculationFlags({CalculationFlags::Occupation_Probability, CalculationFlags::Entropy});

    clock_t t;

    size_t length_squared = length*length;
    std::vector<double> nos(length_squared), nob_def1(length_squared), entrpy(length_squared);
    size_t j{};
    bool flag{};
    for(value_type i{} ; i != ensembleSize ; ++i){

        t = clock();
        sp.reset();
        j = 0;
        while (j < length_squared){
            flag = sp.occupy();
            if(flag) {
                nos[j] += sp.numberOfOccupiedSite();
                nob_def1[j] += sp.numberOfBondsInTheLargestCluster_v2(); // is this or
//                nob_def1[j] += sp.numberOfSitesInTheLargestCluster(); // this
//                entrpy[j] += sp.entropy();  // old method and takes long time
                entrpy[j] += sp.entropy(); // faster method
//                cout << nos[j] << "  " << nob_def1[j] << endl;

                ++j;

            }

        }
//        cout << "back " <<  nos.back() << " " << nob_def1.back() << endl;

        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    // Normalizing
    for(size_t i{}; i!= length_squared ; ++i){
        nos[i] /= double(ensembleSize);
        entrpy[i] /= double(ensembleSize);
        nob_def1[i] /= double(ensembleSize);
    }

    cout << "last " << nob_def1.back() << endl;
    // calculating and writing to a file
    string filename = "sq_lattice_site_percolation_ballistic_deposition_L1_";
    if (periodicity) {
        filename += to_string(length) + "_calculated_" + currentTime();
    }    else{
        filename += "non-periodic_"+to_string(length) + "_calculated_" + currentTime();
    }
    filename += ".txt";
    calculate_and_write_to_file(length, ensembleSize, 2*length_squared, nob_def1, entrpy, filename);

}


/**
 *
 * @param length
 * @param ensembleSize
 * @param length_squared
 * @param nob
 * @param entrpy
 * @param filename
 */
void calculate_and_write_to_file(value_type length, value_type ensembleSize, value_type total_number_of_bonds,
                                 const vector<double> &nob, const vector<double> &entrpy, const string &filename) {
    size_t length_squared = length*length;
    ofstream fout(filename);
    string column_names = "#<p>\t<T>\t<H(p,L)>\t<P<p,L>>\t<C(p,L)>\t<X(p,L)>\n";
    fout << column_names;

    string details = "#p = occupation probability\n";
    details += "#T = temperature = 1-p \n";
    details += "#H(p,L) = Entropy = sum( - u_i * log(u_i))\n";
    details += "#P(p,L) = Order parameter = (number of bonds in largest cluster) / (total number of bonds)\n";
    details += "#C(p,L) = Specific heat = -T dH/dT\n";
    details += "#X(p,L) = Susceptibility = dP/dp \n";
    details += "#u_i = (number of bonds in the i-th cluster) / (total number of bonds)\n";
    fout << details << endl;

    string header = "BEGIN_HEADER\n";
    header += "length\t"+ to_string(length) +"\n";
    header += "ensemble_size\t"+ to_string(ensembleSize) +"\n";
    header += "data_line\t16\n";   // todo change if neccessary;
    header += "END_HEADER\n";
    fout << header << endl;

    double p, H, P, C{}, X{};
    for(size_t i{}; i!= length_squared; ++i){
        p = (i+1)/ double(length_squared);
        H = entrpy[i];                      // Entropy
        P = nob[i] / double(total_number_of_bonds);    // Order parameter

        if(i != 0)
            C = -1 * (1-p) * (entrpy[i] - entrpy[i-1]) * length_squared;    //  minus sign is logical and to get upward curve
        if(i != 0)
            X = (nob[i] - nob[i-1])/2;          // susceptibility

        fout << p << '\t' << (1-p) << '\t' << H << '\t' << P << '\t' << C << '\t' << X << endl;
    }
    fout.close();
}



/**
 * Data required to calculate entropy, specific heat
 * @param length
 * @param ensembleSize
 */
void entropy_order_parameter(value_type length, value_type ensembleSize, value_type impure_site){
    SitePercolation_ps_v9 sp(length, true);

    cout << "Signature " << sp.getSignature() << endl;
    // setting up the calculation flag
//    sp.calculationFlags({CalculationFlags::Occupation_Probability, CalculationFlags::Entropy});

    clock_t t;
    std::vector<double> nos(length*length), nob(length*length), entrpy(length*length);
    size_t j{};
    for(value_type i{} ; i != ensembleSize ; ++i){

        t = clock();
        sp.reset();
        j = 0;
        while (sp.occupy()){
            nos[j] += sp.numberOfOccupiedSite();
            nob[j] += sp.numberOfBondsInTheLargestCluster_v2();
            entrpy[j] += sp.entropy();
            ++j;
        }

        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }


    string common;
    common += to_string(length);
    common += '_';
    common += currentTime();;
    common += ".txt";

    string bond_filename = "entropy_order_parameter_L_" + common;

    ofstream fout(bond_filename);        // the output data file

    fout << sp.signature << endl;
    // here only minimum data is printed so that it would take smaller size on the storage
    fout << "#<Number of Occupied Sites>\t<Entropy>\t<number of bonds in the largest cluster>" << endl << endl;
    fout << "BEGIN_HEADER" << endl;
    fout << "ensemble size " << ensembleSize << endl;
    fout << "length " << length << endl;
    fout << "data_line " << 12 << endl; // todo
    fout << "END_HEADER" << endl;
    fout << "#Percolation on a square lattice" << endl;
    fout << "#Number of occupied site_index_sequence, Entropy, number of bonds in the largest cluster" << endl;
    fout << "#Using the information above we can calculate the followings" << endl;
    fout << "#Occupation probability = (Number of occupied site_index_sequence) / (length * length)" << endl;
    fout << "#Order parameter = (number of bonds in the largest cluster) / (2 * length * length)" << endl;
    fout << "#first column is printed as integer so that we can easily take ensemble average form given data" << endl;

    // writing the normalized data to file
    for(size_t i{}; i!=length*length ; ++i){
        fout << nos[i]/double(ensembleSize) << '\t'
             << entrpy[i]/double(ensembleSize) << '\t'
             << nob[i]/double(ensembleSize) << endl;
    }
    fout.close();

}



/**
 *
 * @param length
 * @param ensemble_size
 * @param impure_sites
 */
void box_counting_square_lattice(value_type length, value_type ensemble_size, value_type impure_sites){
    SitePercolation_ps_v9 sp(length, true); // todo breaks down the program // possible solution: remove redundant variables

    cout << "Signature " << sp.signature << endl;
    clock_t t;
    
    map<value_type, array<value_type, 2>> N_delta;
    // finding all delta to perform box counting
    array<value_type, 2> dummy;
    vector<value_type> delta;
    cout << "Delta :{";
    for(value_type i{1}; i != length; ++i){
        if(length % i == 0){
            delta.push_back(i);
            N_delta[i][0] = 0;
            N_delta[i][1] = 0;
            cout << i << ',';
        }
    }
    cout << "}" << endl;


    for(value_type i{} ; i != ensemble_size ; ++i){
        t = clock();
        sp.reset();

        while(!sp.detectWrapping()){
            sp.occupy();
        }


        for(auto d: delta){
            dummy = sp.box_counting_v2(d);
            N_delta[d][0] += dummy[0];
            N_delta[d][1] += dummy[1];
//            cout << "not defined " << __LINE__ << endl;
        }

        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    string common;
    common += to_string(length);
    common += '_';
    common += currentTime();;
    common += ".txt";

    string time_filename = "sq_lattice_box_counting_L_" + common;

    ofstream fout(time_filename);      // contains the time information of the spanning cluster


    /**
     * First uncomment line is the signature
     * Second uncomment line is the parameter value
     * Third and all other uncomment lines are data
     */

    fout << "#<delta>\t<N(delta) global>\t<N(delta> spanning>" << endl;
    fout << "BEGIN_HEADER" << endl;
    fout << "ensemble_size\t"  << ensemble_size << endl;
    fout << "length\t"         << length << endl;
    fout << "data_line\t"      << 11 << endl;
    fout << "END_HEADER" << endl;
    fout << "#Box counting with different size to get global fractal dimension" << endl;
    fout << "#delta : length of the box" << endl;
    fout << "#N(delta) global: number of points found by box counting (global picture)" << endl;
    fout << "#N(delta) spanning: number of points found by box counting (only the spanning cluster)" << endl;

    cout << "Used Delta :{";
    for(auto d: delta){
        cout << d << ',';
        fout << d << '\t'
             << N_delta[d][0] / double(ensemble_size) << '\t'
             << N_delta[d][1] / double(ensemble_size)<< endl;
    }
    cout << "}" << endl;

    fout.close();
    cout << "file is closed " << endl;
}






