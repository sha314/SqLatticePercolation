//
// Created by shahnoor on 9/30/2017.
//

#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>

#include "cmd_args.h"
#include "../types.h"
#include "percolation_fractal.h"



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
-i      Number of impure sites. Default value is 0
-m      only counted for explosive percolation. default value is 0
        0   ->  sum rule
        1   ->  product rule

-p      place sites for
            0   ->  Just place all the sites and printout the result in the commandline only once to see
                    what's going on there and how much time it requires.

            1   ->  Place sites until first spanning. Without periodicity.

            2   ->  MicroCanonical Ensemble data of Site percolation
                        Number of occupied sites,
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
                        sites_sc = number of sites of the spanning cluster
                        bonds_sc = number of bonds of the spanning cluster

            4   ->  MicroCanonical Ensemble data of Site percolation
                        Number of occupied sites,
                        Entropy,
                        number of bonds in the largest cluster
                    writes only the ensemble average of the above quantity.

            5   ->  explosive percolation
                    Data of the following quantity at threshold, i.e., when spanning occurs.
                        p_c = critical occupation probability
                        id_sc = id of the spanning cluster
                        b_t_sc = birth time of the spanning cluster
                        sites_sc = number of sites of the spanning cluster
                        bonds_sc = number of bonds of the spanning cluster

            6   ->  MicroCanonical Ensemble data of explosive Site percolation
                        Number of occupied sites,
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

            10  ->  Cluster length data by bonds and sites
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

-p      place sites for
            0   ->  Just place all the sites and printout the result in the commandline only once to see
                    what's going on there and how much time it requires.

            1   ->  Critical Point data (L1)
                    Data of the following quantity at threshold, i.e., when spanning occurs.
                        p_c = critical occupation probability
                        id_sc = id of the spanning cluster
                        b_t_sc = birth time of the spanning cluster
                        sites_sc = number of sites of the spanning cluster
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
    SitePercolation_ps_v8 sp(length, impure_site);
    sp.periodicityFlag(true);   // periodicity is turned on
    cout << "Signature " << sp.signature << endl;
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

    fout << sp.signature << endl;
    // here only minimum data is printed so that it would take smaller size on the storage
    fout << "#Percolation data" << endl;
    fout << "#Number of occupied sites, Entropy, number of bonds in the largest cluster" << endl;
    fout << "#length = " << endl << length << endl;
    fout << "#Using the information above we can calculate the followings" << endl;
    fout << "#Occupation probability = (Number of occupied sites) / (length * length)" << endl;
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
            double H3 = sp.entropy_v3();
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
 * Data required to calculate entropy, specific heat for explosive percolation
 * @param length
 * @param ensembleSize
 */
void entropy_order_parameter_explosive(value_type length, value_type ensembleSize, int rule, value_type impure_site){
    SitePercolation_ps_v8 sp(length, impure_site);
    sp.periodicityFlag(true);   // periodicity is turned on
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
        while (sp.placeSite_explosive(rule)){
            nos[j] += sp.numberOfOccupiedSite();
            nob[j] += sp.numberOfBondsInTheLargestCluster_v2();
            entrpy[j] += sp.entropy_v3();
            ++j;
        }

        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    // Normalizing
    for(size_t i{}; i!= length_squared ; ++i){
        nos[i] /= double(ensembleSize);
        entrpy[i] /= double(ensembleSize);
        nob[i] /= double(ensembleSize);
    }

    string common;
    common += to_string(length);
    common += '_';
    common += "explosive_";
    common += currentTime();
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
    fout << "#Explosive Percolation on a square lattice" << endl;
    fout << "#Number of occupied sites, Entropy, number of bonds in the largest cluster" << endl;
    fout << "#Using the information above we can calculate the followings" << endl;
    fout << "#Occupation probability = (Number of occupied sites) / (length * length)" << endl;
    fout << "#Order parameter = (number of bonds in the largest cluster) / (2 * length * length)" << endl;
    fout << "#first column is printed as integer so that we can easily take ensemble average form given data" << endl;

    // writing the normalized data to file
    for(size_t i{}; i!= length_squared ; ++i){
        fout << nos[i] << '\t'
             << entrpy[i] << '\t'
             << nob[i] << endl;
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
    SitePercolation_ps_v8 sp(length, true);  // periodicity is turned on

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
void percolation_data_sq_lattice_explosive(value_type length, value_type ensembleSize, int rule){
    SitePercolation_ps_v8 sp(length, true); // periodic

    cout << "Signature " << sp.signature << endl;
    // setting up the calculation flag
    // sp.calculationFlags({CalculationFlags::Occupation_Probability, CalculationFlags::Entropy});

    clock_t t;
    size_t length_squared = length*length;
    std::vector<double> nos(length_squared), nob(length_squared), entrpy(length_squared);
    size_t j{};
    for(value_type i{} ; i != ensembleSize ; ++i){

        t = clock();
        sp.reset();
        j = 0;
        while (sp.placeSite_explosive(rule)){
            nos[j] += sp.numberOfOccupiedSite();
            nob[j] += sp.numberOfBondsInTheLargestCluster_v2();
            entrpy[j] += sp.entropy_v3();
            ++j;
        }

        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    // Normalizing
    for(size_t i{}; i!= length_squared ; ++i){
        nos[i] /= double(ensembleSize);
        entrpy[i] /= double(ensembleSize);
        nob[i] /= double(ensembleSize);
    }

    // calculating and writing to a file
    string filename = "sq_lattice_site_percolation_explosive_"
                      + to_string(length) + "_calculated_" + currentTime();
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
    SitePercolationBallisticDeposition sp(length, periodicity);

    cout << "Signature " << sp.signature << endl;
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
                entrpy[j] += sp.entropy_v3(); // faster method
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
    SitePercolation_ps_v8 sp(length, impure_site);
    sp.periodicityFlag(true);   // periodicity is turned on
    cout << "Signature " << sp.signature << endl;
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
            entrpy[j] += sp.entropy_v3();
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
    fout << "#Number of occupied sites, Entropy, number of bonds in the largest cluster" << endl;
    fout << "#Using the information above we can calculate the followings" << endl;
    fout << "#Occupation probability = (Number of occupied sites) / (length * length)" << endl;
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
 * Does the work according to the Docs/readme.txt file
 * @param argc
 * @param argv
 */
void cmd_args_site_percolation(int argc, char **argv){
    if(argc == 1){
        print_help();
        exit(0);
    }

    value_type ensembleSize{1}, length{5};
    value_type impure_site{0};
    int place_sites_for{-1};
    value_type L_start{}, L_end{};
    int rule{0}; // 0-> sum rule, 1-> product rule

    char * filename=NULL;

    for(int i{1} ; i < argc;){
        std::cout << i << "\t" << argv[i] << std::endl;

        switch (argv[i][1]){
            case 'm':
                rule = atoi(argv[i+1]);
                i += 2;
                break;
            case 'n':
                // ensemble size
                ensembleSize = value_type(atoi(argv[i+1]));
                i += 2;
                break;

            case 'l':
                // length
                if(L_end != 0){
                    cout << "Range of length parameter is already set" << endl;
                    break;
                }
                length = value_type(atoi(argv[i+1]));
                L_start = length;
                i += 2;
                break;

            case 'r':
                // range of length
                if(length != 0){
                    cout << "Length parameter is already set" << endl;
                    break;
                }
                L_start = value_type(atoi(argv[i+1]));
                L_end = value_type(atoi(argv[i+2]));
                i += 3;
                break;

            case 'i':
                impure_site = value_type(atoi(argv[i+1]));
                i += 2;
                break;

            case 'p':
                place_sites_for = value_type(atoi(argv[i+1]));
                i += 2;
                break;

            case 'h':
                print_help();
                exit(0);

            default:
                print_help();
                exit(0);
        }
    }

    cout << "Performing calculations for" << endl;
    cout << "Length of the Square Lattice " << length << endl;
    cout << "Ensemble size " << ensembleSize << endl;

    if(place_sites_for == 0) {
        cout << "Placing all the sites" << endl;
    }else{
        cout << "Placing sites until the spanning" << endl;
    }

//    run_program(ensembleSize, length, how_many_site_to_place, step);

    if(length != 0){
        cout << "Length parameter is already set" << endl;
        L_end = L_start;
    }

    if(place_sites_for == 0) {
        clock_t t0 = clock();
        SitePercolation_ps_v8 sp(L_start, impure_site);
        cout << "Signature " << sp.signature << endl;
        sp.periodicityFlag(true);
        while (sp.occupy());
        cout << "Time elapsed : " << (clock() - t0) / double(CLOCKS_PER_SEC) << " sec" << endl;
        if (L_start <= 50) {
            sp.viewLatticeExtended();
            sp.viewClusterExtended();
            sp.viewCluster_id_index();
        }
    }else if(place_sites_for == 1) {
        // place sites until first spanning and generate all data possible
        // can generate several gigabytes of data
        spanning_data_for_lengths(L_start, L_end, ensembleSize, impure_site);
    } else if(place_sites_for == 2) {
        percolation_data_all_site_placed(length, ensembleSize, impure_site);
    }else if(place_sites_for == 3) {
        critical_point_data_sq_lattice_site_percolation(length, ensembleSize, impure_site);
    }else if(place_sites_for == 4) {
        entropy_order_parameter(length, ensembleSize, impure_site);
    }else if(place_sites_for == 5) {
        critical_point_data_explosive(length, ensembleSize, rule, impure_site);
    }else if(place_sites_for == 6) {
        entropy_order_parameter_explosive(length, ensembleSize, rule, impure_site);
    }else if(place_sites_for == 7) {
        percolation_data_sq_lattice(length, ensembleSize);
    }else if(place_sites_for == 8) {
        percolation_data_sq_lattice_explosive(length, ensembleSize, rule);
    }else if(place_sites_for == 9) {
        // global fractal dimension data
        // and number of bonds and number of sites on the clusters data at spanning
        box_counting_square_lattice(length, ensembleSize, impure_site);
    }else if(place_sites_for == 10) {
        // and number of bonds and number of sites on the clusters data at spanning
        cluster_length_sq_lattice(length, ensembleSize, impure_site);
    }
    else{
        cout << "process id argument mismatched : line " << __LINE__ << endl;
    }

}



/**
 * Does the work according to the Docs/readme.txt file
 * @param argc
 * @param argv
 */
void cmd_args_site_percolation_ballistic_deposition(int argc, char **argv){
    if(argc == 1){
        print_help_ballistic();
        exit(0);
    }

    value_type ensembleSize{1}, length{5};
    int place_sites_for{-1};
    value_type L_start{}, L_end{};
    int rule{0}; // 0-> sum rule, 1-> product rule

    char * filename=NULL;

    for(int i{1} ; i < argc;){
        std::cout << i << "\t" << argv[i] << std::endl;

        switch (argv[i][1]){
            case 'm':
                rule = atoi(argv[i+1]);
                i += 2;
                break;
            case 'n':
                // ensemble size
                ensembleSize = value_type(atoi(argv[i+1]));
                i += 2;
                break;

            case 'l':
                // length
                if(L_end != 0){
                    cout << "Range of length parameter is already set" << endl;
                    break;
                }
                length = value_type(atoi(argv[i+1]));
                L_start = length;
                i += 2;
                break;

            case 'r':
                // range of length
                if(length != 0){
                    cout << "Length parameter is already set" << endl;
                    break;
                }
                L_start = value_type(atoi(argv[i+1]));
                L_end = value_type(atoi(argv[i+2]));
                i += 3;
                break;

            case 'p':
                place_sites_for = value_type(atoi(argv[i+1]));
                i += 2;
                break;

            case 'h':
                print_help_ballistic();
                exit(0);

            default:
                print_help_ballistic();
                exit(0);
        }
    }

    cout << "Performing calculations for" << endl;
    cout << "Length of the Square Lattice " << length << endl;
    cout << "Ensemble size " << ensembleSize << endl;

    if(place_sites_for == 0) {
        cout << "Placing all the sites" << endl;
    }else{
        cout << "Placing sites until the spanning" << endl;
    }

//    run_program(ensembleSize, length, how_many_site_to_place, step);

    if(length != 0){
        cout << "Length parameter is already set" << endl;
        L_end = L_start;
    }

    if(place_sites_for == 0) {
        clock_t t0 = clock();
        SitePercolationBallisticDeposition sp(length, true);
        cout << "Signature " << sp.signature << endl;
        while (sp.occupy());
        cout << "Time elapsed : " << (clock() - t0) / double(CLOCKS_PER_SEC) << " sec" << endl;
        if (L_start <= 50) {
            sp.viewLatticeExtended();
            sp.viewClusterExtended();
            sp.viewCluster_id_index();
        }
    }else if(place_sites_for == 1) {
        critical_point_data_sq_l_sp_ballistic_deposition_L1(length, ensembleSize);
    }
    else if(place_sites_for == 2){
        percolation_data_sq_lattice_ballistic_L1(true, length, ensembleSize);//peridic
    }
    else if(place_sites_for == 3){
        percolation_data_sq_lattice_ballistic_L1(false, length, ensembleSize);// non-peridic
    }
    else{
        cout << "process id argument mismatched : line " << __LINE__ << endl;
    }

}



/**
 *
 * @param length
 * @param ensemble_size
 * @param impure_sites
 */
void box_counting_square_lattice(value_type length, value_type ensemble_size, value_type impure_sites){
    SitePercolation_ps_v8 sp(length, impure_sites); // todo breaks down the program // possible solution: remove redundant variables
    sp.periodicityFlag(false);

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

        sp.placeSiteUntilFirstSpanning_v2();

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



/**
 *
 * @param length
 * @param ensemble_size
 * @param impure_sites
 */
void cluster_length_sq_lattice(value_type length, value_type ensemble_size, value_type impure_sites){
    SitePercolation_ps_v8 sp(length, impure_sites); // todo breaks down the program // possible solution: remove redundant variables
    sp.periodicityFlag(false);

    cout << "Signature " << sp.signature << endl;
    clock_t t;
    vector<value_type> nos, nob;

    string common;
    common += to_string(length);
    common += '_';
    common += currentTime();;
    common += ".txt";

    string site_filename = "sq_lattice_cluster_length_sites_L_" + common;
    string bond_filename = "sq_lattice_cluster_length_bonds_L_" + common;

    ofstream fout_site(site_filename);        // the output data file
    ofstream fout_bond(bond_filename);        // the output data file

    /**
     * First uncomment line is the signature
     * Second uncomment line is the parameter value
     * Third and all other uncomment lines are data
     */
    fout_site << "#" << sp.signature << endl;
    fout_site << "#cluster length data in site percolation" << endl;
    fout_site << "#clustter length = number of sites" << endl;
    fout_site << "BEGIN_HEADER" << endl;
    fout_site << "ensemble_size\t"  << ensemble_size << endl;
    fout_site << "length\t"         << length << endl;
    fout_site << "data_line\t"      << 10 << endl;
    fout_site << "END_HEADER" << endl;
    fout_site << "#<data...>" << endl;

    fout_bond << "#" << sp.signature << endl;
    fout_bond << "#cluster length data in site percolation" << endl;
    fout_bond << "#clustter length = number of bonds" << endl;
    fout_bond << "BEGIN_HEADER" << endl;
    fout_bond << "ensemble_size\t"  << ensemble_size << endl;
    fout_bond << "length\t"         << length << endl;
    fout_bond << "data_line\t"      << 10 << endl;
    fout_bond << "END_HEADER" << endl;
    fout_bond << "#<data...>" << endl;


    map<value_type, value_type> N_delta;
    // finding all delta to perform box counting
    vector<value_type> delta;
    for(value_type i{1}; i != length; ++i){
        if(length % i == 0){
            delta.push_back(i);
            N_delta[i] = 0;
        }
    }

    for(value_type i{} ; i != ensemble_size ; ++i){
        t = clock();
        sp.reset();

        sp.placeSiteUntilFirstSpanning_v2();


        nos = sp.number_of_site_in_clusters();
        nob = sp.number_of_bonds_in_clusters();

        for(auto a: nos){
            fout_site << a << ' ';   // one space is enough to separate the data
        }

        for(auto b: nob){
            fout_bond << b << ' ';   // one space is enough to separate the data
        }

        fout_site << endl;
        fout_bond << endl;

        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    fout_site.close();
    fout_bond.close();
    cout << "file is closed " << endl;
}
/**
 * give cmd argument as
 * L start
 * L end
 * ensembleSize
 * @param argc
 * @param argv
 */
void cmd_moment_conservation(int argc, char** argv) {
    cout << "Moment conservation data" << endl;
    if(argc != 4){
        cout << "Not correct argument";
        cout << "L_start L_end ensembleSize" << endl;
        cout << "Order must be maintained" << endl;
        return;
    }
    value_type L_start = value_type(atoi(argv[1]));
    value_type L_end = value_type(atoi(argv[2]));
    value_type ensembleSize = value_type(atoi(argv[3]));

    moment_data(L_start, L_end, ensembleSize, 1.789);

}



/**
 *
 * @param ensemble_size
 * @param length
 * @param impure_sites
 */
void critical_point_data_sq_lattice_site_percolation(value_type length, value_type ensemble_size,
                                                     value_type impure_sites){
    SitePercolation_ps_v8 sp(length, impure_sites); // todo the program breaks down // possible solution remove redundant variables
    sp.periodicityFlag(false);

    cout << "Signature " << sp.signature << endl;
    clock_t t;
    vector<value_type> nos, nob;

    string common;
    common += to_string(length);
    common += '_';
    common += currentTime();;
    common += ".txt";

    string time_filename = "critical_point_data_L_" + common;

    ofstream fout_spanning_cluster(time_filename);      // contains the time information of the spanning cluster

    /**
     * First uncomment line is the signature
     * Second uncomment line is the parameter value
     * Third and all other uncomment lines are data
     */

    fout_spanning_cluster << "#<p_c>\t<id_sc>\t<b_t_sc>\t<sites_sc>\t<bonds_sc>" << endl;
    fout_spanning_cluster << "BEGIN_HEADER" << endl;
    fout_spanning_cluster << "ensemble_size\t"  << ensemble_size << endl;
    fout_spanning_cluster << "length\t"         << length << endl;
    fout_spanning_cluster << "data_line\t"      << 11 << endl;
    fout_spanning_cluster << "END_HEADER" << endl;
    fout_spanning_cluster << "#p_c = critical occupation probability" << endl;
    fout_spanning_cluster << "#id_sc = id of the spanning cluster" << endl;
    fout_spanning_cluster << "#b_t_sc = birth time of the spanning cluster" << endl;
    fout_spanning_cluster << "#sites_sc = number of sites of the spanning cluster" << endl;
    fout_spanning_cluster << "#bonds_sc = number of bonds of the spanning cluster" << endl;


//    double p_c{}, id_sc{}, b_t_sc{}, sites{}, bonds{};

    for(value_type i{} ; i != ensemble_size ; ++i){
        t = clock();
        sp.reset();

        sp.placeSiteUntilFirstSpanning_v2();

        fout_spanning_cluster << sp.occupationProbability() << '\t' << sp.firstSpanningClusterID_v2() << '\t'
                              << sp.birthTimeOfSpanningCluster() << '\t'
                              << sp.numberOfSitesInTheSpanningClusters_v2() << '\t'
                              << sp.numberOfBondsInTheSpanningClusters_v2() << endl;

        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    cout << "reached the end " << endl;
    fout_spanning_cluster.close();
    cout << "file is closed " << endl;
}



/**
 *
 * @param length
 * @param ensemble_size
 * @param rule : sum rule or product rule
 * @param impure_sites
 */
void critical_point_data_explosive(value_type length, value_type ensemble_size, int rule, value_type impure_sites){
    SitePercolation_ps_v8 sp(length, impure_sites); // todo breaks down the program // possible solution remove redundant variables
    sp.periodicityFlag(false);

    cout << "Signature " << sp.signature << endl;
    clock_t t;
    vector<value_type> nos, nob;

    string common;
    common += to_string(length);
    common += '_';
    if(rule == 0){
        common += "sum";
    }
    else if(rule == 1){
        common += "product";
    }
    else{
        cout << "Invalid rule : line " << __LINE__ << endl;
    }
    common += '_';
    common += currentTime();;
    common += ".txt";

    string time_filename = "critical_point_data_explosive_L_" + common;

    ofstream fout_spanning_cluster(time_filename);      // contains the time information of the spanning cluster

    /**
     * First uncomment line is the signature
     * Second uncomment line is the parameter value
     * Third and all other uncomment lines are data
     */

    fout_spanning_cluster << "#<p_c>\t<id_sc>\t<b_t_sc>\t<sites_sc>\t<bonds_sc>" << endl;
    fout_spanning_cluster << "#Data for explosive percolation" << endl;
    if(rule == 0){
        fout_spanning_cluster << "#using sum rule" << endl;
    }
    else if(rule == 1){
        fout_spanning_cluster << "#using product rule" << endl;
    }
    fout_spanning_cluster << "BEGIN_HEADER" << endl;
    fout_spanning_cluster << "ensemble_size\t"  << ensemble_size << endl;
    fout_spanning_cluster << "length\t"         << length << endl;
    fout_spanning_cluster << "data_line\t"      << 11 << endl;
    fout_spanning_cluster << "END_HEADER" << endl;
    fout_spanning_cluster << "#p_c = critical occupation probability" << endl;
    fout_spanning_cluster << "#id_sc = id of the spanning cluster" << endl;
    fout_spanning_cluster << "#b_t_sc = birth time of the spanning cluster" << endl;
    fout_spanning_cluster << "#sites_sc = number of sites of the spanning cluster" << endl;
    fout_spanning_cluster << "#bonds_sc = number of bonds of the spanning cluster" << endl;


//    double p_c{}, id_sc{}, b_t_sc{}, sites{}, bonds{};

    for(value_type i{} ; i != ensemble_size ; ++i){
        t = clock();
        sp.reset();

        sp.placeSite_explosive_UntilFirstSpanning(rule);

        fout_spanning_cluster << sp.occupationProbability() << '\t' << sp.firstSpanningClusterID_v2() << '\t'
                              << sp.birthTimeOfSpanningCluster() << '\t'
                              << sp.numberOfSitesInTheSpanningClusters_v2() << '\t'
                              << sp.numberOfBondsInTheSpanningClusters_v2() << endl;

        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    cout << "reached the end " << endl;
    fout_spanning_cluster.close();
    cout << "file is closed " << endl;
}



/**
 *
 * @param ensemble_size
 * @param length
 * @param impure_sites
 */
void critical_point_data_sq_l_sp_ballistic_deposition_L1(value_type length, value_type ensemble_size){
    SitePercolationBallisticDeposition sp(length, false); // todo the program breaks down // possible solution remove redundant variables

    cout << "Signature " << sp.signature << endl;
    clock_t t;
    vector<value_type> nos, nob;

    string common;
    common += to_string(length);
    common += '_';
    common += currentTime();;
    common += ".txt";

    string filename = "critical_point_data_sq_lattice_sp_ballistic_L1_L_" + common;

    ofstream fout_spanning_cluster(filename);      // contains the time information of the spanning cluster

    /**
     * First uncomment line is the signature
     * Second uncomment line is the parameter value
     * Third and all other uncomment lines are data
     */

    fout_spanning_cluster << "#<p_c>\t<id_sc>\t<b_t_sc>\t<sites_sc>\t<bonds_sc>" << endl;
    fout_spanning_cluster << "BEGIN_HEADER" << endl;
    fout_spanning_cluster << "ensemble_size\t"  << ensemble_size << endl;
    fout_spanning_cluster << "length\t"         << length << endl;
    fout_spanning_cluster << "data_line\t"      << 11 << endl;
    fout_spanning_cluster << "END_HEADER" << endl;
    fout_spanning_cluster << "#p_c = critical occupation probability" << endl;
    fout_spanning_cluster << "#id_sc = id of the spanning cluster" << endl;
    fout_spanning_cluster << "#b_t_sc = birth time of the spanning cluster" << endl;
    fout_spanning_cluster << "#sites_sc = number of sites of the spanning cluster" << endl;
    fout_spanning_cluster << "#bonds_sc = number of bonds of the spanning cluster" << endl;


    for(value_type i{} ; i != ensemble_size ; ++i){
        t = clock();
        sp.reset();

        value_type j = 0;
        bool flag{};
        while(true){
            flag = sp.occupy();
            if(flag){
                ++j;
            }

            if(j >= sp.length()){
                break;
            }

            if(sp.detectSpanning_v6(sp.lastPlacedSite())){
                cout << "Spanning at " << sp.occupationProbability() << endl;
                break;
            }
        }

        fout_spanning_cluster << sp.occupationProbability() << '\t' << sp.firstSpanningClusterID_v2() << '\t'
                              << sp.birthTimeOfSpanningCluster() << '\t'
                              << sp.numberOfSitesInTheSpanningClusters_v2() << '\t'
                              << sp.numberOfBondsInTheSpanningClusters_v2() << endl;

        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    cout << "reached the end " << endl;
    fout_spanning_cluster.close();
    cout << "file is closed " << endl;
}


