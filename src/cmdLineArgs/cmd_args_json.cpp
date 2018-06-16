//
// Created by Shahnoor on 6/7/2018.
//

#include <fstream>
#include "cmd_args_json.h"
#include "../percolation/percolation.h"
#include "../util/time_tracking.h"


using namespace std;


/**
 * Calculates the following with periodic bounding condition
 *  1. critical occupation probability, p_c
 *  2. Order Parameter, P(p,L)
 *  3. Entropy, H(p,L)
 * @param length
 * @param ensembleSize
 */
void percolation_sq_lattice_periodic(value_type length, value_type ensembleSize) {
    SitePercolation_ps_v8 sp(length, true);
    cout << "Signature " << sp.signature << endl;

    clock_t t;
    size_t length_squared = length*length;
    std::vector<double> nos(length_squared), nob1(length_squared), nob2(length_squared), entrpy(length_squared);
    size_t j{};
    bool wrapping_occured {false};
    std::vector<double> pcs(ensembleSize);
    for(value_type i{} ; i != ensembleSize ; ++i){

        t = clock();
        sp.reset();
        j = 0;
        wrapping_occured = false;
        bool successful = false;
        while (true){
            successful = sp.occupy();
            if(successful) {
                if (!wrapping_occured && sp.detectWrapping_v1(sp.lastPlacedSite())) {
                    wrapping_occured = true;
                    pcs[i] = sp.occupationProbability();
                }
                if (wrapping_occured) {
                    nob2[j] += sp.numberOfBondsInTheWrappingClusters();
                }
                nos[j] += sp.numberOfOccupiedSite();
                nob1[j] += sp.numberOfBondsInTheLargestCluster_v2();
//            entrpy[j] += sp.entropy();  // old method and takes long time
                entrpy[j] += sp.entropy_v3(); // faster method
                ++j;
            }
            if(j >= length_squared){ // length_squared is the number of site
                break;
            }
        }
        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    // Taking Average
    for(size_t i{}; i!= length_squared ; ++i){
        nos[i] /= double(ensembleSize);
        entrpy[i] /= double(ensembleSize);
        nob1[i] /= double(ensembleSize);
        nob2[i] /= double(ensembleSize);
    }

    // calculating and writing to a file
    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length
                << ", \"ensemble_size\":" << ensembleSize
                << ", \"signature\":" << sp.getSignature()
                << "}" << endl;


    string tm = currentTime();
    string filename = sp.getSignature() + "_" + to_string(length) + "_" + tm;
    filename += ".txt";
    write_percolation_data_json(length, nob1, nob2, entrpy, filename, header_info.str());

    filename = sp.getSignature() + "_critical_" + to_string(length) + "_" + tm;
    filename += ".txt";
    write_critical_data_json(pcs, filename, header_info.str());
}


/**
 * Calculates the following with non-periodic bounding condition
 *  1. critical occupation probability, p_c
 *  2. Order Parameter, P(p,L)
 *  3. Entropy, H(p,L)
 * @param length
 * @param ensembleSize
 */
void percolation_sq_lattice_non_periodic(value_type length, value_type ensembleSize) {
    SitePercolation_ps_v8 sp(length, false);
    cout << "Signature " << sp.signature << endl;

    clock_t t;
    size_t length_squared = length*length;
    std::vector<double> nos(length_squared), nob1(length_squared), nob2(length_squared), entrpy(length_squared);
    size_t j{};
    bool spanning_occured {false};
    std::vector<double> pcs(ensembleSize);
    bool successful{false};
    for(value_type i{} ; i != ensembleSize ; ++i){

        t = clock();
        sp.reset();
        j = 0;
        spanning_occured = false;
        successful = false;
        while (true){
            successful = sp.occupy();
            if(successful) {
                if (!spanning_occured && sp.detectSpanning_v6(sp.lastPlacedSite())) {
                    spanning_occured = true;
                    pcs[i] = sp.occupationProbability();
                }
                if (spanning_occured) {
                    nob2[j] += sp.numberOfBondsInTheSpanningClusters();
                }
                nos[j] += sp.numberOfOccupiedSite();
                nob1[j] += sp.numberOfBondsInTheLargestCluster_v2();
//            entrpy[j] += sp.entropy();  // old method and takes long time
                entrpy[j] += sp.entropy_v3(); // faster method
                ++j;
            }
            if(j >= length_squared){ // length_squared is the number of site
                break;
            }
        }
        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    // Taking Average
    for(size_t i{}; i!= length_squared ; ++i){
        nos[i] /= double(ensembleSize);
        entrpy[i] /= double(ensembleSize);
        nob1[i] /= double(ensembleSize);
        nob2[i] /= double(ensembleSize);
    }

    // calculating and writing to a file
    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length
                << ", \"ensemble_size\":" << ensembleSize
                << ", \"signature\":" << sp.getSignature()
                << "}" << endl;


    string tm = currentTime();
    string filename = sp.getSignature() + "_" + to_string(length) + "_" + tm;
    filename += ".txt";
    write_percolation_data_json(length, nob1, nob2, entrpy, filename, header_info.str());

    filename = sp.getSignature() + "_critical_" + to_string(length) + "_" + tm;
    filename += ".txt";
    write_critical_data_json(pcs, filename, header_info.str());
}

////////////// Generic
/**
 * Calculates the following with periodic bounding condition
 *  1. critical occupation probability, p_c
 *  2. Order Parameter, P(p,L)
 *  3. Entropy, H(p,L)
 * @param length
 * @param ensembleSize
 */
template <class PERCOLATION>
void percolation_sq_lattice_periodic_g(value_type length, value_type ensembleSize) {
    PERCOLATION sp(length, true);
    cout << "Signature " << sp.signature << endl;

    clock_t t;
    size_t length_squared = length*length;
    std::vector<double> nos(length_squared), nob1(length_squared), nob2(length_squared), entrpy(length_squared);
    size_t j{};
    bool wrapping_occured {false};
    std::vector<double> pcs(ensembleSize);
    for(value_type i{} ; i != ensembleSize ; ++i){

        t = clock();
        sp.reset();
        j = 0;
        wrapping_occured = false;
        bool successful = false;
        while (true){
            successful = sp.occupy();
            if(successful) {
                if (!wrapping_occured && sp.detectWrapping_v1(sp.lastPlacedSite())) {
                    wrapping_occured = true;
                    pcs[i] = sp.occupationProbability();
                }
                if (wrapping_occured) {
                    nob2[j] += sp.numberOfBondsInTheWrappingClusters();
                }
                nos[j] += sp.numberOfOccupiedSite();
                nob1[j] += sp.numberOfBondsInTheLargestCluster_v2();
//            entrpy[j] += sp.entropy();  // old method and takes long time
                entrpy[j] += sp.entropy_v3(); // faster method
                ++j;
            }
            if(j >= length_squared){ // length_squared is the number of site
                break;
            }
        }
        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    // Taking Average
    for(size_t i{}; i!= length_squared ; ++i){
        nos[i] /= double(ensembleSize);
        entrpy[i] /= double(ensembleSize);
        nob1[i] /= double(ensembleSize);
        nob2[i] /= double(ensembleSize);
    }

    // calculating and writing to a file
    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length
                << ", \"ensemble_size\":" << ensembleSize
                << ", \"signature\":" << sp.getSignature()
                << "}" << endl;


    string tm = currentTime();
    string filename = sp.getSignature() + "_" + to_string(length) + "_" + tm;
    filename += ".txt";
    write_percolation_data_json(length, nob1, nob2, entrpy, filename, header_info.str());

    filename = sp.getSignature() + "_critical_" + to_string(length) + "_" + tm;
    filename += ".txt";
    write_critical_data_json(pcs, filename, header_info.str());
}


/**
 * Calculates the following with non-periodic bounding condition
 *  1. critical occupation probability, p_c
 *  2. Order Parameter, P(p,L)
 *  3. Entropy, H(p,L)
 * @param length
 * @param ensembleSize
 */
template <class PERCOLATION>
void percolation_sq_lattice_non_periodic_g(value_type length, value_type ensembleSize) {
    PERCOLATION sp(length, false);
    cout << "Signature " << sp.signature << endl;

    clock_t t;
    size_t length_squared = length*length;
    std::vector<double> nos(length_squared), nob1(length_squared), nob2(length_squared), entrpy(length_squared);
    size_t j{};
    bool spanning_occured {false};
    std::vector<double> pcs(ensembleSize);
    for(value_type i{} ; i != ensembleSize ; ++i){

        t = clock();
        sp.reset();
        j = 0;
        spanning_occured = false;
        while (sp.occupy()){
            if(!spanning_occured && sp.detectSpanning_v6(sp.lastPlacedSite())){
                spanning_occured = true;
                pcs[i] = sp.occupationProbability();

//                cout << sp.occupationProbability();
            }
            if(spanning_occured){
                nob2[j] += sp.numberOfBondsInTheSpanningClusters();
            }
            nos[j] += sp.numberOfOccupiedSite();
            nob1[j] += sp.numberOfBondsInTheLargestCluster_v2();
//            entrpy[j] += sp.entropy();  // old method and takes long time
            entrpy[j] += sp.entropy_v3(); // faster method
            ++j;
        }
        cout << "\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    // Taking Average
    for(size_t i{}; i!= length_squared ; ++i){
        nos[i] /= double(ensembleSize);
        entrpy[i] /= double(ensembleSize);
        nob1[i] /= double(ensembleSize);
        nob2[i] /= double(ensembleSize);
    }

    // calculating and writing to a file
    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length
                << ", \"ensemble_size\":" << ensembleSize
                << ", \"signature\":" << sp.getSignature()
                << "}" << endl;


    string tm = currentTime();
    string filename = sp.getSignature() + "_" + to_string(length) + "_" + tm;
    filename += ".txt";
    write_percolation_data_json(length, nob1, nob2, entrpy, filename, header_info.str());

    filename = sp.getSignature() + "_critical_" + to_string(length) + "_" + tm;
    filename += ".txt";
    write_critical_data_json(pcs, filename, header_info.str());
}



/**
 *
 * @param length
 * @param ensembleSize
 * @param pcs
 * @param filename
 */
void write_critical_data_json(vector<double> pcs, string filename, const string &header_info) {
    ofstream fout(filename);
    // JSON formated header
    fout << header_info;
    fout << "#critical occupation probability or pc" << endl;
    fout << "#<pc>\n";

    for(double pc : pcs){
        fout << pc << endl; // writing value of pc to the file
    }
}


/**
 *
 * @param length
 * @param ensembleSize
 * @param total_number_of_bonds
 * @param nob1
 * @param nob2
 * @param entrpy
 * @param filename
 */
void write_percolation_data_json(
        value_type length,
        const vector<double> &nob1,
        const vector<double> &nob2,
        const vector<double> &entrpy,
        const string &filename,
        const string &header_info
)
{
    size_t length_squared = length*length;
    double total_number_of_bonds = 2*length_squared;
    ofstream fout(filename);
    // JSON formated header
    fout << header_info;

    string column_names = "#<p>\t<H(p,L)>\t<P1(p,L)>\t<P2(p,L)>\n";
    fout << column_names;

    string details = "#p = occupation probability\n";
    details += "#T = temperature = 1-p \n";
    details += "#H(p,L) = Entropy = sum( - u_i * log(u_i))\n";
    details += "#P1(p,L) = Order parameter = (number of bonds in largest cluster) / (total number of bonds)\n";
    details += "#P2(p,L) = Order parameter = (number of bonds in spanning or wrapping cluster) / (total number of bonds)\n";
    details += "#C(p,L) = Specific heat = -T dH/dT\n";
    details += "#X(p,L) = Susceptibility = dP/dp \n";
    details += "#u_i = (number of bonds in the i-th cluster) / (total number of bonds)\n";

    fout << details << endl;


    double p, H, P1, P2;
    for(size_t i{}; i!= length_squared; ++i){
        p = (i+1)/ double(length_squared);
        H = entrpy[i];                      // Entropy
        P1 = nob1[i] / total_number_of_bonds;    // Order parameter
        P2 = nob2[i] / total_number_of_bonds;    // Order parameter

        fout << p << '\t' << H << '\t' << P1 << '\t' << P2 << endl;
    }
    fout.close();
}


void print_help_json(){
    std::string hlp = R"***(
Percolation on a Square Lattice

Available Argument

-h      Help
-l      Length of the lattice. Default value is 5
-n      Ensemble size. Default value is 1
-t      Percolation Type. Default value is 0
            0   Regular Site percolation or (L0)

            1   Ballistic Deposition L1 (occupy up to first nearest
                    neighbor if selected site is already occupied. L0 is included here.)

            2   Ballistic Deposition L2 (occupy up to second nearest
                    neighbor along the first nearest neighbor only if
                    the selected site and first nearest neighbor is occupied. L0 and L1 is included here)

-p      place sites for. Default value is 0
            0   ->  Just place all the sites and printout the result in the commandline only once to see
                    what's going on there and how much time it requires.

            1   ->  MicroCanonical Ensemble data of Site percolation (Periodic)
                        Occupation Probability, p
                        Entropy, H = sum u_i log(u_i)
                        Order Parameter, P_1 =(bonds in the largest cluster) / (total bonds)
                        Order Parameter, P_2 =(bonds in the wrapping cluster) / (total bonds)
                    And Critical Occupation Probabilities (pc) data

            2   ->  MicroCanonical Ensemble data of Site percolation (Non-Periodic)
                        Occupation Probability, p
                        Entropy, H = sum u_i log(u_i)
                        Order Parameter, P_1 =(bonds in the largest cluster) / (total bonds)
                        Order Parameter, P_2 =(bonds in the spanning cluster) / (total bonds)
                    And Critical Occupation Probabilities (pc) data

    )***";
    std::cout << hlp << std::endl;
}


/**
 *
 * @param argc
 * @param argv
 */
void cmd_args_json(int argc, char **argv) {
    if(argc == 1){
        print_help_json();
        exit(0);
    }
    value_type ensembleSize{1}, length{5};
    value_type impure_site{0};
    int place_sites_for{-1};
    int percolation_type{0};


    char * filename=NULL;

    for(int i{1} ; i < argc;){
        std::cout << i << "\t" << argv[i] << std::endl;

        switch (argv[i][1]){

            case 'n':
                // ensemble size
                ensembleSize = value_type(atoi(argv[i+1]));
                i += 2;
                break;

            case 'l':
                // length
                length = value_type(atoi(argv[i+1]));
                i += 2;
                break;

            case 't':
                // range of length
                percolation_type = atoi(argv[i+1]);
                i += 2;
                break;

            case 'p':
                place_sites_for = value_type(atoi(argv[i+1]));
                i += 2;
                break;

            case 'h':
                print_help_json();
                exit(0);

            default:
                print_help_json();
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


    if(place_sites_for == 0) {
        clock_t t0 = clock();
        SitePercolation_ps_v8 sp(length, true);
        cout << "Signature " << sp.signature << endl;
        while (sp.occupy());
        cout << "Time elapsed : " << (clock() - t0) / double(CLOCKS_PER_SEC) << " sec" << endl;
        if (length <= 50) {
            sp.viewLatticeExtended();
            sp.viewClusterExtended();
            sp.viewCluster_id_index();
        }
    }else if(place_sites_for == 1) {
        // place sites until first spanning and generate all data possible
        // can generate several gigabytes of data
        percolation_sq_lattice_periodic(length, ensembleSize);
    }
    else if(place_sites_for == 2) {
        // place sites until first spanning and generate all data possible
        // can generate several gigabytes of data
        percolation_sq_lattice_non_periodic(length, ensembleSize);
    }
}


/**
 *
 * @param argc
 * @param argv
 */
void cmd_args_json_g(int argc, char **argv) {
    if(argc == 1){
        print_help_json();
        exit(0);
    }
    value_type ensembleSize{1}, length{5};
    int place_sites_for{-1};
    int percolation_type{0};


    char * filename=NULL;

    for(int i{1} ; i < argc;){
        std::cout << i << "\t" << argv[i] << std::endl;

        switch (argv[i][1]){

            case 'n':
                // ensemble size
                ensembleSize = value_type(atoi(argv[i+1]));
                i += 2;
                break;

            case 'l':
                // length
                length = value_type(atoi(argv[i+1]));
                i += 2;
                break;

            case 't':
                // range of length
                percolation_type = atoi(argv[i+1]);
                i += 2;
                break;

            case 'p':
                place_sites_for = value_type(atoi(argv[i+1]));
                i += 2;
                break;

            case 'h':
                print_help_json();
                exit(0);

            default:
                print_help_json();
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

    if(percolation_type == 0){
        // regular site percolation
        simulate <SitePercolation_ps_v8> (place_sites_for, length, ensembleSize);
    }
    else if(percolation_type == 1){
        simulate <SitePercolationBallisticDeposition_L1> (place_sites_for, length, ensembleSize);
    }
    else if(percolation_type == 2){
        simulate <SitePercolationBallisticDeposition_L2> (place_sites_for, length, ensembleSize);
    }
}


/**
 * Simulate based on Template and argument
 * @tparam PERCOLATION
 * @param p
 * @param length
 * @param ensembleSize
 */
template <class PERCOLATION>
void simulate(int p, value_type length, value_type ensembleSize){

    if(p == 0) {
        clock_t t0 = clock();
        PERCOLATION sp(length, true);
        cout << "Signature " << sp.signature << endl;
        while (sp.occupy());
        cout << "Time elapsed : " << (clock() - t0) / double(CLOCKS_PER_SEC) << " sec" << endl;
        if (length <= 50) {
            sp.viewLatticeExtended();
            sp.viewClusterExtended();
            sp.viewCluster_id_index();
        }
    }else if(p == 1) {
        // place sites until first spanning and generate all data possible
        // can generate several gigabytes of data
        percolation_sq_lattice_periodic_g<PERCOLATION>(length, ensembleSize);
    }
    else if(p == 2) {
        // place sites until first spanning and generate all data possible
        // can generate several gigabytes of data
        percolation_sq_lattice_non_periodic_g<PERCOLATION>(length, ensembleSize); // todo stackdump occurs
    }
}