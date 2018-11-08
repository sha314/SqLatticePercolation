//
//
//
//

#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <thread>
#include <mutex>

#include "lattice/lattice.h"
#include "percolation/percolation.h"
#include "util/time_tracking.h"
#include "util/printer.h"
#include "fractal/fractal_tool.h"
#include "cmdLineArgs/cmd_args.h"
#include "cmdLineArgs/cmd_args.h"
#include "cmdLineArgs/cmd_args_json.h"
#include "tests/test_percolation.h"


using namespace std;





/**
 * data for log(Mass) vs log(Length) curve
 * for spanning sites only
 */
void percolation_to_get_fractal(
        value_type ensemble_size, value_type L_start, value_type l_end, value_type delta_L=20
)
{
    double M_spanning_sites{}, M_sites{}, M_spanning_bonds{}, M_bonds{}, N_cluster{};
    ofstream fout("spanning_data_diff_length.txt");
    fout << "#All data when spanning occurs for different length" << endl;
    fout << "#<Length>\t<Number of cluster>\t<Number of spanning sites>\t<Total number of sites> \t<Number of spanning bonds>" << endl;
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
        value_type ensemble_size, value_type L_start, value_type l_end, value_type delta_L=20
)
{
    double M{};
    ofstream fout("mass_length_data.txt");
    fout << "#Measuring fractal dimension by only spanning sites" << endl;
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
    fout << "#cluster size is measured by number of sites in it" << endl;

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

/****
 *  All the function that is run in main
 * @param argc
 * @param argv
 */
void run_in_main(int argc, char** argv){

//    cmd_moment_conservation(argc, argv);

//    cmd_args_site_percolation(argc, argv); // 2017.12.22
//    cmd_args_site_percolation_ballistic_deposition(argc, argv); // 2018.04.02

//    cmd_args_json_g(argc, argv); // 2018.06.11

//    cout << "Finished cmd_args_site_percolation" << endl;

//    fractalArea_by_bonds(10, 5);
//    percolation_fractalDimension_by_spanning_site(500, 30, 150, 20);
//    percolation_to_get_fractal(500, 30, 150, 20);
//    percolation_fractalDimension_by_site(500, 30, 90, 10);

//    moment_data(50, 100, 100, 1.789);

//    moment_test_2(30);
//    moment_test_2(50);
//    moment_test_2(70);

    //    box_counting_function_check(10, 100);

//    box_counting_test(100, 50, {1, 2, 5, 10});

//    test_Fractal_2d();

//    test_bond_percolation();
//    test_ballistic_deposition(1);
//    test_site_percolation();


//    cluster_size_threaded(argc, argv);
//    cluster_size(argc, argv);
//    measure_entropy_by_site(argc, argv);
//    measure_clusters(argc, argv);
//    weighted_relabeling_test(argc, argv);

//    bond_percolation(argc, argv);

//    simulate_site_percolation(argc, argv);
    simulate_site_percolation_T<SitePercolationBallisticDeposition_L1_v2>(argc, argv);
    simulate_site_percolation_T<SitePercolationBallisticDeposition_L2_v2>(argc, argv);
//    simulate_site_percolation();
//    entropyJumps(argc, argv);
//    site_percolation_new_version(argc, argv);

}




/**************************************
 *  The main function
 *
 ***************************************/
int main(int argc, char** argv) {
    cout << currentTime() << endl;

    cout << "Compiled on " << __DATE__ << "\t at " << __TIME__ << endl;
    std::cout << "Percolation in a Square Lattice" << std::endl;
    auto t_start = std::chrono::system_clock::now();

    time_t seed = time(NULL);
    srand(seed);    // seeding

    run_in_main(argc, argv);

    auto t_end= std::chrono::system_clock::now();
    std::chrono::duration<double> drtion = t_end - t_start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(t_end);
    cout << "Program finished at " << std::ctime(&end_time) << endl;
    std::cout << "Time elapsed "   << getFormattedTime(drtion.count()) << std::endl;
    return 0;
}


