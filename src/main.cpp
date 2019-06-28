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
#include "tests/test_v11.h"


using namespace std;



/****
 *  All the function that is run in main
 * @param argc
 * @param argv
 */
void run_in_main(int argc, char** argv){
    cout << "length = atoi(argv[1])" << endl;
    cout << "ensemble_size = atoi(argv[2])" << endl;
    value_type length = atoi(argv[1]);
    value_type ensemble_size = atoi(argv[2]);

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

//    test_bond_percolation(length, ensemble_size);
//    test_site_percolation(length, ensemble_size);



//    cluster_size_distribution_site(length, ensemble_size); // 2019.06.17
//    cluster_size_distribution_bond(length, ensemble_size); // 2019.06.22
//    measure_entropy_by_site(argc, argv);
//    measure_clusters(argc, argv);
//    weighted_relabeling_test(argc, argv);

//    bond_percolation(argc, argv);

//    simulate_site_percolation(length, ensemble_size);
//    simulate_site_percolation_T<SitePercolationBallisticDeposition_L1_v2>(length, ensemble_size); // 2018.11.03
//    simulate_site_percolation_T<SitePercolationBallisticDeposition_L2_v2>(length, ensemble_size); // 2018.11.03
//    simulate_site_percolation_detailed(length, ensemble_size); // 2018.12.18
//    simulate_bond_percolation_v2(length, ensemble_size);
//    visual(length);
//    simulate_site_percolation();
//    entropyJumps(argc, argv);
//    site_percolation_new_version(argc, argv);

    // date 2019.06.28
    test_sq_lattice_v2(length);
    test_percolation(length);

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

//    time_t seed = time(NULL);
//    srand(seed);    // seeding

    run_in_main(argc, argv);

    auto t_end= std::chrono::system_clock::now();
    std::chrono::duration<double> drtion = t_end - t_start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(t_end);
    cout << "Program finished at " << std::ctime(&end_time) << endl;
    std::cout << "Time elapsed "   << getFormattedTime(drtion.count()) << std::endl;
    return 0;
}


