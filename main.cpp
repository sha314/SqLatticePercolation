//
//
//
//

#include <iostream>
#include <fstream>
#include <ctime>

#include "src/lattice/lattice.h"
#include "src/tests/test_percolation.h"
#include "src/percolation/percolation.h"
#include "src/latticeConfiguration.h"
#include "src/util/time_tracking.h"
#include "src/util/printer.h"
#include "src/fractal/fractal_tool.h"
#include "src/cmdLineArgs/percolation_fractal.h"
#include "src/cmdLineArgs/cmd_args.h"
#include "src/ziff_algorithm.h"
#include "src/cmdLineArgs/cmd_args.h"
#include "src/cmdLineArgs/cmd_args_json.h"


using namespace std;



/**
 *
 */
void fractalDimension();


void test_Fractal_2d();

/**
 * data for log(Mass) vs log(Length) curve
 */
void percolation_fractalDimension_by_site(
        value_type ensemble_size, value_type L_start, value_type l_end, value_type delta_L=20
)
{

    double M{};
    ofstream fout("fractal_dimension_mass_length_data.txt");
    fout << "#Measuring fractal dimension by all sites when spanned" << endl;
    fout << "#<Length>\t<Mass>" << endl;
    clock_t t, t_outer;
    for(value_type len{L_start}; len <= l_end ; len += 10){
        cout << "Length " << len << " : " << endl;
        SitePercolation_ps_v8 sp(len);
        M = 0;
        t_outer = clock();
        for(value_type i{} ; i != ensemble_size ; ++i){
            t = clock();
            sp.reset();
            sp.placeSiteUntilFirstSpanning_v2();
            M += sp.numberOfOccupiedSite();
            cout << "\t\tIteration " << i << " . birthTime " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
        }
        fout << len << '\t' << M / double(ensemble_size) << endl; // writing data to file
        cout << "Length " << len << " . Time taken "
             << getFormattedTime((clock() - t_outer) / double(CLOCKS_PER_SEC)) << endl;
    }
    fout.close();
}



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
        SitePercolation_ps_v8 sp(len);

        M_spanning_sites=0;
        M_sites=0;
        M_spanning_bonds=0;
        M_bonds=0;
        N_cluster=0;

        t_outer = clock();
        for(value_type i{} ; i != ensemble_size ; ++i){
            t = clock();
            sp.reset();
            sp.placeSiteUntilFirstSpanning_v2();
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
        SitePercolation_ps_v8 sp(len);
        M = 0;
        t_outer = clock();
        for(value_type i{} ; i != ensemble_size ; ++i){
            t = clock();
            sp.reset();
            sp.placeSiteUntilFirstSpanning_v2();
            M += sp.numberOfSitesInTheSpanningClusters();
            cout << "\t\tIteration " << i << " . birthTime " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
        }
        fout << len << '\t' << M / double(ensemble_size) << endl; // writing data to file
        cout << "Length " << len << " . Time taken "
             << getFormattedTime((clock() - t_outer) / double(CLOCKS_PER_SEC)) << endl;
    }
    fout.close();
}



void moment_test_2(value_type len, double d_f=1.789){
    //
    SitePercolation_ps_v8 sp(len);
    sp.periodicityFlag(false);
    ofstream fout("moment_conservation_" + to_string(sp.length()) + ".txt");
    fout << "# length is " << sp.length() << endl;
    fout << "# obtained data :" << endl;
    clock_t t;
    vector<value_type> vec;
    double M{}, M_avg{};
    size_t ensemble_size = 100;
    double sum_a{};
//    value_type _possible_neighbors{};     // number of occupied site in a lattice
    for(value_type i{} ; i != ensemble_size ; ++i){
        t = clock();
        sp.reset();
        sp.placeSiteUntilFirstSpanning_v2();
        vec = sp.number_of_site_in_clusters();
        M = 0;
        sum_a = 0;
        for(auto a: vec){
            fout << a << "\t";
            sum_a += a;
            M += pow(a, d_f);
        }
        if(sum_a != sp.numberOfOccupiedSite()){
            cerr << "Warning********-----------********" << endl;
        }
        fout << endl << "Total sites " <<  sum_a << endl;
        fout << "moment = " << M << endl << endl;
        M_avg += M / ensemble_size;
        cout << "\t\tIteration " << i << " . birthTime " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }
    fout << "Finally taking average : M_df = " << M_avg << endl;
}


/**
 * Counting number of blocks in the site when first spanning cluster appears
 * using box counting method
 * @param ensemble_size
 * @param len
 * @param delta
 */
void box_counting_test(value_type ensemble_size, value_type len, vector<value_type> delta){
    SitePercolation_ps_v8 sp(len);
//    sp.configure(LatticeSiteConfiguration::configuration_6x6_0());

//    cout << "number of unit element " << sp.box_counting(1) << endl;
//    cout << "number of unit element " << sp.box_counting(2) << endl;
//    cout << "number of unit element " << sp.box_counting(5) << endl;

    map<value_type, value_type> N_delta;
    clock_t t;
    double nos{};
    for(value_type i{} ; i != ensemble_size ; ++i){
        t = clock();
        cout << "Iteration " << i << " : ";
        sp.reset();
        sp.placeSiteUntilFirstSpanning_v2();
        nos += sp.numberOfOccupiedSite();
        for(auto d: delta){
            N_delta[d] +=  sp.box_counting(d);
            cout << "not defined " << __LINE__ << endl;
        }
        cout << nos <<"==" << N_delta[1] << endl;
        if(nos != double(N_delta[1])){
            cerr << "nos != double(N_delta[1])" << endl;
        }
        cout << " . Time elapsed " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }


    string filename("fractal_dim_");
    filename += to_string(len);
    filename += ".txt";

    ofstream fout(filename);
    fout << "#Length " << len << endl;
    fout << "#Ensemble size " << ensemble_size << endl;
    fout << "#Delta " << delta << endl;
    fout << "#<delta>\t<N(delta)>" << endl;
    for(auto d: delta){

        fout << d << '\t' << N_delta[d] / double(ensemble_size) << endl;

        if(d==1){
            if(double(N_delta[1]) != nos){
                cerr << "N_delta[1] != sp.numberOfOccupiedSite() : " << N_delta[1] << " != " << nos << endl;
            }
        }
    }
    fout.close();

}



/**
 * Counting number of blocks in the site when first spanning cluster appears
 * using box counting method
 * @param ensemble_size
 * @param len
 * @param delta
 */
void box_counting_function_check(value_type ensemble_size, value_type len){
    SitePercolation_ps_v8 sp(len);

    double N_delta{};
    clock_t t;
    for(value_type i{} ; i != ensemble_size ; ++i){
        t = clock();
        cout << "Iteration " << i << " : ";
        sp.reset();
        sp.placeSiteUntilFirstSpanning_v2();

        N_delta =  sp.box_counting(1); // define todo
        cout << "not defined " << __LINE__ << endl;
        cout << N_delta << " == " << sp.numberOfOccupiedSite();
        if(N_delta != sp.numberOfOccupiedSite()){
            cerr << "N_delta[1] != sp.numberOfOccupiedSite() : " <<N_delta << " != " << sp.numberOfOccupiedSite() << endl;
        }
        cout << " . Time elapsed " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }


}

void box_counting_function_check2(value_type len){
    SitePercolation_ps_v8 sp(len);

    sp.placeSiteUntilFirstSpanning_v2();
    sp.viewLatticeByID();

    vector<value_type> delta = get_delta(len);
    cout << "Box counting " << endl;
    cout << "Found delta " << delta << endl;
    array<value_type , 2> tmp;
    for(auto d: delta){
        tmp = sp.box_counting_v2(d);
        cout << d << "->" << tmp[0] << "," << tmp[1] << endl;
    }

}


void set_comparator(){
    set<Site, ClusterIDComparator> sets;
    Site a({0, 0}, 5);
    Site b({0, 2}, 5);
    Site c({0, 1}, 5);
    a.groupID(0);
    b.groupID(1);
    c.groupID(2);

    sets.insert(a);
    sets.insert(b);
    sets.insert(c);
    for(const Site s: sets){
        cout << s << endl;
    }
}

void test(){
    SitePercolation_ps_v8 sp(6);
    sp.periodicityFlag(true);
    sp.logging_flag(true);

    size_t sz = sp.length()*sp.length();
    sz = 22;
    for(size_t i{}; i <= sz ; ++i){
        sp.placeSite_weighted_v8();
//        sp.placeSite_v7();
        cout << i << " th site " << sp.lastPlacedSite() << endl;
        sp.viewLatticeByID();
        sp.viewClusterExtended();
        sp.viewCluster_id_index();
    }

    sp.viewLatticeByID();
    sp.viewClusterExtended();
    sp.viewCluster_id_index();
    sp.placeSite_weighted_v8();
    sp.viewLatticeByID();
    sp.viewClusterExtended();
    sp.viewCluster_id_index();
    sp.logs();
}


void test_spanning(){
    SitePercolation_ps_v8 sp(10);
    sp.periodicityFlag(false);
//    while (sp.placeSiteForSpanning()){
    while (sp.placeSite_explosive_sum_rule_ForSpanning()){

        if(sp.isSpanned()){
            sp.viewLatticeByID();
//            sp.spanningIndices();
            break;
        }
//        sp.viewLatticeByID();
    }

    sp.spanningIndices();
    sp.viewLatticeByID();

}


void test_bond_percolation(){
    BondPercolation_pb_v0 bp(5, false);
    bp.viewClusterExtended();
    bp.viewLatticeByID();
    cout << bp.lastPlacedBond() << endl;

    value_type i{};
    while(bp.occupy()) {
        cout << i << " th bond " << bp.lastPlacedBond() << endl;
        bp.viewLatticeByID();
        bp.viewClusterExtended();
        ++i;
//        if( i == 7) {
//            break;
//        }
    }
//    bp.viewClusterExtended();
}



void test_site_percolation(){
//    srand(1526028791);    // seeding
    SitePercolation_ps_v8 sp(100, false);
//    sp.periodicityFlag(true);
    value_type i = 0;
    while(sp.occupy()) {
        ++i;
//        sp.viewClusterExtended();
//        sp.viewSiteByID();
//        sp.viewSiteByRelativeIndex();
        if(sp.detectSpanning_v6(sp.lastPlacedSite())){
//            cout << "last site " << sp.lastPlacedSite() << endl;
//            cout << "pc = " << sp.occupationProbability() << endl;
            sp.numberOfBondsInTheSpanningClusters();
        }
//        if ( i == 25){
//            break;
//        }
    }
    cout << "spanning cluster ids " << endl;
    sp.numberOfBondsInTheSpanningClusters();
    cout << "Function ended" << endl;
}


void test_ballistic_deposition(size_t seed){
    value_type length = 10;
    srand(seed); ///// seeding
    SitePercolationBallisticDeposition_L2 sp(length, false);

    cout << "Signature : " << sp.signature << endl;

    size_t length_squared = length * length;
    std::vector<double> nos(length_squared), nob1(length_squared), nob2(length_squared), entrpy(length_squared);
    size_t j{};
    bool wrapping_occured{false};

    sp.reset();
    j = 0;
    wrapping_occured = false;
    bool successful = false;
    double pc{};
    while (true) {
        successful = sp.occupy();
        if (successful) {
            cout << "last site " << sp.lastPlacedSite() << endl;
            sp.viewSiteByID();
            if(sp.periodicity()) {
//                sp.viewSiteByRelativeIndex();
//                if(sp.detectWrapping_v1(sp.lastPlacedSite())){
//                    sp.viewSiteByRelativeIndex();
//                    cout << "***********************Wrapping ************************" << endl;
//                    sp.wrappingIndices();
//                    return;
//                }
                if (!wrapping_occured && sp.detectWrapping_v1(sp.lastPlacedSite())) {
                    wrapping_occured = true;
                    pc = sp.occupationProbability();
//                    sp.writeVisualLatticeData("lattice-visual-data-"+to_string(length)+".txt", false);
//                    sp.viewSiteByID();
//                    sp.viewClusterExtended();

                }
                if (wrapping_occured) {
                    nob2[j] += sp.numberOfBondsInTheWrappingClusters();
                }
            }
            else{
                if (!wrapping_occured && sp.detectSpanning_v6(sp.lastPlacedSite())) {
                    wrapping_occured = true;
                    pc = sp.occupationProbability();

                }
                if (wrapping_occured) {
                    nob2[j] += sp.numberOfBondsInTheSpanningClusters_v2();
                }
            }
            nos[j] += sp.numberOfOccupiedSite();
            nob1[j] += sp.numberOfBondsInTheLargestCluster_v2();
//            entrpy[j] += sp.entropy();  // old method and takes long time
            entrpy[j] += sp.entropy_v3(); // faster method
            ++j;
        }
        if (j >= length_squared) { // length_squared is the number of site
            break;
        }
    }
//    return;
//    cout << "p\tH\tP\tP" << endl;
//    for (value_type i{}; i != length_squared; ++i){
//        cout << nos[i] << "\t" << entrpy[i] << "\t" << nob1[i] << "\t" << nob2[i] << endl;
//    }
//    cout << "pc = " << pc << endl;

}


/**************************************
 *  The main function
 *
 ***************************************/
int main(int argc, char** argv) {
    cout << currentTime() << endl;

    cout << "Compiled on " << __DATE__ << "\t at " << __TIME__ << endl;
    std::cout << "Percolation in a Square Lattice" << std::endl;
    clock_t t = clock();
    auto seed = time(NULL);
//    srand(seed);    // seeding

//    cmd_moment_conservation(argc, argv);

//    cmd_args_site_percolation(argc, argv); // 2017.12.22
//    cmd_args_site_percolation_ballistic_deposition(argc, argv); // 2018.04.02

    cmd_args_json_g(argc, argv); // 2018.06.11

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

//    cout << "Seed " << seed << endl;
    std::cout << "Program Ended : Time elapsed "
              << getFormattedTime((clock() - t) / double(CLOCKS_PER_SEC)) << std::endl;
    return 0;
}


