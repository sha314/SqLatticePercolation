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
#include "tests/test_percolation.h"
#include "percolation/percolation.h"
#include "latticeConfiguration.h"
#include "util/time_tracking.h"
#include "util/printer.h"
#include "fractal/fractal_tool.h"
#include "cmdLineArgs/percolation_fractal.h"
#include "cmdLineArgs/cmd_args.h"
#include "ziff_algorithm.h"
#include "cmdLineArgs/cmd_args.h"
#include "cmdLineArgs/cmd_args_json.h"


using namespace std;



/**
 *
 */
void fractalDimension();


void test_Fractal_2d();

void cluster_size_threaded(int argc, char **argv);

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
    a.set_groupID(0);
    b.set_groupID(1);
    c.set_groupID(2);

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
//                if(sp.detectWrapping(sp.lastPlacedSite())){
//                    sp.viewSiteByRelativeIndex();
//                    cout << "***********************Wrapping ************************" << endl;
//                    sp.wrappingIndices();
//                    return;
//                }
                if (!wrapping_occured && sp.detectWrapping()) {
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
            entrpy[j] += sp.entropy(); // faster method
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


void single_thread(
        value_type length,
        value_type ensemble_size,
        value_type thread_id
){
    cout << "Not good : line " <<__LINE__ << endl;
    SitePercolation_ps_v8 sp(length, true);

    sp.simulate_periodic_critical(ensemble_size);

    vector<vector<double>> critical(3);
//    critical[0] = sp.get_pcs();
//    critical[1] = sp.get_spanning_cluster_size_bonds();
//    critical[2] = sp.get_spanning_cluster_size_sites();

    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length
                << ", \"ensemble_size\":" << ensemble_size
                << ", \"signature\":\"" << sp.getSignature() << "\""
                << "}" << endl;

    string tm = currentTime();

    string filename = sp.getSignature() + "_critical_" + to_string(length) + "_thread_" + to_string(thread_id) + '_' + tm;
    filename += ".txt";
    write_critical_data_json_v2(critical, filename, header_info.str());

}


void cluster_size_threaded(int argc, char **argv) {
    value_type length = atoi(argv[1]);
    value_type ensemble_size = atoi(argv[2]);


    cout << length << ' ' << ensemble_size << endl;

//    vector<vector<double>> _critical_data(3);
//    single_thread(length, ensemble_size, _critical_data);

    value_type number_of_threads = thread::hardware_concurrency();
//    value_type number_of_threads = 1;
//
    cout << "number of threads " << number_of_threads << endl;
//
    vector<thread> threads(number_of_threads);
//    vector<value_type> iterations(number_of_threads);
//    vector<vector<vector<double>>> _critical_data(number_of_threads);
//

////    _critical_data[0].resize(ensemble_size);
////    _critical_data[1].resize(ensemble_size);
////    _critical_data[2].resize(ensemble_size);
    value_type loops_per_thread = ensemble_size / number_of_threads;

    for(value_type i{}; i != number_of_threads; ++ i){

        threads[i] = thread(
                single_thread,
                length,
                loops_per_thread,
                i
        );
    }



    for(value_type i{}; i != number_of_threads; ++ i){
        if(threads[i].joinable()){
            cout << " line : " << __LINE__ << endl;
            threads[i].join();
        }
    }

    cout << " line : " << __LINE__ << endl;

    SitePercolation_ps_v8 sp(10, true);

}



void cluster_size(int argc, char **argv) {
    value_type length = atoi(argv[1]);
    value_type ensemble_size = atoi(argv[2]);

    cout << "Not good " << __LINE__ << endl;
    cout << "length " << length << " ensemble_size " << ensemble_size << endl;

    vector<vector<double>> _critical_data(3);

    SitePercolation_ps_v8 sp(length, true);

    sp.simulate_periodic_critical(ensemble_size);

//    _critical_data[0] = sp.get_pcs();
//    _critical_data[1] = sp.get_spanning_cluster_size_bonds();
//    _critical_data[2] = sp.get_spanning_cluster_size_sites();

    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length
                << ", \"ensemble_size\":" << ensemble_size
                << ", \"signature\":\"" << sp.getSignature() << "\""
                << "}" << endl;

    string tm = currentTime();

    string filename = sp.getSignature() + "_critical_" + to_string(length) + '_' + tm;
    filename += ".txt";
    write_critical_data_json_v2(_critical_data, filename, header_info.str());

}

void measure_entropy_by_site(int argc, char** argv){
    value_type length = atoi(argv[1]);
    value_type ensemble_size = atoi(argv[2]);


    cout << "length " << length << " ensemble_size " << ensemble_size << endl;

    value_type length_squared = length*length;
    vector<double> entropy(length_squared);
//    vector<double> entropy2(length_squared);

    SitePercolation_ps_v8 sp(length, true);
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

/**
 * For rank-size distribution
 * @param argc
 * @param argv
 */
void measure_clusters(int argc, char** argv){
    value_type length = atoi(argv[1]);
    value_type ensemble_size = atoi(argv[2]);


    cout << "length " << length << " ensemble_size " << ensemble_size << endl;

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    SitePercolation_ps_v8 sp(length, true);

    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length
                << ", \"ensemble_size\":" << ensemble_size
                << ", \"signature\":\"" << sp.getSignature() << "\""
                << "}" << endl;

    string tm = currentTime();

    string filename_s = sp.getSignature() + "_cluster_by_site_" + to_string(length) + '_' + tm;
    string filename_b = sp.getSignature() + "_cluster_by_bond_" + to_string(length) + '_' + tm;
    filename_s += ".txt";
    filename_b += ".txt";

    ofstream fout_s(filename_s);
    // JSON formated header
    fout_s << header_info.str();
    fout_s << "#each line is an independent realization" << endl;
    fout_s << "#each line contains information about all clusters at critical point" << endl;
    fout_s << "#cluster size is measured by number of sites in it" << endl;

    ofstream fout_b(filename_b);
    // JSON formated header
    fout_b << header_info.str();
    fout_b << "#each line is an independent realization" << endl;
    fout_b << "#each line contains information about all clusters at critical point" << endl;
    fout_b << "#cluster size is measured by number of bonds in it" << endl;


    value_type counter{};
    for(value_type i{} ; i != ensemble_size ; ++i){

        sp.reset();

        bool successful = false;
        auto t_start = std::chrono::system_clock::now();
        counter = 0;
        while (true){
            successful = sp.occupy();
            if(successful) {
                if(sp.detectWrapping()){
                    vector<value_type> site, bond;
                    sp.get_cluster_info(site, bond);

                    for(value_type j{}; j != site.size(); ++j){
                        fout_s << site[j] << ',';
                        fout_b << bond[j] <<',';
                    }


                    fout_s << endl;
                    fout_b << endl;

                    break;
                }

                ++counter;
            }
            if(counter >= length_squared){ // length_squared is the number of site
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


    fout_b.close();
    fout_s.close();
}


void weighted_relabeling_test(int argc, char** argv) {
    value_type length = atoi(argv[1]);
//    value_type ensemble_size = atoi(argv[2]);


//    cout << "length " << length << " ensemble_size " << ensemble_size << endl;


//    srand(4);    // seeding
//    srand(atoi(argv[2]));    // seeding from cmdl

//    value_type length = 500;
    value_type length_squared = length * length;
    SitePercolation_ps_v8 sp(length, true);
    double agv_time{};
    size_t avg_count{};
    bool successful{false};
    cout << "sp.maxIterationLimit() " << sp.maxIterationLimit() << endl;
    size_t size = 1;
    for(size_t e{}; e < size; ++e) {
        cout << "step " << e << " of " << size << endl;
//        srand(e);    // seeding
        sp.reset();
        for (size_t j{}; j < sp.maxIterationLimit(); ++j) {
            successful = sp.occupy();
            if (successful) {
                auto index = sp.lastPlacedSite();

//            cout << j << " th site" << index << endl;
//            sp.viewSiteByID();
//            sp.viewClusterExtended();
            } else {
                cout << "False at " << j << endl;
            }
        }
//    sp.viewSiteByID();
//    sp.viewClusterExtended();
        cout << "number of relabeling done " << sp.relabeling_count() << endl;
        cout << "time required " << sp.get_relabeling_time() << endl;
        agv_time += sp.get_relabeling_time();
        avg_count += sp.relabeling_count();
    }
    cout << "avg number of relabeling done " << avg_count / size << endl;
    cout << "avg time required " << agv_time / size << endl;

}

void bond_percolation(int argc, char** argv) {
    value_type length = atoi(argv[1]);
    value_type ensemble_size = atoi(argv[2]);

//    value_type length = 5;
    BondPercolation_pb_v0 bp(length);
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
#P1(p,L) = Order parameter = (number of sites in largest cluster) / (total number of sites)
#P2(p,L) = Order parameter = (number of sites in spanning or wrapping cluster) / (total number of sites)
#C(p,L) = Specific heat = -T dH/dT
#X(p,L) = Susceptibility = dP/dp
#u_i = (number of sites in the i-th cluster) / (total number of sites))***";

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

void bond_percolation_wrapping(){
    value_type length = 5;
    BondPercolation_pb_v0 bp(length, true);

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    for(size_t i{}; i < 50; ++i) {
//        cout << "Iteration " << i << endl;
        bp.reset();
//        bp.viewByRelativeIndex();
        bool successful = false;
        auto t_start = std::chrono::system_clock::now();
        size_t counter = 0;
        while (true) {
            successful = bp.occupy();
            if (successful) {
//                cout << counter << " th bond " << bp.lastPlacedBond() << endl;
                if (bp.detectWrapping()) {
                    cout << bp.occupationProbability() << endl;
//                    bp.viewByRelativeIndex();

                    break;
                }

                ++counter;
            }
            if (counter >= twice_length_squared) { // length_squared is the number of site
                break;
            }
        }
    }
}

void percolation_wrapping_and_jump(int argc, char **argv) {
    value_type length = atoi(argv[1]);
    value_type ensemble_size = atoi(argv[2]);


    cout << "length " << length << " ensemble_size " << ensemble_size << endl;

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    SitePercolation_ps_v9 lattice_percolation(length, true);

    ostringstream header_info;
    header_info << "{"
                << "\"length\":" << length
                << ", \"ensemble_size\":" << ensemble_size
                << ", \"signature\":\"" << lattice_percolation.getSignature() << "\""
                << "}" ;

    string tm = currentTime();

    string filename_s = lattice_percolation.getSignature() + "_cluster_by_site_" + to_string(length) + '_' + tm;
    string filename_b = lattice_percolation.getSignature() + "_cluster_by_bond_" + to_string(length) + '_' + tm;
    string filename_critical = lattice_percolation.getSignature() + "_critical_" + to_string(length) + '_' + tm;
    filename_s += ".csv";
    filename_b += ".csv";
    filename_critical += ".csv";

    string filename = lattice_percolation.getSignature() + "_entropy-jump_" + to_string(length) + '_' + tm;
    filename += ".csv";

    ofstream fout_jump(filename);
    // JSON formated header
    fout_jump << '#' << header_info.str() << endl;
    fout_jump << "#each line is an independent realization" << endl;
    fout_jump << "#each line contains information about all clusters at critical point" << endl;
    fout_jump << "#cluster size is measured by number of sites in it" << endl;

    ofstream fout_s(filename_s);
    // JSON formated header
    fout_s << '#' << header_info.str() << endl;
    fout_s << "#each line is an independent realization" << endl;
    fout_s << "#each line contains information about all clusters at critical point" << endl;
    fout_s << "#cluster size is measured by number of sites in it" << endl;

    ofstream fout_b(filename_b);
    // JSON formated header
    fout_b << '#' << header_info.str() << endl;
    fout_b << "#each line is an independent realization" << endl;
    fout_b << "#each line contains information about all clusters at critical point" << endl;
    fout_b << "#cluster size is measured by number of bonds in it" << endl;

    ofstream fout_critical(filename_critical);
    fout_critical << '#' << header_info.str() << endl;
    fout_critical << "#critical occupation probability or pc" << endl;
    fout_critical << "#<pc>" << endl;

    value_type counter{};
    for(value_type i{} ; i != ensemble_size ; ++i){

        lattice_percolation.reset();

        bool successful = false;
        auto t_start = std::chrono::system_clock::now();
        counter = 0;
        while (true){
            successful = lattice_percolation.occupy();
            if(successful) {
                lattice_percolation.entropy();
                lattice_percolation.jump();
                if(lattice_percolation.detectWrapping()){
                    fout_critical << lattice_percolation.occupationProbability() << endl;
                    vector<value_type> site, bond;

                    lattice_percolation.get_cluster_info(site, bond);

                    for(value_type j{}; j != site.size(); ++j){
                        fout_s << site[j] << ',';
                    }
                    for(value_type j{}; j != bond.size(); ++j){
                        fout_b << bond[j] <<',';
                    }


                    fout_s << endl;
                    fout_b << endl;

                    break;
                }

                ++counter;
            }
            if(counter >= lattice_percolation.maxIterationLimit()){ // twice_length_squared is the number of bonds
                break;
            }
        }
        fout_jump << lattice_percolation.largestEntropyJump() << "," << lattice_percolation.largestEntropyJump_pc() << endl;
        {
            auto t_end = std::chrono::system_clock::now();
            cout << "Iteration " << i
                 //                 << " . Thread " << std::this_thread::get_id()
                 << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << endl;
        }
//        cout << "Relabeling time " << lattice_percolation.get_relabeling_time() << endl;
    }

    fout_b.close();
    fout_s.close();
    fout_critical.close();
    fout_jump.close();
}


void entropyJumps(int argc, char** argv){
    value_type length = atoi(argv[1]);
    value_type ensemble_size = atoi(argv[2]);


    cout << "length " << length << " ensemble_size " << ensemble_size << endl;

    value_type length_squared = length*length;
    value_type twice_length_squared = 2 * length_squared;

    BondPercolation_pb_v0 bp(length, true);

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


void site_percolation_new_version(int argc, char** argv){
    value_type length = 1000;
    length = atoi(argv[1]);
    SitePercolation_ps_v9 lattice_percolation(length);
//    lattice_percolation.reset();
    value_type j{};
    while(lattice_percolation.occupy()) {
        j+=1;
//        cout << j << "-th site " << lattice_percolation.lastPlacedBond() << endl;
//        if(lattice_percolation.detectWrapping()){
//            cout << "Wrapping " << lattice_percolation.occupationProbability() << endl;
//            break;
//        }

//        lattice_percolation.viewLatticeByID();
//        lattice_percolation.viewSiteByRelativeIndex();
//        lattice_percolation.viewClusterExtended();

        if (j >= lattice_percolation.maxIterationLimit()){
            break;
        }
    }
    cout << j << "-th site " << lattice_percolation.lastPlacedSite() << endl;
//    lattice_percolation.viewSiteByID();
//    lattice_percolation.viewSiteByRelativeIndex();
//    lattice_percolation.viewByRelativeIndex();
//    lattice_percolation.viewClusterExtended();
    cout << "relabeling count " << lattice_percolation.relabeling_count() << endl;
    cout << "time " << lattice_percolation.get_relabeling_time() << endl;
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

//    cout << "Seed " << seed << endl;


//    cluster_size_threaded(argc, argv);
//    cluster_size(argc, argv);
//    measure_entropy_by_site(argc, argv);
//    measure_clusters(argc, argv);
//    weighted_relabeling_test(argc, argv);

//    bond_percolation(argc, argv);

//    percolation_wrapping_and_jump(argc, argv);
//    percolation_wrapping_and_jump();
//    entropyJumps(argc, argv);
    site_percolation_new_version(argc, argv);

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


