//
// Created by shahnoor on 3/21/21.
//

#include "run_v13.h"
#include "percolation_v13.h"
#include "percolation_rsbd_v13.h"
#include "percolation_v13_bond.h"

using namespace std;

void plus_equal(vector<double> &left_vec, const vector<double> &right_vec){
    if (left_vec.size() != right_vec.size()){
        cout << "Size mismatch : " << __FILE__ << " : " << __LINE__ << endl;
        cout << left_vec.size() << " != " << right_vec.size() << endl;
        return;
    }
    for(size_t i{}; i < right_vec.size(); ++i){
//        cout << left_vec[i] << " + " << right_vec[i];
        left_vec[i] += right_vec[i];
//        cout << " = " << left_vec[i] << endl;
    }

}

void run_v13_rsbd_L0(int length, int ensemble_size){
    std::cout << "length " << length << " ensemble_size " << ensemble_size << std::endl;

    size_t length_squared = length*length;
    size_t twice_length_squared = 2 * length_squared;

    SitePercolationL0_v13 percolation(length, 0, true);
    percolation.setRandomState(0, true);
//    percolation.init();




    // simulation starts here
    std::vector<double> pcs(ensemble_size), sites_pc(ensemble_size), bonds_pc(ensemble_size);
    value_type counter{};
    std::vector<double> entropy(percolation.maxIterationLimit());
    std::vector<double> order_wraping(percolation.maxIterationLimit()),
            order_largest(percolation.maxIterationLimit()),
            mean_cluster_sz(percolation.maxIterationLimit());

    vector<double> tmp;

    for(value_type i{} ; i != ensemble_size ; ++i){
        auto t_start = std::chrono::system_clock::now();

        percolation.reset();
//        SitePercolationL0_v13 percolation(length, 0, true);
//        percolation.viewCluster(1);
        percolation.run_once_v2();
//        percolation.viewCluster(1);

        pcs[i] = percolation.get_pc();
        sites_pc[i] = percolation.get_wrapping_cluster_site_count_at_pc();
        bonds_pc[i] = percolation.get_wrapping_cluster_bond_count_at_pc();

//        cout << "entropy " << endl;
        tmp = percolation.get_entropy_array();
        plus_equal(entropy, tmp);

//        cout << "order param " << endl;
        tmp = percolation.get_order_param_wrapping_array();
        plus_equal(order_wraping, tmp);

//        cout << "order " << endl;
        tmp = percolation.get_order_param_largest_array();
        plus_equal(order_largest, tmp);

        tmp = percolation.get_mean_cluster_size_array();
        plus_equal(mean_cluster_sz, tmp);

        auto t_end = std::chrono::system_clock::now();
        std::cout << "Iteration " << i
                  //                 << " . Thread " << std::this_thread::get_id()
                  << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << std::endl;

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
                << R"(,"En":)" << ensemble_size
                << R"(,"random_seed":)" << percolation.getRandomState()
                << R"(,"signature":")" << percolation.get_signature() << "\""
                << R"(,"datetime":")" << tm << "\""
                //                << R"(,"cols":)" << R"(["S","n_S"])"
                << "}";


    std::string extension = "_L" + std::to_string(length) + '_' + tm + ".txt";
    std::string filename_critical = percolation.get_signature() + "_critical" + extension;
    std::string filename_entropy_order_parameter = percolation.get_signature() + "_entropy-order" + extension;


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
    for(size_t i{}; i < percolation.maxIterationLimit(); ++i){

        fout << (i + 1) / double(percolation.maxIterationLimit()) << delimiter;
        fout << entropy[i] / double(ensemble_size) << delimiter;
        fout << order_largest[i] / double(ensemble_size) << delimiter;
        fout << order_wraping[i] / double(ensemble_size) << delimiter;
        fout << mean_cluster_sz[i] / double(ensemble_size) ;
        fout << std::endl;
    }
    fout.close();
}

/**
 * @brief Bond percolation testing
 * 
 * @param length 
 * @param ensemble_size 
 */
void run_v13_bond(int length, int ensemble_size){
    std::cout << "length " << length << " ensemble_size " << ensemble_size << std::endl;

    size_t length_squared = length*length;
    size_t twice_length_squared = 2 * length_squared;

    SitePercolationL0_v13 percolation(length, 0, true);
    percolation.setRandomState(0, true);
//    percolation.init();




    // simulation starts here
    std::vector<double> pcs(ensemble_size), sites_pc(ensemble_size), bonds_pc(ensemble_size);
    value_type counter{};
    std::vector<double> entropy(percolation.maxIterationLimit());
    std::vector<double> order_wraping(percolation.maxIterationLimit()),
            order_largest(percolation.maxIterationLimit()),
            mean_cluster_sz(percolation.maxIterationLimit());

    vector<double> tmp;

    for(value_type i{} ; i != ensemble_size ; ++i){
        auto t_start = std::chrono::system_clock::now();

        percolation.reset();
//        SitePercolationL0_v13 percolation(length, 0, true);
//        percolation.viewCluster(1);
        percolation.run_once_v2();
//        percolation.viewCluster(1);

        pcs[i] = percolation.get_pc();
        sites_pc[i] = percolation.get_wrapping_cluster_site_count_at_pc();
        bonds_pc[i] = percolation.get_wrapping_cluster_bond_count_at_pc();

//        cout << "entropy " << endl;
        tmp = percolation.get_entropy_array();
        plus_equal(entropy, tmp);

//        cout << "order param " << endl;
        tmp = percolation.get_order_param_wrapping_array();
        plus_equal(order_wraping, tmp);

//        cout << "order " << endl;
        tmp = percolation.get_order_param_largest_array();
        plus_equal(order_largest, tmp);

        tmp = percolation.get_mean_cluster_size_array();
        plus_equal(mean_cluster_sz, tmp);

        auto t_end = std::chrono::system_clock::now();
        std::cout << "Iteration " << i
                  //                 << " . Thread " << std::this_thread::get_id()
                  << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << std::endl;

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
                << R"(,"En":)" << ensemble_size
                << R"(,"random_seed":)" << percolation.getRandomState()
                << R"(,"signature":")" << percolation.get_signature() << "\""
                << R"(,"datetime":")" << tm << "\""
                //                << R"(,"cols":)" << R"(["S","n_S"])"
                << "}";


    std::string extension = "_L" + std::to_string(length) + '_' + tm + ".txt";
    std::string filename_critical = percolation.get_signature() + "_critical" + extension;
    std::string filename_entropy_order_parameter = percolation.get_signature() + "_entropy-order" + extension;


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
    for(size_t i{}; i < percolation.maxIterationLimit(); ++i){

        fout << (i + 1) / double(percolation.maxIterationLimit()) << delimiter;
        fout << entropy[i] / double(ensemble_size) << delimiter;
        fout << order_largest[i] / double(ensemble_size) << delimiter;
        fout << order_wraping[i] / double(ensemble_size) << delimiter;
        fout << mean_cluster_sz[i] / double(ensemble_size) ;
        fout << std::endl;
    }
    fout.close();
}

void run_v13_rsbd_L1(int length, int ensemble_size){
    std::cout << "length " << length << " ensemble_size " << ensemble_size << std::endl;

    size_t length_squared = length*length;
    size_t twice_length_squared = 2 * length_squared;

    SitePercolationL1_v13 percolation(length, 0, true);
    percolation.setRandomState(0, true);
//    percolation.init();




    // simulation starts here
    std::vector<double> pcs(ensemble_size), sites_pc(ensemble_size), bonds_pc(ensemble_size);
    value_type counter{};
    std::vector<double> entropy(percolation.maxIterationLimit());
    std::vector<double> order_wraping(percolation.maxIterationLimit()),
            order_largest(percolation.maxIterationLimit()),
            mean_cluster_sz(percolation.maxIterationLimit());

    vector<double> tmp;

    for(value_type i{} ; i != ensemble_size ; ++i){
        auto t_start = std::chrono::system_clock::now();

        percolation.reset();
//        SitePercolationL0_v13 percolation(length, 0, true);
//        percolation.viewCluster(1);
        percolation.run_once_v2();
//        percolation.viewCluster(1);

        pcs[i] = percolation.get_pc();
        sites_pc[i] = percolation.get_wrapping_cluster_site_count_at_pc();
        bonds_pc[i] = percolation.get_wrapping_cluster_bond_count_at_pc();

//        cout << "entropy " << endl;
        tmp = percolation.get_entropy_array();
        plus_equal(entropy, tmp);

//        cout << "order param " << endl;
        tmp = percolation.get_order_param_wrapping_array();
        plus_equal(order_wraping, tmp);

//        cout << "order " << endl;
        tmp = percolation.get_order_param_largest_array();
        plus_equal(order_largest, tmp);

        tmp = percolation.get_mean_cluster_size_array();
        plus_equal(mean_cluster_sz, tmp);

        auto t_end = std::chrono::system_clock::now();
        std::cout << "Iteration " << i
                  //                 << " . Thread " << std::this_thread::get_id()
                  << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << std::endl;

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
                << R"(,"En":)" << ensemble_size
                << R"(,"random_seed":)" << percolation.getRandomState()
                << R"(,"signature":")" << percolation.get_signature() << "\""
                << R"(,"datetime":")" << tm << "\""
                //                << R"(,"cols":)" << R"(["S","n_S"])"
                << "}";


    std::string extension = "_L" + std::to_string(length) + '_' + tm + ".txt";
    std::string filename_critical = percolation.get_signature() + "_critical" + extension;
    std::string filename_entropy_order_parameter = percolation.get_signature() + "_entropy-order" + extension;


//    std::ofstream fout_jump(filename);
//    // JSON formated header
//    fout_jump << '#' << header_info.str() << std::endl;
//    fout_jump << "#each line is an independent realization" << std::endl;
//    fout_jump << "#each line contains information about all clusters at critical point" << std::endl;
//    fout_jump << "#cluster size is measured by number of bonds in it" << std::endl;
//    fout_jump << "#<entropy jump>,<occupation probability at which the jump occurs>" << std::endl;

//    std::ofstream fout_s(filename_s);
//    // JSON formated header
//    fout_s << '#' << header_info.str() << std::endl;
//    fout_s << "#each line is an independent realization" << std::endl;
//    fout_s << "#each line contains information about all clusters at critical point" << std::endl;
//    fout_s << "#cluster size is measured by number of site_index_sequence in it" << std::endl;

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
    fout << "#<p>,<H(p,L)>,<P1(p,L)>,<P2(p,L)>,<S(p,L)>" << std::endl;
    fout << "#p = occupation probability" << std::endl;
    fout << "#H(p,L) = Entropy = sum( - u_i * log(u_i))" << std::endl;
    fout << "#P1(p,L) = Order parameter = (number of bonds in largest cluster) / (total number of bonds)" << std::endl;
    fout << "#P2(p,L) = Order parameter = (number of bonds in spanning or wrapping cluster) / (total number of bonds)" << std::endl;
    fout << "#C(p,L) = Specific heat = -T dH/dT" << std::endl;
    fout << "#X(p,L) = Susceptibility = dP/dp" << std::endl;
    fout << "#S(p,L) = mean cluster size = ratio of 2nd moment to first moment = (sum_s s)/(sum_s s^2)" << std::endl;
    fout << "#u_i = (number of bonds in the i-th cluster) / (total number of bonds)" << std::endl;
    for(size_t i{}; i < percolation.maxIterationLimit(); ++i){

        fout << (i + 1) / double(percolation.maxIterationLimit()) << delimiter;
        fout << entropy[i] / double(ensemble_size) << delimiter;
        fout << order_largest[i] / double(ensemble_size) << delimiter;
        fout << order_wraping[i] / double(ensemble_size) << delimiter;
        fout << mean_cluster_sz[i] / double(ensemble_size) ;
        fout << std::endl;
    }
    fout.close();
}


void run_v13_rsbd_L1_Ptilda_test(int length, int ensemble_size){
    std::cout << "length " << length << " ensemble_size " << ensemble_size << std::endl;

    size_t length_squared = length*length;
    size_t twice_length_squared = 2 * length_squared;

    SitePercolationL1_v13_Ptilda percolation(length, 0, true);
    percolation.setRandomState(0, true);
//    percolation.init();


    // simulation starts here
    std::vector<double> pcs(ensemble_size), sites_pc(ensemble_size), bonds_pc(ensemble_size);
    value_type counter{};
    std::vector<double> entropy(percolation.maxIterationLimit());
    std::vector<double> order_wraping(percolation.maxIterationLimit()),
            order_largest(percolation.maxIterationLimit());
    vector<double> tmp;

    for(value_type i{} ; i != ensemble_size ; ++i){
        auto t_start = std::chrono::system_clock::now();

        percolation.reset();
//        SitePercolationL0_v13 percolation(length, 0, true);
//        percolation.viewCluster(1);
        percolation.run_once_v2();
//        percolation.viewCluster(1);


        auto t_end = std::chrono::system_clock::now();
        std::cout << "Iteration " << i
                  //                 << " . Thread " << std::this_thread::get_id()
                  << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << std::endl;

//        cout << "Relabeling time " << lattice_percolation.get_relabeling_time() << endl;
    }
    auto Ptilda = percolation.get_Ptilda();

    std::string tm = getCurrentTime();
    const char delimiter = '\t';
    std::ostringstream header_info;
    header_info << "{"
                << R"("length":)" << length
                << R"(,"ensemble_size":)" << ensemble_size
                << R"(,"En":)" << ensemble_size
                << R"(,"random_seed":)" << percolation.getRandomState()
                << R"(,"signature":")" << percolation.get_signature() << "\""
                << R"(,"datetime":")" << tm << "\""
                //                << R"(,"cols":)" << R"(["S","n_S"])"
                << "}";


    std::string extension = "_L" + std::to_string(length) + '_' + tm + ".txt";

    std::string filename_Ptilda = percolation.get_signature() + extension;
    cout << "got here " << __LINE__ << endl;


    std::ofstream fout_critical(filename_Ptilda);
    fout_critical << '#' << header_info.str() << std::endl;
    fout_critical << "#Ptilda=Probability that at least one of the four neighbor site is empty=p(1-p)/4 theoritically?" << std::endl;
    fout_critical << "#Ptilda= is an array of small numbers. so let's do Ptilda*L^2" << std::endl;
    fout_critical << "#<p>,<Ptilda>" << std::endl;
    for(size_t i{}; i < Ptilda.size(); ++i){
        fout_critical << i/double(length_squared) << delimiter
                      << Ptilda[i]/double(ensemble_size) << std::endl;
    }
    fout_critical.close();
}


void run_v13_bond_percolation(int length, int ensemble_size){
    std::cout << "length " << length << " ensemble_size " << ensemble_size << std::endl;

    size_t length_squared = length*length;
    size_t twice_length_squared = 2 * length_squared;

    BondPercolation_v13 percolation(length, 0, false);
    // BondPercolationExplosive_v13 percolation(length, 5, 0, false);

    percolation.setRandomState(0, true);
//    percolation.init();




    // simulation starts here
    std::vector<double> pcs(ensemble_size), sites_pc(ensemble_size), bonds_pc(ensemble_size);
    value_type counter{};
    std::vector<double> entropy(percolation.maxIterationLimit());
    std::vector<double> order_wraping(percolation.maxIterationLimit()),
            order_largest(percolation.maxIterationLimit()),
            mean_cluster_sz(percolation.maxIterationLimit());

    vector<double> tmp;

    for(value_type i{} ; i != ensemble_size ; ++i){
        auto t_start = std::chrono::system_clock::now();

        percolation.reset();
//        SitePercolationL0_v13 percolation(length, 0, true);
//        percolation.viewCluster(1);
        percolation.run_once();
//        percolation.viewCluster(1);

        pcs[i] = percolation.get_pc();
        sites_pc[i] = percolation.get_wrapping_cluster_site_count_at_pc();
        bonds_pc[i] = percolation.get_wrapping_cluster_bond_count_at_pc();

//        cout << "entropy " << endl;
        tmp = percolation.get_entropy_array();
        plus_equal(entropy, tmp);

//        cout << "order param " << endl;
        tmp = percolation.get_order_param_wrapping_array();
        plus_equal(order_wraping, tmp);

//        cout << "order " << endl;
        tmp = percolation.get_order_param_largest_array();
        plus_equal(order_largest, tmp);

        // tmp = percolation.get_mean_cluster_size_array();
        // plus_equal(mean_cluster_sz, tmp);

        auto t_end = std::chrono::system_clock::now();
        std::cout << "Iteration " << i
                  //                 << " . Thread " << std::this_thread::get_id()
                  << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << std::endl;

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
                << R"(,"En":)" << ensemble_size
                << R"(,"random_seed":)" << percolation.getRandomState()
                << R"(,"signature":")" << percolation.get_signature() << "\""
                << R"(,"datetime":")" << tm << "\""
                //                << R"(,"cols":)" << R"(["S","n_S"])"
                << "}";


    std::string extension = "_L" + std::to_string(length) + '_' + tm + ".txt";
    std::string filename_critical = percolation.get_signature() + "_critical" + extension;
    std::string filename_entropy_order_parameter = percolation.get_signature() + "_entropy-order" + extension;


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
    for(size_t i{}; i < percolation.maxIterationLimit(); ++i){

        fout << (i + 1) / double(percolation.maxIterationLimit()) << delimiter;
        fout << entropy[i] / double(ensemble_size) << delimiter;
        fout << order_largest[i] / double(ensemble_size) << delimiter;
        fout << order_wraping[i] / double(ensemble_size);
        // fout << mean_cluster_sz[i] / double(ensemble_size) ;
        fout << std::endl;
    }
    fout.close();
}

void run_v13_bond_percolation_explosive(int length, int ensemble_size, int M_value){
    std::cout << "length " << length << " ensemble_size " << ensemble_size << std::endl;

    size_t length_squared = length*length;
    size_t twice_length_squared = 2 * length_squared;

    BondPercolationExplosive_v13 percolation(length, M_value, 0, false);

    // percolation.setRandomState(0, false);
//    percolation.init();

    // simulation starts here
    std::vector<double> pcs(ensemble_size), sites_pc(ensemble_size), bonds_pc(ensemble_size);
    value_type counter{};
    std::vector<double> entropy(percolation.maxIterationLimit());
    std::vector<double> order_wraping(percolation.maxIterationLimit()),
            order_largest(percolation.maxIterationLimit()),
            mean_cluster_sz(percolation.maxIterationLimit()), 
            dHs(ensemble_size),
            dPs(ensemble_size);

    vector<double> tmp;

    for(value_type i{} ; i != ensemble_size ; ++i){
        auto t_start = std::chrono::system_clock::now();

        percolation.reset();
        // percolation.check_system();
//        SitePercolationL0_v13 percolation(length, 0, true);
//        percolation.viewCluster(1);
        percolation.run_once();
//        percolation.viewCluster(1);

        pcs[i] = percolation.get_pc();
        sites_pc[i] = percolation.get_wrapping_cluster_site_count_at_pc();
        bonds_pc[i] = percolation.get_wrapping_cluster_bond_count_at_pc();
        dHs[i] = percolation.jump_entropy();
        dPs[i] = percolation.jump_largest_cluster();

//        cout << "entropy " << endl;
        tmp = percolation.get_entropy_array();
        plus_equal(entropy, tmp);

//        cout << "order param " << endl;
        tmp = percolation.get_order_param_wrapping_array();
        plus_equal(order_wraping, tmp);

//        cout << "order " << endl;
        tmp = percolation.get_order_param_largest_array();
        plus_equal(order_largest, tmp);

        // tmp = percolation.get_mean_cluster_size_array();
        // plus_equal(mean_cluster_sz, tmp);


        auto t_end = std::chrono::system_clock::now();
        std::cout << "Iteration " << i
                  //                 << " . Thread " << std::this_thread::get_id()
                  << " . Elapsed time " << std::chrono::duration<double>(t_end - t_start).count() << " sec" << std::endl;

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
                << R"(,"En":)" << ensemble_size
                << R"(,"random_seed":)" << percolation.getRandomState()
                << R"(,"signature":")" << percolation.get_signature() << "\""
                << R"(,"datetime":")" << tm << "\""
                //                << R"(,"cols":)" << R"(["S","n_S"])"
                << "}";


    std::string extension = "_L" + std::to_string(length) + '_' + tm + ".txt";
    std::string filename_critical = percolation.get_signature() + "_critical" + extension;
    std::string filename_entropy_order_parameter = percolation.get_signature() + "_entropy-order" + extension;


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
    for(size_t i{}; i < percolation.maxIterationLimit(); ++i){

        fout << (i + 1) / double(percolation.maxIterationLimit()) << delimiter;
        fout << entropy[i] / double(ensemble_size) << delimiter;
        fout << order_largest[i] / double(ensemble_size) << delimiter;
        fout << order_wraping[i] / double(ensemble_size);
        // fout << mean_cluster_sz[i] / double(ensemble_size) ;
        fout << std::endl;
    }
    fout.close();



    std::string filename_jump = percolation.get_signature() + "_jump_" + extension;

    ofstream fout_jump(filename_jump);
    fout_jump << '#' << header_info.str() << std::endl;
    fout_jump << "#<entropy jump><order parameter jump>" << endl;
    for(size_t i{}; i < dHs.size(); ++i){
        fout_jump << dHs[i] << "\t" << dPs[i] << endl;
    }
    fout_jump.close();
}
