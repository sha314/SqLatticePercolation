//
// Created by shahnoor on 12/4/2017.
//

#include "percolation_fractal.h"

using namespace std;


/**
 *
 * @param L_start
 * @param L_end
 * @param ensembleSize
 * @param df
 */
void moment_data(value_type  L_start, value_type  L_end, value_type ensembleSize, double df){
    string filename = "conservation of moment_";
    filename += currentTime();
    filename += "_.txt";
    ofstream fout(filename);
    fout << "#Moment for different length" << endl;
    fout << "#fractal dimension used, df = " << df << endl;
    fout << "#ensemble size" << endl << endl;
    fout << ensembleSize << endl << endl;
    fout << "#<Length>\t<M_{d_f}" << endl;
    double m{};
    clock_t t_outer;
    for(value_type L=L_start; L <= L_end ; L += 10){
        t_outer = clock();
        cout << "Length " << L << " : " << endl;
        m = moment_for_percolation_by_site(ensembleSize, L, df);
        fout << L << '\t' << m << endl;
        cout << "Length " << L << " . Time taken "
             << getFormattedTime((clock() - t_outer) / double(CLOCKS_PER_SEC)) << endl;
    }
    fout.close();
}



/**
 *
 * @param L_start
 * @param L_end
 * @param ensembleSize
 * @param df
 */
void spanning_data_for_lengths(value_type L_start, value_type L_end, value_type ensembleSize, value_type impure_site){

    clock_t t_outer;
    for(value_type L=L_start; L <= L_end ; L += 10){
        t_outer = clock();
        cout << "Length " << L << " : " << endl;
        spanning_data(ensembleSize, L, impure_site);
        cout << "Length " << L << " . Total Time taken "
             << getFormattedTime((clock() - t_outer) / double(CLOCKS_PER_SEC)) << endl;
    }
}



/**
 *
 * Moment = sum(i, x_i)
 *  here, x_i = number of site in the i-th cluster
 * @param ensemble_size
 * @param length
 * @param d_f               Fractal Dimension. d_f = 2.092 obtained
 * @return
 */
double moment_for_percolation_by_site(value_type ensemble_size, value_type length, double d_f){
    SitePercolation_ps_v8 sp(length);
    sp.periodicityFlag(false);

    clock_t t;
    vector<value_type> vec;
    double M{}, M_avg{};
    string filename = "cluster_length_data_L_";
    filename += to_string(length);
    filename += '_';
    filename += currentTime();
    filename += "_.txt";
    ofstream fout(filename);
    fout << "#cluster length data" << endl;
    fout << "#in site percolation" << endl;
    fout << "#number of sites is the cluster length" << endl;
    fout << "#Length = " << endl << length << endl;
    fout << "#data..." << endl;
//    value_type _possible_neighbors{};     // number of occupied site in a lattice
    for(value_type i{} ; i != ensemble_size ; ++i){
        t = clock();
        sp.reset();
        sp.placeSiteUntilFirstSpanning_v2();
        vec = sp.number_of_site_in_clusters();
        M = 0;
//        _possible_neighbors = 0;

        for(auto a: vec){
//            M += pow(a, d_f);
            fout << a << ' ';   // one space is enough to separate the data
            M += pow(a / double(length * length), d_f);
//            _possible_neighbors += a;             //
        }
        fout << endl;
        M_avg += M/ensemble_size;
        cout << "\t\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }
    fout.close();
    return M_avg;
}

/**
 * Write number of sites and number of bonds to files
 *
 * @param ensemble_size -> number of times to iterate
 * @param length        -> length of the lattice
 * @return
 */
void spanning_data(value_type ensemble_size, value_type length, value_type impure_sites){
    SitePercolation_ps_v8 sp(length, impure_sites);
    sp.periodicityFlag(false);

    cout << "Signature " << sp.signature << endl;
    clock_t t;
    vector<value_type> nos, nob;

    string common;
    common += to_string(length);
    common += '_';
    common += currentTime();;
    common += ".txt";

    string site_filename = "cluster_length_data_sites_L_" + common;
    string bond_filename = "cluster_length_data_bonds_L_" + common;
    string lattice_filename = "lattice_L_" + common;
    string time_filename = "time_L_" + common;

    ofstream fout_site(site_filename);        // the output data file
    ofstream fout_bond(bond_filename);        // the output data file
    ofstream fout_lattice(lattice_filename);    // write the lattice to a file
    ofstream fout_spanning_cluster(time_filename);      // contains the time information of the spanning cluster

    /**
     * First uncomment line is the signature
     * Second uncomment line is the parameter value
     * Third and all other uncomment lines are data
     */
    fout_site << sp.signature << endl;
    fout_site << "#cluster length data in site percolation" << endl;
    fout_site << "#clustter length = number of sites" << endl;
    fout_site << "#length = " << endl << length << endl;
    fout_site << "#<data...>" << endl;

    fout_bond << sp.signature << endl;
    fout_bond << "#cluster length data in site percolation" << endl;
    fout_bond << "#clustter length = number of bonds" << endl;
    fout_bond << "#length = " << endl << length << endl;
    fout_bond << "#<data...>" << endl;

    fout_lattice << "#Structure of a lattice when spanning occurs for the first time" << endl;
    fout_lattice << "#only id of the clusters are printed and separated by spaces" << endl;
    fout_lattice << "#Length = " << endl << length << endl;
    fout_lattice << "#data..." << endl;
    fout_lattice << "#Each new data point starts with \'>\'" << endl;
    fout_lattice << "#Spanning cluster id = is written after \'>\' and a space character \' \'" << endl ;
    fout_lattice << "#Number of times > symbol occurs is the number of lattice present in one file" << endl;

    fout_spanning_cluster << "#Creation or Birth Time of spanning cluster" << endl;
    fout_spanning_cluster << "#Total Time spent = number of sites in the lattice" << endl;
    fout_spanning_cluster << "#Length = " << endl << length << endl;
    fout_spanning_cluster << "#<Birth Time of spanning cluster>\t<Id of the spanning cluster>\t"
            "<Total Time spent>\t<Number of Sites in Spanning cluster>\t"
            "<Number of bonds in Spanning cluster>" << endl;

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
        fout_lattice << '>' << ' ' << sp.firstSpanningClusterID()<<endl;    // indicates new lattice of same length
        fout_lattice << sp.getLatticeIDs() << endl;  // '-' indicates that a new iteration is begining
        fout_spanning_cluster << sp.birthTimeOfSpanningCluster() <<'\t'
                              << sp.firstSpanningClusterID() << '\t' << sp.numberOfOccupiedSite() << '\t'
                              << sp.numberOfSitesInTheSpanningClusters_v2() << '\t'
                              << sp.numberOfBondsInTheSpanningClusters_v2()
                              << endl;

        cout << "\t\tIteration " << i << " . Time " << (clock() - t) / double(CLOCKS_PER_SEC) << " sec" << endl;
    }

    fout_site.close();
    fout_bond.close();
    fout_lattice.close();
    fout_spanning_cluster.close();

}

