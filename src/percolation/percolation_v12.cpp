//
// Created by shahnoor on 11/22/19.
//

#include <algorithm>
#include "percolation_v12.h"

using namespace std;

/**
 *
 * @param length
 */
Percolation_v12::Percolation_v12(int length):_length{length} {
    if(_length < 2)
        _length = 2;
    _lattice = Lattice_v12(length);


    _max_number_of_sites = static_cast<size_t>(_length * _length);
    _max_number_of_bonds = 2*_max_number_of_sites;
}


void Percolation_v12::setRandomState(size_t seed, bool generate_seed) {
//    size_t seed = 0;
//    cerr << "automatic seeding is commented : line " << __LINE__ << endl;
    _random_state = seed;
    if(generate_seed) {
        std::random_device _rd;
        _random_state = _rd();
    }else{
        cerr << "generate_seed = false : line " << __LINE__ << endl;
    }
    _random.seed(_random_state); // seeding
    cout << "seeding with " << _random_state << endl;
}


size_t Percolation_v12::getRandomState() {
//    size_t seed = 0;
//    cerr << "automatic seeding is commented : line " << __LINE__ << endl;
//    std::random_device _rd;
//    auto seed = _rd();
//    _random.seed(seed); // seeding
//    cout << "seeding with " << seed << endl;
    return _random_state;
}

void Percolation_v12::viewCluster() {
    cout << "Percolation_v12::viewCluster" << endl;
    for(size_t i{}; i < _clusters.size(); ++i){
        if(_clusters[i].empty()) continue;
        cout << "cluster[" << i << "] {" << endl;

        auto bonds = _clusters[i].getBondIDs();
        auto sites = _clusters[i].getSiteIDs();
        cout << "  Sites (" << sites.size()  << ") : {";
        for(auto s: sites){
            cout << s << ",";
        }
        cout << "}" << endl;
        cout << "  Bonds (" << bonds.size()  << ") : {";
        for(auto s: bonds){
            cout << s << ",";
        }
        cout << "}" << endl;
        cout << "}" << endl;
    }
}

void Percolation_v12::reset() {
    _clusters.clear();
    _number_of_bonds_in_the_largest_cluster = 0;
    _number_of_sites_in_the_largest_cluster = 0;
}

size_t Percolation_v12::numberOfBondsInTheLargestCluster() {
    return _number_of_bonds_in_the_largest_cluster;
}

size_t Percolation_v12::numberOfSitesInTheLargestCluster() {
    return _number_of_sites_in_the_largest_cluster;
}




