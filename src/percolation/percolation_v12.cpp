//
// Created by shahnoor on 11/22/19.
//

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


void SqLatticeRegularSite::init() {
    //

    // activate bonds and initialize cluster
    _clusters.resize(maxBonds());
    auto bonds = _lattice.getBonds();
    for(size_t i{}; i < _clusters.size(); ++i){

        _clusters[i].addBond()
    }

}
