//
// Created by shahnoor on 10/31/2017.
//

#include "percolation.h"


#include <cstdlib>
#include <climits>
#include <unordered_set>

#include "percolation.h"

#include "../util/printer.h"
#include "site_position.h"


using namespace std;

/**
 *  variables
 */
double total_time_for_find_index_for_placing_new_bonds{};


/**
 *
 * @param length
 */
BondPercolation_pb_v0::BondPercolation_pb_v0(value_type length, bool periodicity)
        : SqLatticePercolation(length) {
    std::cout << "Constructing BondPercolation_pb_v0 object : line " << __LINE__ << endl;

    _periodicity = periodicity;
    _index_sequence_position = 0;
    _lattice = SqLattice(length, false,
                         true);   // since it is a site percolation all bonds will be activated by default

    index_sequence.reserve(_max_bonds);

    // there are 2*L*L cluster initially but only clusters with size larger than 1 should be counted
    _clusters = vector<Cluster_v2>();

    initialize_index_sequence();
    initialize();

//    initialize_indices();

    randomize();  // randomize the untouched_site_indices
}

/**
 * Initialize the Class Object
 */
void BondPercolation_pb_v0::initialize() {
    if (_length <= 2) {
        /*
         * Because if _length=2
         * there are total of 4 distinct bond. But it should have been 8, i.e, (2 * _length * _length = 8)
         */
        cerr << "_length <= 2 does not satisfy _lattice properties for percolation : line" << __LINE__ << endl;
    }

    // to improve performence
    number_of_sites_to_span.reserve(_length_squared);
    number_of_bonds_to_span.reserve(_length_squared);

//    _top_edge.reserve(_length);
//    _bottom_edge.reserve(_length);
//    _left_edge.reserve(_length);
//    _right_edge.reserve(_length);

    randomized_index_sequence = index_sequence;
}


/**
 * Called only once when the object is constructed for the first time
 */
void BondPercolation_pb_v0::initialize_index_sequence() {
    value_type row{}, col{};

    if (!_periodicity) {
        value_type limit = _length_squared - 1; // since if row==4 and col ==4 we don't want that bond
        for (value_type i{}; i < limit; ++i) {
            /*if (row == max_index && col == max_index){
                continue; // it is the last step. so decrease the loop limit by 1 will do the job
            }else */
            if (row == max_index) {
                index_sequence.push_back({BondType::Horizontal, row, col});
            } else if (col == max_index) {
                index_sequence.push_back({BondType::Vertical, row, col});
            } else {
                index_sequence.push_back({BondType::Horizontal, row, col});
                index_sequence.push_back({BondType::Vertical, row, col});
            }

            ++col;
            if (col == _length) {
                col = 0;
                ++row;
            }
        }
    } else {
        for (value_type i{}; i < _length_squared; ++i) {

            index_sequence.push_back({BondType::Horizontal, row, col});
            index_sequence.push_back({BondType::Vertical, row, col});

            ++col;
            if (col == _length) {
                col = 0;
                ++row;
            }
        }
    }
}

/**
 *
 */
void BondPercolation_pb_v0::initialize_indices() {
    indices = vector<value_type>(index_sequence.size());
    for (value_type i{}; i != indices.size(); ++i) {
        indices[i] = i; // assign index first
    }
}


/**
 * Randomize the site placing order
 * Takes 3.031000 sec for 1000 x 1000 sites
 */
void BondPercolation_pb_v0::randomize() {
    randomized_index_sequence = index_sequence;
    value_type len = randomized_index_sequence.size();
    value_type j{};
    BondIndex tmp;
    for (value_type i{}; i != len; ++i) {
        // select a j from the array. which must come from the ordered region
        j = i + std::rand() % (len - i);
        // perform the swapping with i-th and j-th value
        tmp = randomized_index_sequence[i];
        randomized_index_sequence[i] = randomized_index_sequence[j]; // todo problematic
        randomized_index_sequence[j] = tmp;
    }

//    cout << "Index sequence : " << randomized_index_sequence.size()
//         << " : " << randomized_index_sequence << endl;

}


/**
 * Reset all calculated values and then call initiate()
 * to initiallize for reuse
 *
 * caution -> it does not erase _calculation_flags, for it will be used for calculation purposes
 */
void BondPercolation_pb_v0::reset() {
    _total_number_of_active_bonds = 0;
    _clusters.clear();
    _entropy.clear();
    _number_of_occupied_bonds = 0;
    number_of_bonds_to_span.clear();
    _sites_required_to_min_span = 0;

    initialize();
}


/**
 * Create a custome configuration for the lattice
 */
void BondPercolation_pb_v0::configure(std::vector<BondIndex> bond_indices) {
    cout << "Entry -> configure() : line " << endl;

}


///**
// *
// */
//void BondPercolation_pb_v0::calculate_occupation_probability() {
//    // occupation probability = number of site present / total number of site
//    _number_of_occupied_bonds.push_back(_total_number_of_active_bonds);
//}

/**
 *
 */
void BondPercolation_pb_v0::calculate_entropy() {
    double p_i{}; // probability for the cluster
    double S{}; // entropyDistribution
    double counter{}; // counts how many elements _clusters contains
    for (value_type i{}; i != _clusters.size(); ++i) {
        counter += _clusters[i].numberOfBonds();
        p_i = _clusters[i].numberOfBonds() / max_cluster_length;
        S += p_i * log(p_i);
    }
    // for cluster with numberOfBonds 1
    S += ((max_cluster_length - counter) / max_cluster_length) * log(1.0 / max_cluster_length);
    S *= -1;    // since S = - p_i * log(p_i)
    _entropy.push_back(S);
}

/**
 *
 * @param cf
 */
void BondPercolation_pb_v0::calculationFlags(std::vector<CalculationFlags> cf) {
    if (cf.empty()) {
        cout << "No calculation flag available : line " << __LINE__ << endl;
    }

    _calculation_flags.insert(_calculation_flags.begin(), cf.begin(), cf.end());

}


/**
 *
 * @param clstr
 * @param id
 */
void BondPercolation_pb_v0::relabel_sites(const Cluster_v2 &clstr, int id) {

    const vector<Index> sites = clstr.getSiteIndices();
    for (auto a: sites) {
        _lattice.getSite(a).groupID(id);
    }

}

/**
 *
 * @param clstr
 * @param id
 */
void BondPercolation_pb_v0::relabel_bonds(const Cluster_v2 &clstr, int id) {
    vector<BondIndex> bonds = clstr.getBondIndices();
    for (auto a: bonds) {
        _lattice.getBond(a).groupID(id);
    }
}


/**
 * All site placing method in one place
 *
 * @return true if operation is successfull
 */
bool BondPercolation_pb_v0::occupy() {
    if (_index_sequence_position >= _max_bonds) {
        return false;
    }

    value_type v = placeBond_v0();

    return v != ULONG_MAX;
}


/**
 * Place one site at a time. Actually activates one site and 4 bonds at a time.
 * Also label and relabel each bond and site as per their corresponding cluster
 *
 * @return merged_cluster_index -> the index of _cluster
 *                                      where new site and bond indicis are placed
 *                                      or the merged cluster index
 *
 */
value_type BondPercolation_pb_v0::placeBond_v0() {

    if (_index_sequence_position == randomized_index_sequence.size()) {
        return ULONG_MAX;// unsigned long int maximum value
    }

    BondIndex current_bond = randomized_index_sequence[_index_sequence_position];
    _last_placed_bond = current_bond;
//    cout << "placing bond " << current_bond << " : " << __LINE__ << endl;

    _lattice.activateBond(current_bond);

    ++_number_of_occupied_bonds;
    ++_index_sequence_position;

    // find the bonds for this site
    vector<BondIndex> bonds;
    vector<Index> sites;

//    connection_v1(current_bond, sites, bonds); // problem
    connection_v2(current_bond, sites, bonds);
//    cout << "Found bonds : " << bonds << endl;
//    cout << "Found sites : " << sites << endl;
    mark_sites(sites);

    // find one of hv_bonds in _clusters and add ever other value to that place. then erase other position
    set<value_type> found_index_set = find_index_for_placing_new_bonds_v1(bonds);


//    cout << "Found indices " << found_index_set << endl;

//    subtract_entropy_for(found_index_set);  // tracking entropy change
    value_type merged_cluster_index = manage_clusters_v1(
            found_index_set, sites, current_bond
    );
//    add_entropy_for(merged_cluster_index); // tracking entropy change

    // running tracker
//    track_numberOfBondsInLargestCluster(); // tracking number of bonds in the largest cluster

    return merged_cluster_index;

}


/**
 * For marking sites that is a part of cluster with size > 1
 * @param sites
 */
void BondPercolation_pb_v0::mark_sites(vector<Index> &sites) {
    for (auto it{sites.begin()}; it < sites.end(); ++it) {
        if (_lattice.getSite(*it).groupID() >= 0) {
            // if id >= 0 then it is counted in the cluster so erase it
//            sites.erase(it); // commenting for now.
        } else {
            _lattice.getSite(*it).groupID(0); // else set id == 0
        }
    }
}

/***********************
 * new methods
 *
 ********************************************/


/**
 * Functions that will give correct value
 *      entropy()
 *      entropy_v3()
 * @param found_index_set
 * @param sites
 * @param site
 * @return
 */
value_type BondPercolation_pb_v0::manage_clusters_v1(
        const set<value_type> &found_index_set,
        vector<Index> &sites,
        BondIndex &bond
) {
    vector<value_type> found_index(found_index_set.begin(), found_index_set.end());
    value_type merged_cluster_index{};

    if (found_index_set.size() > 0) {
        unsigned long &base = found_index[0];
        _clusters[base].addBondIndex(bond);
        _lattice.getBond(bond).groupID(_clusters[base].ID()); // relabeling for 1 site

        // put_values_to_the_cluster new values in the 0-th found index
        _clusters[base].insert(sites); // todo
        // merge clusters with common values from all other cluster
        int tmp_id;
        for (value_type k{1}; k != found_index.size(); ++k) {

            tmp_id = _clusters[base].ID();

            // erase the redundant cluster ids
            int id_to_be_deleted = _clusters[found_index[k]].ID();
            _cluster_index_from_id.erase(id_to_be_deleted); // first erase previous keys

            // perform relabeling on the sites
            relabel_bonds(_clusters[found_index[k]], tmp_id);


            // since we use cluster id to relabel cluster when merging, cluster also need to be updated
            _clusters[found_index[k]].ID(tmp_id);

            /// since all cluster will eventually get merged to cluster with
            /// index found_index[0] whatever the id is index must be found_index[0]

            // store values of other found indices to the cluster
            _clusters[base].insert(_clusters[found_index[k]]);

            // delete the merged cluster
            auto it = _clusters.begin() + found_index[k];
            _clusters.erase(it);
            // after erasing found_index cannot remain same. come on.
            // If you don't do this we will be deleting wrond clusters
            for (value_type m{k + 1}; m != found_index.size(); ++m) {
                // only the values of found_index that will be used next Time must be reduced
                --found_index[m];
            }
//            cout << "merging cluster completed : line " << __LINE__ << endl;
        }
        merged_cluster_index = base;

    } else {
        // create new element for the cluster
        _clusters.push_back(Cluster_v2(_cluster_id));
        merged_cluster_index = _clusters.size() - 1;  // this new cluster index
//        _cluster_index_from_id[_cluster_id] = _clusters.size() - 1; // keeps track of cluster id and cluster index
        _cluster_index_from_id.insert(_cluster_id); // new version
        _lattice.getBond(bond).groupID(_cluster_id); // relabeling for 1 site
        _cluster_id++;  // increase the cluster id for next round
        _clusters.back().insert(sites);
        _clusters[merged_cluster_index].addBondIndex(bond);

    }


    // data for short cut calculation
    _index_last_modified_cluster = merged_cluster_index;
    _bonds_in_cluster_with_size_two_or_more += sites.size();

    return merged_cluster_index;
}


/**
 * Take a bond index only if the corresponding site is active
 *
 * @param site
 * @param site_neighbor
 * @param bond_neighbor
 */
void
BondPercolation_pb_v0::connection_v1(BondIndex bond, vector<Index> &site_neighbor, vector<BondIndex> &bond_neighbor) {
    // for a bond there is two site. each site connects to maximum of 3 bonds.
    // select the site if one of the three bonds is active
    value_type next_column = (bond.column_ + 1) % _length;
    value_type prev_column = (bond.column_ + 1 + _length) % _length;
    value_type prev_row = (bond.row_ + 1 + _length) % _length;
    value_type next_row = (bond.row_ + 1) % _length;
    if (!_periodicity) {
        cout << "Testing ... ... .. .. .. . . . : line " << __LINE__ << endl;
        if (bond.horizontal()) {
            site_neighbor.resize(2);
            site_neighbor[0] = {bond.row_, bond.column_}; // left site
            site_neighbor[1] = {bond.row_, next_column}; // right site

            // without periodicity
            if (bond.row_ == min_index) { // top edge including corners
                if (bond.column_ == min_index) {
                    // upper left corner

                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, bond.row_, bond.column_}); // left end of bond
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, bond.row_, next_column}); // right end bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, next_column}); // right end bond
                    }

                } else if (bond.column_ == max_index) {
                    // upper right corner
                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, bond.row_, bond.column_}); // left end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column}); // left end bond
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, bond.row_, next_column}); // right end bond
                    }

                } else {
                    // top edge excluding corners

                    bond_neighbor.reserve(4);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, bond.row_, bond.column_}); // left end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column}); // left end bond
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, bond.row_, next_column}); // right end bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, next_column}); // right end bond
                    }
                }
            } else if (bond.row_ == max_index) { // bottom edge including corners

                if (bond.column_ == min_index) {
                    // lower left corner
                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_}); // left end of bond
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, prev_row, next_column}); // right end bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, next_column}); // right end bond
                    }

                } else if (bond.column_ == max_index) {
                    // lower right corner
                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_}); // left end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column}); // left end of bond
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, prev_row, next_column}); // right end bond
                    }

                } else {
                    // bottom edge excluding corners
                    //  bottom edge
                    bond_neighbor.reserve(4);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_}); // left end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column}); // left end of bond
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        bond_neighbor.push_back({BondType::Vertical, prev_row, next_column}); // right end bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, next_column}); // right end bond
                    }

                }
            }
                /* site.x_ > min_index && site.x_ < max_index &&  is not possible anymore*/
            else if (bond.column_ == min_index) { // left edge excluding corners
                bond_neighbor.reserve(5);
                if (_lattice.getSite(site_neighbor[0]).isActive()) {
                    bond_neighbor.push_back({BondType::Vertical, bond.row_, bond.column_}); // left end of bond
                    bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_}); // left end of bond
                }
                if (_lattice.getSite(site_neighbor[1]).isActive()) {
                    bond_neighbor.push_back({BondType::Vertical, prev_row, next_column}); // right end bond
                    bond_neighbor.push_back({BondType::Vertical, bond.row_, next_column}); // right end bond
                    bond_neighbor.push_back({BondType::Horizontal, bond.row_, next_column}); // right end bond
                }

            } else if (bond.column_ == max_index) {
                // right edge no corners
                bond_neighbor.reserve(5);
                if (_lattice.getSite(site_neighbor[0]).isActive()) {
                    bond_neighbor.push_back({BondType::Vertical, bond.row_, bond.column_}); // left end of bond
                    bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_}); // left end of bond
                    bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column}); // left end of bond
                }
                if (_lattice.getSite(site_neighbor[1]).isActive()) {
                    bond_neighbor.push_back({BondType::Vertical, prev_row, next_column}); // right end bond
                    bond_neighbor.push_back({BondType::Vertical, bond.row_, next_column}); // right end bond
                }

            }

        } else if (bond.vertical()) {
            // without periodicity
            site_neighbor.resize(2);
            site_neighbor[0] = {bond.row_, bond.column_}; // top end of bond
            site_neighbor[1] = {next_row, bond.column_}; // bottom end of bond

            if (bond.row_ == min_index) { // top edge including corners
                if (bond.column_ == min_index) {
                    // upper left corner
                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        // top end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, bond.column_});
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        // bottom end of bond
                        bond_neighbor.push_back({BondType::Horizontal, next_row, bond.column_});
                        bond_neighbor.push_back({BondType::Vertical, next_row, bond.column_});
                    }

                } else if (bond.column_ == max_index) {
                    // upper right corner

                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        // top end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column});
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        // bottom end of bond
                        bond_neighbor.push_back({BondType::Horizontal, next_row, prev_column});
                        bond_neighbor.push_back({BondType::Vertical, next_row, bond.column_});
                    }

                } else {
                    // top edge excluding corners

                    bond_neighbor.reserve(5);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        // top end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, bond.column_});
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column});
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        // bottom end of bond
                        bond_neighbor.push_back({BondType::Horizontal, next_row, prev_column});
                        bond_neighbor.push_back({BondType::Horizontal, next_row, bond.column_});
                        bond_neighbor.push_back({BondType::Vertical, next_row, bond.column_});
                    }

                }
            } else if (bond.row_ == max_index) { // bottom edge including corners
                if (bond.column_ == min_index) {
                    // lower left corner
                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        // top end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, bond.column_});
                        bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_});
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        // bottom end of bond
                        bond_neighbor.push_back({BondType::Horizontal, next_row, bond.column_});

                    }

                } else if (bond.column_ == max_index) {
                    // lower right corner
                    bond_neighbor.reserve(3);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        // top end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column});
                        bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_});
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        // bottom end of bond
                        bond_neighbor.push_back({BondType::Horizontal, next_row, prev_column});
                    }

                } else {
                    // bottom edge excluding corners
                    //  bottom edge
                    bond_neighbor.reserve(5);
                    if (_lattice.getSite(site_neighbor[0]).isActive()) {
                        // top end of bond
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, bond.column_});
                        bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column});
                        bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_});
                    }
                    if (_lattice.getSite(site_neighbor[1]).isActive()) {
                        // bottom end of bond
                        bond_neighbor.push_back({BondType::Horizontal, next_row, bond.column_});
                        bond_neighbor.push_back({BondType::Horizontal, next_row, prev_column});
                    }
                }
            }
                /* site.x_ > min_index && site.x_ < max_index &&  is not possible anymore*/
            else if (bond.column_ == min_index) {
                // left edge not in the corners
                bond_neighbor.reserve(4);
                if (_lattice.getSite(site_neighbor[0]).isActive()) {
                    // top end of bond
                    bond_neighbor.push_back({BondType::Horizontal, bond.row_, bond.column_});
                    bond_neighbor.push_back({BondType::Vertical, bond.row_, bond.column_});
                }
                if (_lattice.getSite(site_neighbor[1]).isActive()) {
                    // bottom end of bond
                    bond_neighbor.push_back({BondType::Horizontal, next_row, bond.column_});
                    bond_neighbor.push_back({BondType::Vertical, next_row, bond.column_});
                }

            } else if (bond.column_ == max_index) {
                // right edge no corners
                bond_neighbor.reserve(4);
                if (_lattice.getSite(site_neighbor[0]).isActive()) {
                    // top end of bond
                    bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column});
                    bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_});
                }
                if (_lattice.getSite(site_neighbor[1]).isActive()) {
                    // bottom end of bond
                    bond_neighbor.push_back({BondType::Vertical, next_row, bond.column_});
                    bond_neighbor.push_back({BondType::Horizontal, next_row, prev_column});
                }
            }
        }
        return;
    }
    // A bond can have 2 site neighbor and 6 bond neighbor
    site_neighbor = vector<Index>(2);
    site_neighbor[0] = {bond.row_, bond.column_};
    // as our BondIndex suggests neighbor[1] should be downward or rightward, we set it below.


    // now we have all the neighboring site
    // it's time to find the 6 neighboring bond_neighbor
    // without periodicity
    if (bond.horizontal()) {
        // increase column index for the right neighbor
        site_neighbor[1] = {bond.row_, next_column};

        bond_neighbor.reserve(6);
        if (_lattice.getSite(site_neighbor[0]).isActive()) {
            // left end of bond
            bond_neighbor.push_back({BondType::Vertical, bond.row_, bond.column_});
            bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_});
            bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column});
        }
        if (_lattice.getSite(site_neighbor[1]).isActive()) {
            // right end bond
            bond_neighbor.push_back({BondType::Vertical, prev_row, next_column});
            bond_neighbor.push_back({BondType::Vertical, bond.row_, next_column});
            bond_neighbor.push_back({BondType::Horizontal, bond.row_, next_column});
        }
    } else if (bond.vertical()) {
        // increase row index
        site_neighbor[1] = {next_row, bond.column_};

        bond_neighbor.reserve(6);
        if (_lattice.getSite(site_neighbor[0]).isActive()) {
            // top end of bond
            bond_neighbor.push_back({BondType::Horizontal, bond.row_, bond.column_});
            bond_neighbor.push_back({BondType::Horizontal, bond.row_, prev_column});
            bond_neighbor.push_back({BondType::Vertical, prev_row, bond.column_});
        }
        if (_lattice.getSite(site_neighbor[1]).isActive()) {
            // bottom end of bond
            bond_neighbor.push_back({BondType::Horizontal, next_row, bond.column_});
            bond_neighbor.push_back({BondType::Horizontal, next_row, prev_column});
            bond_neighbor.push_back({BondType::Vertical, next_row, bond.column_});
        }

    }

}


/**
 * Take a bond index only if the corresponding site is active
 *
 * @param site
 * @param site_neighbor
 * @param bond_neighbor
 */
void
BondPercolation_pb_v0::connection_v2(BondIndex bond, vector<Index> &site_neighbor, vector<BondIndex> &bond_neighbor) {
    // for a bond there is two site. each site connects to maximum of 3 bonds.
    // select the site if one of the three bonds is active
    value_type next_column = (bond.column_ + 1) % _length;
    value_type prev_column = (bond.column_ - 1 + _length) % _length;
    value_type prev_row = (bond.row_ - 1 + _length) % _length;
    value_type next_row = (bond.row_ + 1) % _length;
    if (!_periodicity) {

        cout << "Testing ... ... .. .. .. . . . : line " << __LINE__ << endl;

        if (bond.horizontal()) {
            // without periodicity
            connection_2_horizontal_no_periodicity(bond, site_neighbor, bond_neighbor,
                                                   next_column, prev_column, prev_row);

        } else if (bond.vertical()) {
            // without periodicity
            connection_2_vertical_no_peridicity(bond, site_neighbor, bond_neighbor,
                                                prev_column, prev_row, next_row);
        }

    } else {
        // with periodicity
        connection_2_periodic(bond, site_neighbor, bond_neighbor,
                              prev_column, next_column, prev_row, next_row);

    }
}

/**
 * Can also be used for any bond that is not in the boundary.
 *
 * Take a site index only if the corresponding bond is active. // todo
 *
 * @param bond
 * @param site_neighbor
 * @param bond_neighbor
 * @param next_column
 * @param prev_column
 * @param prev_row
 * @param next_row
 */
void BondPercolation_pb_v0::connection_2_periodic(const BondIndex &bond, vector<Index> &site_neighbor,
                                                  vector<BondIndex> &bond_neighbor, value_type prev_column,
                                                  value_type next_column, value_type prev_row, value_type next_row) {
    // A bond can have 2 site neighbor and 6 bond neighbor
    site_neighbor.reserve(2);

//    vector<BondIndex> bond_tmp(6);
    // as our BondIndex suggests neighbor[1] should be downward or rightward, we set it below.


    // now we have all the neighboring site
    // it's time to find the 6 neighboring bond_neighbor
    // without periodicity
    if (bond.horizontal()) {
        // increase column index for the right neighbor

        bond_neighbor.resize(6);

        // left end of bond
        bond_neighbor[0] = {BondType::Vertical, bond.row_, bond.column_};
        bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_};
        bond_neighbor[2] = {BondType::Horizontal, bond.row_, prev_column};
        if(!(_lattice.getBond(bond_neighbor[0]).isActive() ||
           _lattice.getBond(bond_neighbor[1]).isActive() ||
           _lattice.getBond(bond_neighbor[2]).isActive())
                ){
            site_neighbor.push_back({bond.row_, bond.column_});
        }

        // right end bond
        bond_neighbor[3] = {BondType::Vertical, prev_row, next_column};
        bond_neighbor[4] = {BondType::Vertical, bond.row_, next_column};
        bond_neighbor[5] = {BondType::Horizontal, bond.row_, next_column};
        if(!(_lattice.getBond(bond_neighbor[3]).isActive() ||
           _lattice.getBond(bond_neighbor[4]).isActive() ||
           _lattice.getBond(bond_neighbor[5]).isActive())
                ){
            site_neighbor.push_back({bond.row_, next_column});
        }


    } else if (bond.vertical()) {
        // increase row index

        bond_neighbor.resize(6);

        // top end of bond
        bond_neighbor[0] = {BondType::Horizontal, bond.row_, bond.column_};
        bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column};
        bond_neighbor[2] = {BondType::Vertical, prev_row, bond.column_};
        if(!(_lattice.getBond(bond_neighbor[0]).isActive() ||
             _lattice.getBond(bond_neighbor[1]).isActive() ||
             _lattice.getBond(bond_neighbor[2]).isActive())
                ){
            site_neighbor.push_back({bond.row_, bond.column_});
        }

        // bottom end of bond
        bond_neighbor[3] = {BondType::Horizontal, next_row, bond.column_};
        bond_neighbor[4] = {BondType::Horizontal, next_row, prev_column};
        bond_neighbor[5] = {BondType::Vertical, next_row, bond.column_};
        if(!(_lattice.getBond(bond_neighbor[3]).isActive() ||
             _lattice.getBond(bond_neighbor[4]).isActive() ||
             _lattice.getBond(bond_neighbor[5]).isActive())
                ){
            site_neighbor.push_back({next_row, bond.column_});
        }

    }
}

/**
 * Non-Periodic Connection for vertical bond
 *
 * @param bond : vertical bond
 * @param site_neighbor : neighboring sites
 * @param bond_neighbor : neighboring bonds
 * @param prev_column
 * @param prev_row
 * @param next_row
 */
void BondPercolation_pb_v0::connection_2_vertical_no_peridicity(const BondIndex &bond, vector<Index> &site_neighbor,
                                                                vector<BondIndex> &bond_neighbor,
                                                                value_type prev_column, value_type prev_row,
                                                                value_type next_row) {
    bool top_site{false}, bottom_site{false};

    if (bond.row_ == min_index) { // top edge including corners
        if (bond.column_ == min_index) {
            // upper left corner
            bond_neighbor.resize(3);
            // top end of bond
            bond_neighbor[0] = {BondType::Horizontal, bond.row_, bond.column_};

            top_site = !_lattice.getBond(bond_neighbor[0]).isActive();

            // bottom end of bond
            bond_neighbor[1] = {BondType::Horizontal, next_row, bond.column_};
            bond_neighbor[2] = {BondType::Vertical, next_row, bond.column_};

            bottom_site = !(_lattice.getBond(bond_neighbor[1]).isActive() ||
                          _lattice.getBond(bond_neighbor[2]).isActive());

        } else if (bond.column_ == max_index) {
            // upper right corner

            bond_neighbor.resize(3);
            // top end of bond
            bond_neighbor[0] = {BondType::Horizontal, bond.row_, prev_column};
            top_site = !_lattice.getBond(bond_neighbor[0]).isActive();

            // bottom end of bond
            bond_neighbor[1] = {BondType::Horizontal, next_row, prev_column};
            bond_neighbor[2] = {BondType::Vertical, next_row, bond.column_};

            bottom_site = !(_lattice.getBond(bond_neighbor[1]).isActive() ||
                          _lattice.getBond(bond_neighbor[2]).isActive());
        } else {
            // top edge excluding corners

            bond_neighbor.resize(5);

            // top end of bond
            bond_neighbor[0] = {BondType::Horizontal, bond.row_, bond.column_};
            bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column};

            top_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                       _lattice.getBond(bond_neighbor[1]).isActive());

            // bottom end of bond
            bond_neighbor[2] = {BondType::Horizontal, next_row, prev_column};
            bond_neighbor[3] = {BondType::Horizontal, next_row, bond.column_};
            bond_neighbor[4] = {BondType::Vertical, next_row, bond.column_};

            bottom_site = !(_lattice.getBond(bond_neighbor[2]).isActive() ||
                          _lattice.getBond(bond_neighbor[3]).isActive() ||
                          _lattice.getBond(bond_neighbor[4]).isActive());

        }
    } else if (bond.row_ == max_index) { // bottom edge including corners
        if (bond.column_ == min_index) {
            // lower left corner
            bond_neighbor.resize(3);

            // top end of bond
            bond_neighbor[0] = {BondType::Horizontal, bond.row_, bond.column_};
            bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_};

            top_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                       _lattice.getBond(bond_neighbor[1]).isActive());

            // bottom end of bond
            bond_neighbor[2] = {BondType::Horizontal, next_row, bond.column_};

            bottom_site = !_lattice.getBond(bond_neighbor[2]).isActive();


        } else if (bond.column_ == max_index) {
            // lower right corner
            bond_neighbor.resize(3);

            // top end of bond
            bond_neighbor[0] = {BondType::Horizontal, bond.row_, prev_column};
            bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_};

            top_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                       _lattice.getBond(bond_neighbor[1]).isActive());

            // bottom end of bond
            bond_neighbor[2] = {BondType::Horizontal, next_row, prev_column};

            bottom_site = !_lattice.getBond(bond_neighbor[2]).isActive();

        } else {
            // bottom edge excluding corners
            //  bottom edge
            bond_neighbor.resize(5);

            // top end of bond
            bond_neighbor[0] = {BondType::Horizontal, bond.row_, bond.column_};
            bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column};
            bond_neighbor[2] = {BondType::Vertical, prev_row, bond.column_};

            top_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                       _lattice.getBond(bond_neighbor[1]).isActive() ||
                       _lattice.getBond(bond_neighbor[2]).isActive());

            // bottom end of bond
            bond_neighbor[3] = {BondType::Horizontal, next_row, bond.column_};
            bond_neighbor[4] = {BondType::Horizontal, next_row, prev_column};

            bottom_site = !(_lattice.getBond(bond_neighbor[3]).isActive() ||
                          _lattice.getBond(bond_neighbor[4]).isActive());
        }
    }
        /* site.x_ > min_index && site.x_ < max_index &&  is not possible anymore*/
    else if (bond.column_ == min_index) {
        // left edge not in the corners
        bond_neighbor.resize(4);
        // top end of bond
        bond_neighbor[0] = {BondType::Horizontal, bond.row_, bond.column_};
        bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_};

        top_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                   _lattice.getBond(bond_neighbor[1]).isActive());

        // bottom end of bond
        bond_neighbor[2] = {BondType::Horizontal, next_row, bond.column_};
        bond_neighbor[3] = {BondType::Vertical, next_row, bond.column_};

        bottom_site = !(_lattice.getBond(bond_neighbor[2]).isActive() ||
                      _lattice.getBond(bond_neighbor[3]).isActive());

    } else if (bond.column_ == max_index) {
        // right edge no corners
        bond_neighbor.resize(4);

        // top end of bond
        bond_neighbor[0] = {BondType::Horizontal, bond.row_, prev_column};
        bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_};

        top_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                   _lattice.getBond(bond_neighbor[1]).isActive());

        // bottom end of bond
        bond_neighbor[2] = {BondType::Vertical, next_row, bond.column_};
        bond_neighbor[3] = {BondType::Horizontal, next_row, prev_column};

        bottom_site = !(_lattice.getBond(bond_neighbor[2]).isActive() ||
                      _lattice.getBond(bond_neighbor[3]).isActive());
    }else{
        // one level inside the lattice
        bond_neighbor.resize(6);

        // top end of bond
        bond_neighbor[0] = {BondType::Horizontal, bond.row_, bond.column_};
        bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column};
        bond_neighbor[2] = {BondType::Vertical, prev_row, bond.column_};

        top_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                   _lattice.getBond(bond_neighbor[1]).isActive() ||
                   _lattice.getBond(bond_neighbor[2]).isActive());

        // bottom end of bond
        bond_neighbor[3] = {BondType::Horizontal, next_row, bond.column_};
        bond_neighbor[4] = {BondType::Horizontal, next_row, prev_column};
        bond_neighbor[5] = {BondType::Vertical, next_row, bond.column_};

        bottom_site = !(_lattice.getBond(bond_neighbor[3]).isActive() ||
                      _lattice.getBond(bond_neighbor[4]).isActive() ||
                      _lattice.getBond(bond_neighbor[5]).isActive());
    }

    site_neighbor.reserve(2);
    if (top_site) {
        site_neighbor.push_back({bond.row_, bond.column_});
    }
    if (bottom_site) {
        site_neighbor.push_back({next_row, bond.column_});
    }
}

void BondPercolation_pb_v0::connection_2_horizontal_no_periodicity(const BondIndex &bond,
                                                                   vector<Index> &site_neighbor,
                                                                   vector<BondIndex> &bond_neighbor,
                                                                   value_type next_column,
                                                                   value_type prev_column,
                                                                   value_type prev_row) {
    bool left_site{false}, right_site{false};

    // without periodicity
    if (bond.row_ == min_index) {
        // top edge including corners
        if (bond.column_ == min_index) {
            // upper left corner

            bond_neighbor.resize(3);
            // left end of bond
            bond_neighbor[0] = {BondType::Vertical, bond.row_, bond.column_};
            // right end bond
            bond_neighbor[1] = {BondType::Vertical, bond.row_, next_column};
            bond_neighbor[2] = {BondType::Horizontal, bond.row_, next_column};

            left_site = !_lattice.getBond(bond_neighbor[0]).isActive();

            right_site = !(_lattice.getBond(bond_neighbor[1]).isActive() ||
                         _lattice.getBond(bond_neighbor[2]).isActive());

        } else if (bond.column_ == max_index) {
            // upper right corner
            bond_neighbor.resize(3);
            // left end of bond
            bond_neighbor[0] = {BondType::Vertical, bond.row_, bond.column_};
            bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column};

            // right end bond
            bond_neighbor[3] = {BondType::Vertical, bond.row_, next_column};

            left_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                        _lattice.getBond(bond_neighbor[1]).isActive());

            right_site = !_lattice.getBond(bond_neighbor[3]).isActive();

        } else {
            // top edge excluding corners

            bond_neighbor.resize(4);
            // left end of bond
            bond_neighbor[0] = {BondType::Vertical, bond.row_, bond.column_};
            bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column}; // left end bond

            left_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                        _lattice.getBond(bond_neighbor[1]).isActive());

            // right end bond
            bond_neighbor[2] = {BondType::Vertical, bond.row_, next_column};
            bond_neighbor[3] = {BondType::Horizontal, bond.row_, next_column}; // right end bond

            right_site = !(_lattice.getBond(bond_neighbor[2]).isActive() ||
                         _lattice.getBond(bond_neighbor[3]).isActive());
        }
    } else if (bond.row_ == max_index) {
        // bottom edge including corners

        if (bond.column_ == min_index) {
            // lower left corner
            bond_neighbor.resize(3);
            // left end of bond
            bond_neighbor[0] = {BondType::Vertical, prev_row, bond.column_};
            left_site = !_lattice.getBond(bond_neighbor[0]).isActive();

            // right end bond
            bond_neighbor[1] = {BondType::Vertical, prev_row, next_column}; // right end bond
            bond_neighbor[2] = {BondType::Horizontal, bond.row_, next_column};

            right_site = !(_lattice.getBond(bond_neighbor[1]).isActive() ||
                         _lattice.getBond(bond_neighbor[2]).isActive());

        } else if (bond.column_ == max_index) {
            // lower right corner
            bond_neighbor.resize(3);
            // left end of bond
            bond_neighbor[0] = {BondType::Vertical, prev_row, bond.column_};
            bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column}; // left end of bond

            left_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                        _lattice.getBond(bond_neighbor[1]).isActive());

            // right end bond
            bond_neighbor[2] = {BondType::Vertical, prev_row, next_column};

            right_site = !_lattice.getBond(bond_neighbor[2]).isActive();

        } else {
            // bottom edge excluding corners
            //  bottom edge
            bond_neighbor.resize(4);

            bond_neighbor[0] = {BondType::Vertical, prev_row, bond.column_}; // left end of bond
            bond_neighbor[1] = {BondType::Horizontal, bond.row_, prev_column}; // left end of bond

            left_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                        _lattice.getBond(bond_neighbor[1]).isActive());

            bond_neighbor[2] = {BondType::Vertical, prev_row, next_column}; // right end bond
            bond_neighbor[3] = {BondType::Horizontal, bond.row_, next_column}; // right end bond

            right_site = !(_lattice.getBond(bond_neighbor[2]).isActive() ||
                         _lattice.getBond(bond_neighbor[3]).isActive());


        }
    }
        /* site.x_ > min_index && site.x_ < max_index &&  is not possible anymore*/
    else if (bond.column_ == min_index) { // left edge excluding corners
        bond_neighbor.resize(5);

        bond_neighbor[0] = {BondType::Vertical, bond.row_, bond.column_}; // left end of bond
        bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_}; // left end of bond

        left_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                    _lattice.getBond(bond_neighbor[1]).isActive());

        bond_neighbor[2] = {BondType::Vertical, prev_row, next_column}; // right end bond
        bond_neighbor[3] = {BondType::Vertical, bond.row_, next_column}; // right end bond
        bond_neighbor[4] = {BondType::Horizontal, bond.row_, next_column}; // right end bond

        right_site = !(_lattice.getBond(bond_neighbor[2]).isActive() ||
                     _lattice.getBond(bond_neighbor[3]).isActive() ||
                     _lattice.getBond(bond_neighbor[4]).isActive());

    } else if (bond.column_ == max_index) {
        // right edge no corners
        bond_neighbor.resize(5);

        bond_neighbor[0] = {BondType::Vertical, bond.row_, bond.column_}; // left end of bond
        bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_}; // left end of bond
        bond_neighbor[2] = {BondType::Horizontal, bond.row_, prev_column}; // left end of bond

        left_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                    _lattice.getBond(bond_neighbor[1]).isActive() ||
                    _lattice.getBond(bond_neighbor[2]).isActive());

        bond_neighbor[3] = {BondType::Vertical, prev_row, next_column}; // right end bond
        bond_neighbor[4] = {BondType::Vertical, bond.row_, next_column}; // right end bond

        right_site = !(_lattice.getBond(bond_neighbor[3]).isActive() ||
                     _lattice.getBond(bond_neighbor[4]).isActive());

    }else{
        // one level inside the lattice

        bond_neighbor.resize(6);

        // left end of bond
        bond_neighbor[0] = {BondType::Vertical, bond.row_, bond.column_};
        bond_neighbor[1] = {BondType::Vertical, prev_row, bond.column_};
        bond_neighbor[2] = {BondType::Horizontal, bond.row_, prev_column};

        left_site = !(_lattice.getBond(bond_neighbor[0]).isActive() ||
                    _lattice.getBond(bond_neighbor[1]).isActive() ||
                    _lattice.getBond(bond_neighbor[2]).isActive());

        // right end bond
        bond_neighbor[3] = {BondType::Vertical, prev_row, next_column};
        bond_neighbor[4] = {BondType::Vertical, bond.row_, next_column};
        bond_neighbor[5] = {BondType::Horizontal, bond.row_, next_column};

        right_site = !(_lattice.getBond(bond_neighbor[3]).isActive() ||
                     _lattice.getBond(bond_neighbor[4]).isActive() ||
                     _lattice.getBond(bond_neighbor[5]).isActive());
    }

    site_neighbor.reserve(2);
    if (left_site) {
        site_neighbor.push_back({bond.row_, bond.column_});
    }
    if (right_site) {
        site_neighbor.push_back({bond.row_, next_column});
    }
}

/**
 * Takes 0.093 sec time when total time is 9.562 sec
 *
 * Find one row from _cluster to place 4 or less new bonds
 * Also remove the matched index values, because they will be inserted later.
 * This gives an advantage, i.e., you don't need to perform a checking.
 * todo takes so much time
 * @param hv_bonds
 * @return a set
 */
set<value_type>
BondPercolation_pb_v0::find_index_for_placing_new_bonds_v1(const vector<BondIndex> &neighbors) {

    value_type total_found_bonds{}; // for checking program running Time only

    set<value_type> found_index_set;    // use set to prevent repeated index
    value_type x{};

    if (!_cluster_index_from_id.empty()) {
        for (auto n: neighbors) {
            if (_cluster_index_from_id.count(_lattice.getBond(n).groupID()) > 0) {
                x = _cluster_index_from_id[_lattice.getBond(n).groupID()];
                found_index_set.insert(x);

            }
        }
    }

    return found_index_set;
}


void BondPercolation_pb_v0::calculate_spanning_probability() {
    calculate_spanning_probability_by_largest_cluster();
}


/**
 * Use Group_ID to identify Bond and Site in the same cluster
 */
void BondPercolation_pb_v0::calculate_spanning_probability_by_largest_cluster() {
    if (debug_calculate_spanning_probability_by_largest_cluster_3) {
        cout << "calculate_spanning_probability_by_largest_cluster() : line " << __LINE__ << endl;
    }
    // find the largest cluster
    value_type l_largest_cluster{0}; // number of site of the largest cluster
    double b;
    for (auto a: _clusters) {
        b = a.numberOfSites();
        if (b < _length) // no percolation yet
            continue;
        if (b > l_largest_cluster) {
            l_largest_cluster = b;
        }
    }
//    double spanning_prob = l_largest_cluster/_number_of_occupied_sites;

    number_of_bonds_to_span.push_back(l_largest_cluster);

}


/**
 * Counts the number of active sites in the lattice
 */
value_type BondPercolation_pb_v0::count_number_of_active_site() {
    value_type counter{};
    for (value_type i{}; i != _length; ++i) {
        for (value_type j{}; j != _length; ++j) {
            if (_lattice.getSite({i, j}).isActive())
                ++counter;
        }
    }
    cout << "Number of active site : " << counter << endl;
    return counter;
}


//
///**
// * Give a site Index and it will give 4 bonds corresponding to that site
// * @param index
// * @return
// */
//vector<Index> BondPercolation_pb_v0::get_Sites_for_bonds(BondIndex index) {
//    if (debug_get_Bond_for_site) {
//        cout << "Entry -> get_Bond_for_site() : line " << __LINE__ << endl;
//    }
//    vector<Index> connected_sites(2);
//
//    if(index.horizontal()) {
//        // for horizontal bond, row remains the same
//        connected_sites[0] = _lattice.getSite({index.x_, index.column_}).ID();
//        auto c = (index.column_ + 1) % _length;
//        connected_sites[1] = _lattice.getSite({index.x_, c}).ID();
//    }else{
//        // for vertical bond, column remains the same
//        connected_sites[0] = _lattice.getSite({index.x_, index.column_}).ID();
//        auto r = (index.x_ + 1) % _length;
//        connected_sites[1] = _lattice.getSite({r, index.column_}).ID();
//    }
//
//
//    if (debug_get_Bond_for_site) {
//        cout << "Bonds for untouched site " << index << endl;
//        for (auto it = connected_sites.begin(); it != connected_sites.end(); ++it) {
//            cout << "ID " << (*it) << " = " << *it << endl;
//        }
//    }
//    return connected_sites;
//}



value_type BondPercolation_pb_v0::number_of_site_in_spanning_clusters(unordered_set<int> g_ids) {
    value_type nos{};   // number of sites
    for (auto b: g_ids) {
        for (auto a: _clusters) {
            if (a.ID() == b) {
                nos += a.numberOfSites();
            }
        }
    }
    return nos;
}


/***
 *
 * @return
 */
bool BondPercolation_pb_v0::detectSpanning() {
    if (debug_4_detectSpanning) {
        cout << "Entry -> detectSpanning() : line " << __LINE__ << endl;
    }

    // if any of the clusters does not have sites > _length --> no spanning
    value_type x{};
    for (auto a: _clusters) {
        if (x < a.numberOfSites())
            x = a.numberOfSites();
    }
    if (x < _length) {
        if (debug_4_detectSpanning) {
            cout << "not enough sites in any cluster to span : line " << __LINE__ << endl;
        }
        number_of_bonds_to_span.push_back(0);
        return false;
    }

    // vertical spanning
    std::unordered_set<int> first_row_group_ids;
    std::unordered_set<int> last_row_group_ids;
    // horizontal spanning
    std::unordered_set<int> first_column_group_ids;
    std::unordered_set<int> last_column_group_ids;

    for (value_type i{}; i != _length; ++i) {
        // for rows
        if (_lattice.getSite({0, i}).isActive()) {
            first_row_group_ids.insert(_lattice.getSite({0, i}).groupID());
        }
        if (_lattice.getSite({_length - 1, i}).isActive()) {
            last_row_group_ids.insert(_lattice.getSite({_length - 1, i}).groupID());
        }

        // for columns
        if (_lattice.getSite({i, 0}).isActive()) {
            first_column_group_ids.insert(_lattice.getSite({i, 0}).groupID());
        }
        if (_lattice.getSite({i, _length - 1}).isActive()) {
            last_column_group_ids.insert(_lattice.getSite({i, _length - 1}).groupID());
        }

    }

//    cout << '{';
//    for(auto a: first_row_group_ids){
//        cout << a << ',';
//    }
//    cout << '}' << endl;
//
//    cout << '{';
//    for(auto b: last_row_group_ids){
//        cout << b << ',';
//    }
//    cout << '}' << endl;

    std::vector<int> vertical_spanning_id;
    std::vector<int> horizontal_spanning_id;




    // to prevent repetitionn, use set
    unordered_set<int> spanning_ids;
    for (auto a: first_row_group_ids) {
        for (auto b: last_row_group_ids) {
            if (a == b) {
                spanning_ids.insert(a);
                vertical_spanning_id.push_back(a);
            }
        }
    }

    for (auto a: first_column_group_ids) {
        for (auto b: last_column_group_ids) {
            if (a == b) {
                spanning_ids.insert(a);
                horizontal_spanning_id.push_back(a);
            }
        }
    }


    value_type number_of_sites_in_spanning_cluster = number_of_site_in_spanning_clusters(spanning_ids);

    number_of_bonds_to_span.push_back(number_of_sites_in_spanning_cluster);

    if (debug_4_detectSpanning) {
        cout << "Number of sites to span " << number_of_sites_in_spanning_cluster << ','
             << "Number of active sites " << _total_number_of_active_bonds << ','
             << "Total number of sites " << _length_squared << endl;
        cout << "Number of cluters " << _clusters.size() << endl;
        cout << "Spanning Info : line " << __LINE__ << endl;
        cout << "Vertical spanning for cluster_ID   : " << vertical_spanning_id << endl;
        cout << "Horizontal spanning for cluster_ID : " << horizontal_spanning_id << endl;
        viewClusterExtended();
        viewLatticeExtended();
    }

    return !spanning_ids.empty();
}


std::vector<double> BondPercolation_pb_v0::spanningProbability() const {
    vector<double> x(number_of_bonds_to_span.size());
    for (value_type i{}; i != x.size(); ++i) {
        x[i] = number_of_bonds_to_span[i] / double(_length_squared);
    }
    return x;
}


double BondPercolation_pb_v0::occupationProbability() const {
    return _number_of_occupied_bonds / double(_length_squared);
}


/**
 *
 */
void BondPercolation_pb_v0::periodicity_status() {
    cout << "Periodicity = (" << _periodicity << ")" << endl;
}





