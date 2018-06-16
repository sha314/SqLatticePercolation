//
// Created by shahnoor on 04-May-18.
//

#ifndef PERCOLATION_SITE_PERCOLATION_H
#define PERCOLATION_SITE_PERCOLATION_H

#include <cstdlib>

void calculate_and_write_to_file_v2(value_type length, value_type ensembleSize, value_type total_number_of_bonds,
                                 const std::vector<double> &nob1, const std::vector<double> &nob2,
                                 const std::vector<double> &entrpy, const std::string &filename);

// todo make a template that will do all
template <class TPercolation>
void sq_lattice_site_percolation_periodic(size_t length, size_t ensemble_size);


void sq_lattice_site_percolation_periodic(size_t length, size_t ensemble_size);
void sq_lattice_site_percolation_non_periodic(size_t length, size_t ensemble_size);


void sq_lattice_site_percolation_explosive_periodic(size_t length, size_t ensemble_size);
void sq_lattice_site_percolation_explosive_non_periodic(size_t length, size_t ensemble_size);


/*********************
 * Bond Percolation
 */
void sq_lattice_bond_percolation_periodic(size_t length, size_t ensemble_size);
void sq_lattice_bond_percolation_non_periodic(size_t length, size_t ensemble_size);


/***********
 * Ballistic deposition in site percolation in square lattice
 */
void sq_lattice_site_percolation_ballistic_deposition_L1_periodic(size_t length, size_t ensemble_size);
void sq_lattice_site_percolation_ballistic_deposition_L1_non_periodic(size_t length, size_t ensemble_size);

void sq_lattice_site_percolation_ballistic_deposition_L2_periodic(size_t length, size_t ensemble_size);
void sq_lattice_site_percolation_ballistic_deposition_L2_non_periodic(size_t length, size_t ensemble_size);


#endif //PERCOLATION_SITE_PERCOLATION_H
