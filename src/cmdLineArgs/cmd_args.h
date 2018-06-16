//
// Created by shahnoor on 9/30/2017.
//

#ifndef PERCOLATION_2_CMD_ARGS_H
#define PERCOLATION_2_CMD_ARGS_H

#include "../types.h"

void print_help();
void print_help_ballistic(); // defined on 02 Apr 2018

void cmd_args_site_percolation(int argc, char **argv);
void cmd_args_site_percolation_ballistic_deposition(int argc, char **argv);


void cmd_moment_conservation(int argc, char** argv);
void cmd_cluster_data(int argc, char** argv);
void percolation_data_canonical(value_type length, value_type ensembleSize);


void critical_point_data_sq_lattice_site_percolation(value_type length, value_type ensemble_size,
                                                     value_type impure_sites = 0);
void critical_point_data_explosive(value_type length, value_type ensemble_size, int rule, value_type impure_sites=0);
// defined on 02 Apr 2018
void critical_point_data_sq_l_sp_ballistic_deposition_L1(value_type length, value_type ensemble_size);


// defined on 26 Jan 2018
void entropy_order_parameter(value_type length, value_type ensembleSize, value_type impure_site=0);
void entropy_order_parameter_explosive(value_type length, value_type ensembleSize, int rule, value_type impure_site=0);


// defined on 31 Jan 2018
void calculate_and_write_to_file(value_type length, value_type ensembleSize, value_type total_number_of_bonds,
                                 const std::vector<double> &nob1,
                                 const std::vector<double> &entrpy, const std::string &filename);

// defined on 31 Jan 2018
void percolation_data_sq_lattice(value_type length, value_type ensembleSize);
void percolation_data_sq_lattice_explosive(value_type length, value_type ensembleSize, int rule);

// defined on 02 Apr 2018
void percolation_data_sq_lattice_ballistic(value_type length, value_type ensembleSize);
// defined on 22 Arp 2018
void percolation_data_sq_lattice_ballistic_non_periodic(value_type length, value_type ensembleSize);


// defined on 04 Feb, 2018
void box_counting_square_lattice(value_type length, value_type ensemble_size, value_type impure_site);
void cluster_length_sq_lattice(value_type length, value_type ensemble_size, value_type impure_sites);


#endif //PERCOLATION_2_CMD_ARGS_H
