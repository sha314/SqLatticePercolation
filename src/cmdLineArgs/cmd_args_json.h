//
// Created by Shahnoor on 6/7/2018.
//

#ifndef PERCOLATION_CMD_ARGS_JSON_H
#define PERCOLATION_CMD_ARGS_JSON_H

/**
 * Write in JSON format.
 * good for parsing in python.
 */

#include "../types.h"

#include <vector>
#include <string>

// defined on 07 June 2018
void percolation_sq_lattice_periodic(value_type length, value_type ensembleSize);
void percolation_sq_lattice_non_periodic(value_type length, value_type ensembleSize);

void write_percolation_data_json(
        value_type length,
        const std::vector<double> &nob1,
        const std::vector<double> &nob2,
        const std::vector<double> &entrpy,
        const std::string &filename,
        const std::string &header_info
);

void write_critical_data_json(
        std::vector<double> pcs,
        std::string filename,
        const std::string &header_info
);


void print_help_json();
void cmd_args_json(int argc, char** argv);


// template fuctions
template <class PERCOLATION>
void percolation_sq_lattice_periodic_g(value_type length, value_type ensembleSize);
template <class PERCOLATION>
void percolation_sq_lattice_non_periodic_g(value_type length, value_type ensembleSize);

template <class PERCOLATION>
void simulate(int p, value_type length, value_type ensembleSize);

void cmd_args_json_g(int argc, char **argv);


#endif //PERCOLATION_CMD_ARGS_JSON_H
