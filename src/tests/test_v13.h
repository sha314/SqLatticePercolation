//
// Created by shahnoor on 2/15/21.
//

#ifndef SQLATTICEPERCOLATION_TEST_V13_H
#define SQLATTICEPERCOLATION_TEST_V13_H

#include "../types.h"

void test_lattice();
void test_detect_wrapping();

void test_v13(int argc, char** argv);
void run_ensemble_v13(int argc, char** argv);

void percolation_fractalDimension_by_spanning_site_v13(value_type ensemble_size, value_type L_start, value_type l_end, value_type delta_L);
#endif //SQLATTICEPERCOLATION_TEST_V13_H
