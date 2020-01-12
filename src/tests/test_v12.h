//
// Created by shahnoor on 12/23/19.
//

#ifndef SQLATTICEPERCOLATION_TEST_V12_H
#define SQLATTICEPERCOLATION_TEST_V12_H

/*
 * Date : 2019.12.23
 * Author : Shahnoor
 */
void test_v12(int argc, char *argv[]);

void test_v12_lattice(int length);
void test_v12_percolation_site(int length);
void test_v12_percolation_bond(int length);
void run_v12_regular_site(int length, int ensemble_size) ;
void run_v12_regular_bond(int length, int ensemble_size) ;

void demarcationLine_v12(int length);

#endif //SQLATTICEPERCOLATION_TEST_V12_H
