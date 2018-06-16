//
// Created by shahnoor on 12/4/2017.
//

#ifndef PERCOLATION_PERCOLATION_FRACTAL_H
#define PERCOLATION_PERCOLATION_FRACTAL_H

#include "../types.h"
#include <iostream>
#include <fstream>
#include "../percolation/percolation.h"
#include "../util/time_tracking.h"


void moment_data(value_type  L_start, value_type  L_end, value_type ensembleSize=100, double df=1.789);
double moment_for_percolation_by_site(value_type ensemble_size, value_type length, double d_f=1.789);

void spanning_data_for_lengths(value_type L_start, value_type L_end, value_type ensembleSize,
                               value_type impure_site = 0);

void spanning_data(value_type ensemble_size, value_type length, value_type impure_site = 0);


#endif //PERCOLATION_PERCOLATION_FRACTAL_H
