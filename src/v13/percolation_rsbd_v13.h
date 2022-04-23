//
// Created by shahnoor on 2/15/21.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_RSBD_H
#define SQLATTICEPERCOLATION_PERCOLATION_RSBD_H


#include "lattice_v13.h"
#include "cluster.h"
#include "percolation_v13.h"
#include "../types.h"
#include "../../../../../../usr/include/c++/11/string"
#include <set>
#include <algorithm>
#include <random>



class SitePercolationL1_v13: public SitePercolationL0_v13{
    // int x_occupied = 0;

    // most of the times we select the originally chosen random site but sometimes that 
    // site is not empty and we select 1st nearest neighbor. The ratio of 
    // 1st nn to total site in a particular lattice might be an interesting quantity.
    // also the trend of doing it.
    int nn_1st_count = 0;
    int original_count = 0;
public:

    SitePercolationL1_v13(int length, value_type seed=0, bool generate_seed=true);
    SitePercolationL1_v13(SitePercolationL1_v13&&) = default;
    std::string get_signature() override { return SitePercolation_v13::get_signature() + "L1_";};

    void reset() {
        SitePercolationL0_v13::reset(); 
        // x_occupied=0;
        nn_1st_count = 0;
        original_count = 0;
    }


    P_STATUS select_site() override;

    double get_1st_nn_to_total_site_ratio(){return double(nn_1st_count) / (original_count + nn_1st_count);}

//     void run_once_v2() {
//         double max_ratio = 0;
//         // SitePercolationL0_v13::run_once_v2();
//         //# sq_lattice_p.viewLattice(3)
// //# sq_lattice_p.viewCluster()
//     double p, H, P1, P2;

//     while (place_one_site()) {
//         if (status != P_STATUS::SUCESS) continue;
//         detect_wrapping();
//         if (first_run) {
//             p = occupation_prob();
//             occupation_prob_list.push_back(p);
//         }else{
//             double tmp = get_1st_nn_to_total_site_ratio();
//             if (tmp > max_ratio) max_ratio = tmp;
//             // std::cout << "a/(a+b) : " << tmp << std::endl;
//         }
//         H = entropy();
//         P1 = order_param_wrapping();
//         P2 = order_param_largest_clstr();

//         entropy_list.push_back(H);
//         order_wrapping_list.push_back(P1);
//         order_largest_list.push_back(P2);

//     }

//     first_run = false;

//         std::cout << "a/(a+b) : " << get_1st_nn_to_total_site_ratio() << " . max ratio " << max_ratio << std::endl;
//     }
};



class SitePercolationL2_v13: public SitePercolationL0_v13{

public:

    SitePercolationL2_v13(int length, value_type seed=0, bool generate_seed=true);
    SitePercolationL2_v13(SitePercolationL2_v13&&) = default;
    std::string get_signature() override { return SitePercolation_v13::get_signature() + "L2_";};

    int get_2nn_in_1nn_direction(int central_X_id, int nn_1_id);

    P_STATUS select_site();

};



#endif //SQLATTICEPERCOLATION_PERCOLATION_RSBD_H
