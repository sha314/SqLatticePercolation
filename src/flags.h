//
// Created by shahnoor on 10/2/2017.
//

#ifndef SITEPERCOLATION_FLAGS_H
#define SITEPERCOLATION_FLAGS_H

// debug flags


const bool debug_connect                    = false;
const bool debug_placeSite                  = false;
const bool debug_check_all_cluster          = false;
const bool debug_check_cluster              = false;
const bool debug_check_cluster_v2           = false;
const bool debug_get_Bond_for_site          = false;
const bool debug_activate_site              = false;


// birthTime flags
const bool time_check_all_cluster     = false;
const bool time_check_cluster         = false;
const bool time_check_cluster_v2      = false;
const bool time_find_index_for_placing_new_bonds = false;


// version 3
const bool debug_3_placeSiteWithRelabeling = false;
const bool debug_calculate_spanning_probability_by_largest_cluster_3 = false;



// version 4
const bool debug_4_placeSiteWithRelabeling = false;
const bool debug_4_detectSpanning = false;
const bool debug_4_get_Bond_for_site_no_periodicity = false;
const bool debug_4_get_Bond_for_site_no_periodicity_2 = false;
const bool debug_4_configure = true;
const bool debug_4_bond_set_index = false;
const bool debug_time_4_relabel_sites_old = false;
const bool debug_time_4_relabel_sites = false;
const bool debug_time_4_put_values_to_the_cluster = false;
const bool debug_time_4_put_values_to_the_cluster_weighted_relabeling = false;


// version 5
const bool debug_5_placeSiteWithRelabeling = false;
const bool debug_5_detectSpanning = false;
const bool debug_5_get_Bond_for_site_no_periodicity_2 = false;
const bool debug_5_configure = true;
const bool debug_5_bond_set_index = false;
const bool debug_time_5_relabel_sites_old = false;
const bool debug_time_5_relabel_sites = false;
const bool debug_time_5_put_values_to_the_cluster = true;
const bool debug_time_5_put_values_to_the_cluster_weighted_relabeling = true;
const bool debug_5_placeSiteWeightedRelabeling = false;
const bool debug_5_calculate_spanning_probability_by_largest_cluster = false;
const bool debug_5_get_Bond_for_site = false;
const bool debug_time_find_index_for_placing_new_bonds_v2 = false;
const bool debug_5_time_relabelMap = false;

#endif //SITEPERCOLATION_FLAGS_H
