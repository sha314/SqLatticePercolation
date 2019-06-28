//
// Created by shahnoor on 6/28/19.
//

#ifndef SQLATTICEPERCOLATION_PERCOLATION_SITE_EXPLOSIVE_H
#define SQLATTICEPERCOLATION_PERCOLATION_SITE_EXPLOSIVE_H

#include "../percolation.h"
#include "percolation_site_v10.h"


/******************************************************************************
 * Explosive site percolation in square lattice with sum rule and product rule
 */
class SitePercolationExplosive: public SitePercolation_ps_v9{

public:
    ~SitePercolationExplosive() = default;
    SitePercolationExplosive(value_type length);
    std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_explosive_";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        s += std::to_string(length());
        return s;
    }

};

/*******************************************************************************
 * Site Percolation Ballistic Deposition
 *
 * *************************************************************/
/****************************************
 *
 * Extended from SitePercolation_ps_v9
 * definition:
 *
 *
 * specification:
 *
 *
 * performence:
 *      1. takes power law iteration with length with exponent ~2.2 in debug mode
 *          iteration_needed = constant*Length ^ (2.2)
 *      2. time ?? todo
 */
class SitePercolationBallisticDeposition_v2: public SitePercolation_ps_v9{
protected:
    // elements of @indices_tmp will be erased if needed but not of @indices
    std::vector<value_type> indices;
    std::vector<value_type> indices_tmp;
    value_type _search_position{};
public:
    static constexpr const char* signature = "SitePercolation_BallisticDeposition_v2";
    virtual ~SitePercolationBallisticDeposition_v2(){
        indices.clear();
        indices_tmp.clear();
    };
    SitePercolationBallisticDeposition_v2(value_type length, bool periodicity);

    virtual bool occupy();

    /************************************
     * Site selection methods
     */
    Index select_site(std::vector<Index> &sites, std::vector<BondIndex> &bonds);
    Index select_site_upto_1nn(std::vector<Index> &sites, std::vector<BondIndex> &bonds);
    Index select_site_upto_2nn(std::vector<Index> &sites, std::vector<BondIndex> &bonds);

    Index select_site_upto_1nn_v2(std::vector<Index> &sites, std::vector<BondIndex> &bonds);
    Index select_site_upto_2nn_v2(std::vector<Index> &sites, std::vector<BondIndex> &bonds);


    void reset(); // todo
    void initialize_indices();
//    void randomize_index();

    virtual std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_ballistic_deposition_";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        return s;
    }


    /***********************************
     * occupy upto 1st nearset neighbor.
     * If the randomly selected site is occupied then select one of the nearest neighor randomly
     * If it is also occupied skip the rest setps and start next iteration Else occupy it
     *
     *
     */

    value_type placeSite_1nn_v2();


    /*********************************
     * occupy upto 2nd nearest neighbor.
     * If the randomly selected site is occupied then select one of the nearest neighor randomly
     * If it is also occupied, select the next neighbor in the direction of motion Else occupy it.
     * If the 2nd nearest neighbor in the direction of motion is also occupied then skip the rest of the steps
     *      and start the next iteration
     *
     */

    value_type placeSite_2nn_v1();

};

/***********
 * Only L1
 */
class SitePercolationBallisticDeposition_L1_v2: public SitePercolationBallisticDeposition_v2{
public:
    ~SitePercolationBallisticDeposition_L1_v2() = default;
    SitePercolationBallisticDeposition_L1_v2(value_type length, bool periodicity)
            :SitePercolationBallisticDeposition_v2(length, periodicity){}

    bool occupy() {
        // if no site is available then return false

        if(_number_of_occupied_sites == maxSites()){
            return false;
        }

        try {
//        value_type v = placeSite_1nn_v0(); // debugging version
            value_type v = placeSite_1nn_v2();
            _occuption_probability = occupationProbability(); // for super class
            return v != ULLONG_MAX;
        }catch (OccupiedNeighbor& on){
//        on.what();
//        cout << "line : " << __LINE__ << endl;
            return false;
        }

    }

    std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_ballistic_deposition_L1";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        return s;
    }

};

/*********************
 *
 */
class SitePercolationBallisticDeposition_L2_v2: public SitePercolationBallisticDeposition_v2{
public:
    ~SitePercolationBallisticDeposition_L2_v2() = default;
    SitePercolationBallisticDeposition_L2_v2(value_type length, bool periodicity)
            :SitePercolationBallisticDeposition_v2(length, periodicity){}

    bool occupy() {
        // if no site is available then return false

        if(_number_of_occupied_sites == maxSites()){
            return false;
        }

        try {

//            value_type v = placeSite_2nn_v0();
            value_type v = placeSite_2nn_v1();
            _occuption_probability = occupationProbability(); // for super class

            return v != ULLONG_MAX;
        }catch (OccupiedNeighbor& on){
//        on.what();
//        cout << "line : " << __LINE__ << endl;
            return false;
        }

    }

    std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_ballistic_deposition_L2";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        return s;
    }

};



/************************************
 * Extended from SitePercolation_ps_v10
 *
 * definition:
 *
 *
 * specification:
 *
 *
 * performence:
 *      1. takes power law iteration with length with exponent ~2.2 in debug mode
 *          iteration_needed = constant*Length ^ (2.2)
 *      2. time ?? todo
 */
class SitePercolationBallisticDeposition_v3: public SitePercolation_ps_v10{
protected:
    // elements of @indices_tmp will be erased if needed but not of @indices
    std::vector<value_type> indices;
    std::vector<value_type> indices_tmp;
    value_type _search_position{};
public:
    static constexpr const char* signature = "SitePercolation_BallisticDeposition_v2";
    virtual ~SitePercolationBallisticDeposition_v3(){
        indices.clear();
        indices_tmp.clear();
    };
    SitePercolationBallisticDeposition_v3(value_type length, bool periodicity);

    virtual bool occupy();

    /************************************
     * Site selection methods
     */
    Index select_site(std::vector<Index> &sites, std::vector<BondIndex> &bonds);

    Index select_site_upto_1nn_v2(std::vector<Index> &sites, std::vector<BondIndex> &bonds);
    Index select_site_upto_2nn_v2(std::vector<Index> &sites, std::vector<BondIndex> &bonds);


    void reset(); // todo
    void initialize_indices();
//    void randomize_index();

    virtual std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_ballistic_deposition_";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        return s;
    }


    /***********************************
     * occupy upto 1st nearset neighbor.
     * If the randomly selected site is occupied then select one of the nearest neighor randomly
     * If it is also occupied skip the rest setps and start next iteration Else occupy it
     *
     *
     */

    value_type placeSite_1nn_v2();


    /*********************************
     * occupy upto 2nd nearest neighbor.
     * If the randomly selected site is occupied then select one of the nearest neighor randomly
     * If it is also occupied, select the next neighbor in the direction of motion Else occupy it.
     * If the 2nd nearest neighbor in the direction of motion is also occupied then skip the rest of the steps
     *      and start the next iteration
     *
     */

    value_type placeSite_2nn_v1();

};

/***********
 * Only L1
 */
class SitePercolationBallisticDeposition_L1_v3: public SitePercolationBallisticDeposition_v3{
public:
    ~SitePercolationBallisticDeposition_L1_v3() = default;
    SitePercolationBallisticDeposition_L1_v3(value_type length, bool periodicity)
            :SitePercolationBallisticDeposition_v3(length, periodicity){}

    bool occupy() {
        // if no site is available then return false

        if(_number_of_occupied_sites == maxSites()){
            return false;
        }

        try {
//        value_type v = placeSite_1nn_v0(); // debugging version
            value_type v = placeSite_1nn_v2();
            _occuption_probability = occupationProbability(); // for super class
            return v != ULLONG_MAX;
        }catch (OccupiedNeighbor& on){
//        on.what();
//        cout << "line : " << __LINE__ << endl;
            return false;
        }

    }

    std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_ballistic_deposition_L1";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        return s;
    }

};

/*********************
 *
 */
class SitePercolationBallisticDeposition_L2_v3: public SitePercolationBallisticDeposition_v3{
public:
    ~SitePercolationBallisticDeposition_L2_v3() = default;
    SitePercolationBallisticDeposition_L2_v3(value_type length, bool periodicity)
            :SitePercolationBallisticDeposition_v3(length, periodicity){}

    bool occupy() {
        // if no site is available then return false

        if(_number_of_occupied_sites == maxSites()){
            return false;
        }

        try {

//            value_type v = placeSite_2nn_v0();
            value_type v = placeSite_2nn_v1();
            _occuption_probability = occupationProbability(); // for super class

            return v != ULLONG_MAX;
        }catch (OccupiedNeighbor& on){
//        on.what();
//        cout << "line : " << __LINE__ << endl;
            return false;
        }

    }

    std::string getSignature() {
        std::string s = "sq_lattice_site_percolation_ballistic_deposition_L2";
        if(_periodicity)
            s += "_periodic_";
        else
            s += "_non_periodic_";
        return s;
    }

};


#endif //SQLATTICEPERCOLATION_PERCOLATION_SITE_EXPLOSIVE_H
