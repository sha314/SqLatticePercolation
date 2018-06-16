//
// Created by shahnoor on 10/2/2017.
//

#ifndef SITEPERCOLATION_HELPER_H
#define SITEPERCOLATION_HELPER_H

#include <algorithm>
#include <vector>
#include <memory>
#include <map>

#include "types.h"
#include "index/index.h"
#include "lattice/bond.h"
#include "flags.h"

/**
 * Compare using second values in a pair or type <int, value_type>
 */
struct Comparator{
    bool operator()(std::pair<int, value_type> elem1 ,std::pair<int, value_type> elem2)
    {
        return elem1.second < elem2.second;
    };
};


/**
 *
 * @tparam A
 * @tparam B
 * @param p
 * @return
 */
template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
    return std::pair<B,A>(p.second, p.first);
}

/**
 *
 * @tparam A
 * @tparam B
 * @param src
 * @return
 */
template<typename A, typename B>
std::map<B,A> flip_map(const std::map<A,B> &src)
{
    std::map<B,A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()),
                   flip_pair<A,B>);
    return dst;
}


// todo Delete comments after Feb 2018

/**
 *
 */
//namespace Lattice_Helper{
//    bool is_neoghbor(Index pt1, Index pt2, value_type length);
//
//    std::vector<Index> get_neighbor(Index pt1, value_type length);

    /*
     * Checking if they are horizontal bond.
     * if needed rearrange the bond index ?? or do it when constructing Bond object ?
     */
//    bool is_h_bond(Bond bond, value_type length);   // todo // unnessesary
//    bool is_v_bond(Bond bond, value_type length);   // todo // unnessesary

//    // much slowest, as it checks all values one by one
//    void check_all_cluster(std::vector<std::vector<Bond>> &cluster);
//    void check_all_cluster(std::vector<std::vector<std::shared_ptr<Bond>>> &cluster);
//    // version 1 : slow, but faster than check_all_cluster()
//    // for L = 50  taken birthTime is 419.78 sec or 6.98 min
//    // for L = 100 taken birthTime is 25,929.2 sec or 7hrs 12 min
//    void check_cluster_v1(std::vector<std::vector<Bond>> &cluster, value_type first_row);
//
//    // version 2 is much faster than version 1, i.e, check_cluster_v1()
//    // for L = 50  taken birthTime is 2.249 sec
//    // for L = 100 taken birthTime is 29.999 sec
//    void check_cluster_v2(std::vector<std::vector<Bond>>& cluster,
//                                          value_type first_row, const std::vector<Bond>& bonds);
//}


template <typename T>
void view(const T& vec){
    std::cout << '{';
    for(auto a: vec){
        std::cout << a << ',';
    }
    std::cout << '}' << std::endl;
}


template <typename T>
void sort(std::vector<std::vector<T>>& vec){
    for(value_type i{} ; i != vec.size() ; ++i){
        std::sort(vec[i].begin(), vec[i].end());
    }
}


#endif //SITEPERCOLATION_HELPER_H
