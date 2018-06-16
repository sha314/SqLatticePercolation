//
// Created by shahnoor on 10/2/2017.
//

#include <iostream>
#include <fstream>

#include "helper.h"

using namespace std;

double total_time_for_cluster{};

// todo Delete comments after Feb 2018
//
///**
// *
// * @param pt1
// * @param length
// * @return
// */
//std::vector<Index> Lattice_Helper::get_neighbor(Index pt1, value_type length)
//{
//    std::vector<Index> neighbor(4);
//    value_type x = (pt1.x_+1) % length;
//    neighbor[0] = {x, pt1.y_};
//    x = pt1.x_ == 0 ? length-1 : pt1.x_ - 1 < 0 ? 0 : pt1.x_-1;
////        row_ = (id.row_ - 1) < 0 ? 0 : ID.row_ - 1;
//    neighbor[1] = {x, pt1.y_};
//    value_type y = (pt1.y_ + 1) % length;
//    neighbor[2] = {pt1.x_, y};
//    y = pt1.y_ == 0 ? length-1 : pt1.y_ - 1 < 0 ? 0 : pt1.y_-1;
////        col_ = (id.col_ - 1) < 0 ? 0 : ID.col_ - 1;
//    neighbor[3] = {pt1.x_, y};
//    return neighbor;
//}
//
///**
// *
// * @param pt1
// * @param pt2
// * @param length
// * @return
// */
//bool Lattice_Helper::is_neoghbor(Index pt1, Index pt2, value_type length)
//{
//    std::vector<Index> neighbor = get_neighbor(pt1, length);
//    if(neighbor.size() != 4){
//        cout << "_neighbor.numberOfBonds() != 4 : line " << __LINE__ << endl;
//    }
//    for(size_t i{} ; i != neighbor.size() ; ++i){
//        if(pt2.x_ == neighbor[i].x_ && pt2.y_ == neighbor[i].y_)
//            return true;
//    }
//    return false;
//}
//
//
//// (0,0) <--> (0,1) and (0,1) <--> (0,0) are same h bond with ID (0,0)
//// (0,0) <--> (0,L-1) and (0,L-1) <--> (0,0)  are also same h bond with ID (0,L-1)
//// because ID is also used to index them
//// can we use std::map on this??
///**
// *
// * @param bond
// * @param length
// * @return
// */
//bool Lattice_Helper::is_h_bond(Bond bond, value_type length)
//{
//    if(bond._end1.x_ == bond._end2.x_){
//        bool b1 =  ((bond._end1.y_ + 1)%length) == bond._end2.y_;
//        bool b2 =  bond._end1.y_ == ((bond._end2.y_ + 1)%length);
//            return b1 || b2;
//    }
//
//    return false;
//}
//
///**
// *
// * @param bond
// * @param length
// * @return
// */
//bool Lattice_Helper::is_v_bond(Bond bond, value_type length)
//{
//
//    if(bond._end1.y_  == bond._end2.y_){
//        bool b1 = ((bond._end1.x_ + 1)%length) == bond._end2.x_;
//        bool b2 =  bond._end1.x_ == ((bond._end2.x_ + 1)%length);
//        return b1 || b2;
//    }
//    return false;
//}
//
//
//
//void Lattice_Helper::check_all_cluster(std::vector<std::vector<Bond>> &cluster)
//{
//    clock_t t = clock();
//    if(debug_check_all_cluster){
//        std::cout << "check_all_cluster" << std::endl;
//    }
//    // we must use less than symbol here. because numberOfBonds of cluster might be decreasing.
//    for(value_type i{} ; i < cluster.size() ; ++i){
//        for(value_type j{} ; j < cluster[i].size() ; ++j){
//            for(value_type k{j+1} ; k < cluster[i].size() ; ++k){
//                if(cluster[i][j] == cluster[i][k]){
//                    if(debug_check_all_cluster) {
//                        std::cout << "Repeated [" << i << "][" << j << "] and ["
//                                  << i << "][" << k << "]" << cluster[i][j] << std::endl;
//                    }
//                    cluster[i].erase(cluster[i].begin() + k);
//                }
//            }
//        }
//    }
//    if(time_check_all_cluster){
//        double t2 = (clock() - t)/ double(CLOCKS_PER_SEC);
//        total_time_for_cluster += t2;
//        cout << "birthTime required for cluster checking = " << t2 << " sec" << endl;
//        cout << "Total birthTime till now " << total_time_for_cluster << " sec" << endl;
//    }
//}
//
//
//void Lattice_Helper::check_all_cluster(std::vector<std::vector<shared_ptr<Bond>>> &cluster)
//{
//    std::cout << "not implemented in helper.cpp : line " << __LINE__ << std::endl;
////    clock_t t = clock();
////    if(debug_check_all_cluster){
////        std::cout << "check_all_cluster" << std::endl;
////    }
////    // we must use less than symbol here. because numberOfBonds of cluster might be decreasing.
////    for(value_type i{} ; i < cluster.numberOfBonds() ; ++i){
////        for(value_type j{} ; j < cluster[i].numberOfBonds() ; ++j){
////            for(value_type k{j+1} ; k < cluster[i].numberOfBonds() ; ++k){
////                if(cluster[i][j] == cluster[i][k]){
////                    if(debug_check_all_cluster) {
////                        std::cout << "Repeated [" << i << "][" << j << "] and ["
////                                  << i << "][" << k << "]" << cluster[i][j] << std::endl;
////                    }
////                    cluster[i].erase(cluster[i].begin() + k);
////                }
////            }
////        }
////    }
////    if(time_check_all_cluster){
////        double t2 = (clock() - t)/ double(CLOCKS_PER_SEC);
////        total_time_for_cluster += t2;
////        cout << "birthTime required for cluster checking = " << t2 << " sec" << endl;
////        cout << "Total birthTime till now " << total_time_for_cluster << " sec" << endl;
////    }
//}
//
//
//// reduces birthTime. but not significantly than check_all_cluster()
//void Lattice_Helper::check_cluster_v1(std::vector<std::vector<Bond>> &cluster, value_type first_row)
//{
//
//    clock_t t = clock();
//    if(debug_check_cluster){
//        std::cout << "check_cluster_v1" << std::endl;
//    }
//    // we must use less than symbol here. because numberOfBonds of cluster might be decreasing.
//    for (value_type j{}; j < cluster[first_row].size(); ++j) {
//        for (value_type k{j + 1}; k < cluster[first_row].size(); ++k) {
//            if (cluster[first_row][j] == cluster[first_row][k]) {
//                if (debug_check_cluster) {
//                    std::cout << "Repeated [" << first_row << "][" << j << "] and ["
//                              << first_row << "][" << k << "]" << cluster[first_row][j] << std::endl;
//                }
//                cluster[first_row].erase(cluster[first_row].begin() + k);
//            }
//        }
//    }
//
//    if(time_check_cluster){
//
//        double t2 = (clock() - t)/ double(CLOCKS_PER_SEC);
//        total_time_for_cluster += t2;
//        cout << "birthTime required for cluster checking = " << t2 << " sec" << endl;
//        cout << "Total birthTime till now " << total_time_for_cluster << " sec" << endl;
//    }
//}
//
//
//
//// reduces birthTime significantly than check_all_cluster()
//// version 2
//void Lattice_Helper::check_cluster_v2(std::vector<std::vector<Bond>>& cluster,
//                                      value_type first_row, const std::vector<Bond>& bonds)
//{
//    clock_t t = clock();
//    if(debug_check_cluster_v2){
//        std::cout << "check_cluster_v1" << std::endl;
//    }
//    if(bonds.size() > 4){
//        cout << "bonds.numberOfBonds() > 4 : error : line " << __LINE__ << endl;
//    }
//    if(debug_check_cluster_v2) {
//        cout << "Bonds : ";
//        view(bonds);
//    }
//    // we must use less than symbol here. because numberOfBonds of cluster might be decreasing.
//    std::vector<value_type> index;
//    for (value_type j{}; j != bonds.size(); ++j) {
//        index.clear();
//        for (value_type k{}; k != cluster[first_row].size(); ++k) {
//            if (bonds[j] == cluster[first_row][k]) {
//                if (debug_check_cluster_v2) {
//                    std::cout << "Repeated [" << first_row << "][" << k << "] "
//                              << cluster[first_row][k] << std::endl;
//                }
//                index.push_back(k);
//            }
//        }
//        if(index.size() > 1){
//            // if index.numberOfBonds() == 1, it means only one is present. so we don't need to erase it
//            // if index.numberOfBonds() == 2, we need to erase one
//            // if index.numberOfBonds() > 2,
//            if(debug_check_cluster_v2) {
//                cout << "index.numberOfBonds() = " << index.size() << " : line " << __LINE__ << endl;
//
//                view(index);
//                cout << "cluster before : " << endl;
//                view(cluster[first_row]);
//            }
//            for(value_type m{1} ; m!= index.size() ; ++m){
//                cluster[first_row].erase(cluster[first_row].begin() + index[m]);
//            }
//            if(debug_check_cluster_v2) {
//                cout << "cluster after : " << endl;
//                view(cluster[first_row]);
//            }
//
//        }
//
//    }
//
//
//
//    if(time_check_cluster_v2){
//        double t2 = (clock() - t)/ double(CLOCKS_PER_SEC);
////        fout << t2 << endl;
//        total_time_for_cluster += t2;
//        cout << "birthTime required for cluster checking = " << t2 << " sec" << endl;
//        cout << "Total birthTime till now " << total_time_for_cluster << " sec" << endl;
//    }
//}


