//
// Created by shahnoor on 2/15/21.
//

#ifndef SQLATTICEPERCOLATION_BOND_H
#define SQLATTICEPERCOLATION_BOND_H

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include "element_v13.h"

class Bond_v13: public Element_v13 {
    int _row = -1;
    int _col = -1;
    int _type = -1; // 0 means horizontal and 1 means vertical
    std::vector<int> connected_site_ids;
    std::string classname = "Bond_v13";

public:
    Bond_v13() = default;
    Bond_v13(int row, int col, int type):Element_v13(){
        _row = row;
        _col = col;
        _type = type;
    }

    void reset(){
        Element_v13::reset();
        connected_site_ids.clear();
    }

    std::string get_str(int formatt=0){
        std::stringstream ss;
        if (formatt == 1){
            ss << "[" << std::setw(5) << get_gid() << "," << std::setw(5) << get_id() << "]";
        } else if(formatt == 2){
            ss << std::setw(5) << get_id() << "(" << std::setw(5) << _row << "," << std::setw(5) << _col << ")";
        }else{
            ss << "(" << std::setw(5) << _row << "," << std::setw(5) << _col << ", " << _type << ")";
        }

        return ss.str();
    }

    void add_connected_site(int site_id){
        if (connected_site_ids.size() >= 2){
            std::cout << "A bond connects only two sites" << std::endl;
            return;
        }
        connected_site_ids.push_back(site_id);
//        if (typeid(site_id).hash_code() == typeid(connected_site_ids).hash_code()){
//            connected_site_ids.insert(connected_site_ids.end(), site_id);
//        }else{
//
//        }
    }

    std::vector<int> get_connected_sites(){ return connected_site_ids;}

};



#endif //SQLATTICEPERCOLATION_BOND_H
