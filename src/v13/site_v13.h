//
// Created by shahnoor on 2/15/21.
//

#ifndef SQLATTICEPERCOLATION_SITE_H
#define SQLATTICEPERCOLATION_SITE_H

#include <string>
#include <vector>
#include <iostream>
#include "element_v13.h"
#include "index_v13.h"

class Site_v13: public Element_v13{
    std::string classname = "Site_v13";
    std::vector<int> connecting_bond_ids;
    Index_v13 index;
    RelativeIndex_v13 relativeIndex;

public:
    Site_v13() = default;
    Site_v13(int row, int col){
        index = Index_v13(row, col);
        relativeIndex = RelativeIndex_v13();
    }
    void reset(){
        Element_v13::reset();
//        connecting_bond_ids.clear();
        relativeIndex = RelativeIndex_v13();
    }

    Index_v13 get_index(){ return index;}
    RelativeIndex_v13 get_relative_index(){ return relativeIndex;}
    int get_gid(){ return Element_v13::get_gid();}
    int get_id(){ return Element_v13::get_id();}
    std::string get_str(int formatt=0){
        std::stringstream ss;
        if (formatt == 1){
            ss << "[" << std::setw(5) << get_gid() << "," << std::setw(5) << get_id() << "]";
        }
        else{
            ss << std::setw(5) << get_id() <<  "(" << std::setw(5) << index.row() << ","
               << std::setw(5) << index.col() << ")";
        }

        return ss.str();
    }

    void add_connecting_bond(int bond_id){
        if (connecting_bond_ids.size() >= 4){
            std::cout << "A site connects only four bonds" << std::endl;
            return;
        }
        connecting_bond_ids.push_back(bond_id);
    }

    void add_connecting_bonds(std::vector<int>& bonds){ // TODO ***
        if (connecting_bond_ids.size() >= 4){
            std::cout << "before : A site connects only four bonds" << std::endl;
            return;
        }

        connecting_bond_ids.insert(connecting_bond_ids.end(), bonds.begin(), bonds.end());
        if (connecting_bond_ids.size() > 4){
            std::cout << "after  : A site connects only four bonds" << std::endl;
            return;
        }
    }

    std::vector<int> get_connecting_bonds(){ return connecting_bond_ids;}

    void init_relative_index(){relativeIndex = RelativeIndex_v13(0, 0);}

    void set_relative_index(RelativeIndex_v13 ri){relativeIndex = ri;}

    bool is_root(){return relativeIndex.x_coord() == 0 && relativeIndex.y_coord() == 0;}
};


#endif //SQLATTICEPERCOLATION_SITE_H
