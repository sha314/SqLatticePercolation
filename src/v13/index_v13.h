//
// Created by shahnoor on 2/15/21.
//

#ifndef SQLATTICEPERCOLATION_INDEX_H
#define SQLATTICEPERCOLATION_INDEX_H


#include <string>
#include <sstream>
#include <iomanip>

class Index_v13{
protected:
    int _component_1=-1;
    int _component_2=-1;

public:
    Index_v13() = default;

    Index_v13(int row, int col){
        _component_1 = row;
        _component_2 = col;
    }

    Index_v13(Index_v13& index) = default;

    Index_v13(Index_v13&& index) = default;

    Index_v13& operator=(const Index_v13&) = default;

    Index_v13 operator+(const Index_v13& other){
        int del_row = _component_1 + other._component_1;
        int del_col = _component_2 + other._component_2;
        return {del_row, del_col};
    }

    Index_v13 operator-(Index_v13 const& second){
        int del_row = _component_1 - second._component_1;
        int del_col = _component_2 - second._component_2;
        return {del_row, del_col};
    }

    Index_v13 subtract(const Index_v13& second){
        int del_row = _component_1 - second._component_1;
        int del_col = _component_2 - second._component_2;
        return Index_v13(del_row, del_col);
    }

    Index_v13 add(const Index_v13 & second){
        int del_row = _component_1 + second._component_1;
        int del_col = _component_2 + second._component_2;
        return {del_row, del_col};
    }

    int row(){ return _component_1;}
    int col(){ return _component_2;}

    virtual std::string get_str(){
        std::stringstream ss;
        ss << "(" << std::setw(3) << _component_1 << "," << std::setw(3) << _component_2 << ")";
        return ss.str();
    }
};

class RelativeIndex_v13: public Index_v13{


public:
    RelativeIndex_v13() : Index_v13(-1, -1){}

    RelativeIndex_v13(int x_rel, int y_rel):Index_v13(x_rel, y_rel){

    }

    RelativeIndex_v13(RelativeIndex_v13&& relI):Index_v13(relI.x_coord(), relI.y_coord()) {}
    RelativeIndex_v13(RelativeIndex_v13& relI):Index_v13(relI.x_coord(), relI.y_coord()) {}
    RelativeIndex_v13& operator=(const RelativeIndex_v13& relI) = default;

    explicit RelativeIndex_v13(Index_v13 index):Index_v13(index){

    }

    int x_coord() {
        return row();
    }

    int y_coord() {
        return col();
    }

    std::string get_str() override {
        std::stringstream ss;
        ss << "<" << std::setw(3) << x_coord() << "," << std::setw(3) << y_coord() << ">";
        return ss.str();
    }

    RelativeIndex_v13 operator-(RelativeIndex_v13& other){
        int del_row = x_coord() - other.x_coord();
        int del_col = y_coord() - other.y_coord();
        return {del_row, del_col};
    }

};


#endif //SQLATTICEPERCOLATION_INDEX_H
