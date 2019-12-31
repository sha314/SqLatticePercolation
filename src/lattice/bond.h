//
// Created by shahnoor on 10/2/2017.
//

#ifndef SITEPERCOLATION_BOND_H
#define SITEPERCOLATION_BOND_H


#include <ostream>
#include <iostream>
#include <sstream>
#include <strings.h>

#include "../index/index.h"
#include "../types.h"
#include "../index/id.h"



/**
 * A bond has two end
 * say a 5x5 lattice bond between end1 (0,0) and end2 (0,1)
 * if _status is false -> bond is not there
 *
 */
struct Bond{
    // check if active or not
    bool _status{false};
    value_type _length;
    int _group_id{-1};

    BondType bondType;
    //relative distance from the root site. {0,0} if it is the root site
    //very useful for detecting wrapping // todo
    IndexRelative _relative_index{0,0};

//    bool horizontal{false};
//    bool vertical{false};

    Index _end1;
    Index _end2;

    BondIndex _id;

    ~Bond() = default;
    Bond() = default;
    Bond(Index end1, Index end2, value_type length){
        _end1.row_ = end1.row_;
        _end1.column_ = end1.column_;
        _end2.row_ = end2.row_;
        _end2.column_ = end2.column_;

        // correct the bond here
        _length = length;
        // check if the bond is valid??
        if(_end1.row_ == _end2.row_){
            bondType = BondType::Horizontal;
//            horizontal = true;
            // means x_ values are equal
            if(_end1.column_ > _end2.column_){
                // case (0,0)<->(0,4) for L=5 and
                // case (0,4)<->(0,0) are equal for bond like this end1=(0,4) end2=(0,0)
                // for all other cases lower index is end1
                if(_end1.column_ == _length-1 && _end2.column_ ==0){
                    // do nothing
//                    std::cout << "_end1.y_ == _length-1 && _end2.y_ ==0 : line " << __LINE__ << std::endl;
                }
                else{
                    // sort them out
                    _end1.column_ = end2.column_;
                    _end2.column_ = end1.column_;
                }
            }
            else if(_end1.column_ < _end2.column_){
                if(_end1.column_ == 0 && _end2.column_ == _length-1){
                    _end1.column_ = end2.column_;
                    _end2.column_ = end1.column_;
                }
            }
//            _id = set_ID<3>({0, _end1.x_, _end1.y_});
        }
        else if(_end1.column_ == _end2.column_){
            bondType = BondType::Vertical;
//            vertical = true;
            // means y_ values are equal
            if(_end1.row_ > _end2.row_){
                // case (0,0)<->(4,0) for L=5 and
                // case (4,0)<->(0,0) are equal for bond like this end1=(4,0) end2=(0,0)
                // for all other cases lower index is end1
                if(_end1.row_ == _length-1 && _end2.row_ ==0){
                    // do nothind
//                    std::cout << "_end1.x_ == _length-1 && _end2.x_ ==0 : line " << __LINE__ << std::endl;
                }
                else{
                    // sort them out
                    _end1.row_ = end2.row_;
                    _end2.row_ = end1.row_;
                }
            }
            else if(_end1.row_ < _end2.row_){
                if(_end1.row_ == 0 && _end2.row_ == _length-1){
                    _end1.row_ = end2.row_;
                    _end2.row_ = end1.row_;
                }
            }
//            _id = set_ID<3>({1, _end1.x_, _end1.y_});
        }
        else{
            std::cout << '(' << _end1.row_ << ',' << _end1.column_ << ')' << "<->"
                    << '(' << _end2.row_ << ',' << _end2.column_ << ')'
                    << " is not a valid bond : line " << __LINE__ << std::endl;
//            _id = set_ID<3>({-1, -1, -1});
        }

//        _id = {(horizontal) ? 0ul : 1ul, _end1.x_, end1.y_};  // unsigned long
        _id = BondIndex(bondType, _end1.row_, _end1.column_);  // unsigned long

    }

    Bond(BondIndex bondIndex, value_type length){
        _id = bondIndex;

    }

    std::vector<Index> getSites() const { return {_end1, _end2};}

    Index id() const {
        return _end1;
    }

    BondIndex ID() const {
        return _id;
    }

    void activate() {_status = true;}
    void deactivate() {
        _relative_index = {0,0};
        _group_id = -1;
        _status = false;
    }
    bool isActive() const { return _status;}
/*
* Group get_ID is the set_ID of the cluster they are in
*/
    int get_groupID() const {return _group_id;}
    void set_groupID(int g_id) {_group_id = g_id;}

//    std::stringstream getBondString() const {
//        std::stringstream ss;
//        if(isActive()) {
//            // place '-' for horizontal bond and '|' for vertical bong
//            if(horizontal) {
//                ss << '(' << _end1 << "<->" << _end2 << ')';
//            }
//            else {
//                ss << '(' << _end1 << "<|>" << _end2 << ')';
//            }
//        }
//        else
//            ss << "(**)";
//        return ss;
//    }

    std::stringstream getBondString() const {
        std::stringstream ss;
        if(isActive()) {
            // place '-' for horizontal bond and '|' for vertical bong
            if(bondType == BondType::Horizontal) {
                ss << '(' << _end1 << "<->" << _end2 << ')';
            }
            else {
                ss << '(' << _end1 << "<|>" << _end2 << ')';
            }
        }
        else
            ss << "(**)";
        return ss;
    }

    bool isHorizontal() const { return bondType == BondType ::Horizontal;}
    bool isVertical()   const { return bondType == BondType ::Vertical;}

    void relativeIndex(IndexRelative r){
        _relative_index = r;
    }

    void relativeIndex(int x, int y){
        _relative_index = {x,y};
    }

    IndexRelative relativeIndex() const {return _relative_index;}
};


 /**
  * Created on : 2019.06.29
  * created by : shahnoor
  * version 2 of Bond
  */
struct Bond_v2{
    // check if active or not
    bool _status{false};
    int _group_id{-1};

    BondIndex _id{};
    ~Bond_v2() = default;
    Bond_v2() = default;
    Bond_v2(BondIndex bondIndex){
        _id = bondIndex;
    }
    Index getSiteA()const;
    Index getSiteB(size_t length)const;
    std::vector<Index> getSiteIndex(size_t length) const ;

    BondIndex getIndex() const {
        return _id;
    }

    void activate() {_status = true;}
    void deactivate() {
        _group_id = -1;
        _status = false;
    }
    bool isActive() const { return _status;}
/*
* Group get_ID is the set_ID of the cluster they are in
*/
    int get_groupID() const {return _group_id;}
    void set_groupID(int g_id) {_group_id = g_id;}


    bool isHorizontal() const { return _id.bondType == BondType ::Horizontal;}
    bool isVertical()   const { return _id.bondType == BondType ::Vertical;}

};



/*************************
 * Date : 2019.12.23
 * Author : Shahnoor
 *
 * Specially build for v12
 *
 * All bonds are identical. No distinction between horizontal and vertical.
 */
class Bond_v12{

    Link _index;
    int _id{-1};



    // check if active or not
    bool _status{false};
    int _group_id{-1};


public:

    ~Bond_v12() = default;
    Bond_v12() = default;
    Bond_v12(Link index, int id){
        _id = id;
        _index = index;
    }

    void reset(){
        _group_id = -1;
        _status = false;
    }

    int get_id() const {return _id;}
    Link getIndex() const {
        return _index;
    }

    std::vector<int> connectedSites() const {return {_index.nodeA(), _index.nodeB()};}
    void activate() {_status = true;}
    void deactivate() {
        _group_id = -1;
        _status = false;
    }
    bool isActive() const { return _status;}
/*
* Group get_ID is the set_ID of the cluster they are in
*/
    int get_groupID() const {return _group_id;}
    void set_groupID(int g_id) {_group_id = g_id;}


};



std::ostream&   operator<<(std::ostream& os, const Bond& bond);
std::ostream&   operator<<(std::ostream& os, const Bond_v2& bond);
std::ostream&   operator<<(std::ostream& os, const Bond_v12& bond);
bool            operator==(Bond a, Bond b);
bool            operator<(const Bond& bond1, const Bond& bond2);
bool            operator>(const Bond& bond1, const Bond& bond2);


#endif //SITEPERCOLATION_BOND_H
