//
// Created by shahnoor on 11/2/2017.
//


#include "bond.h"

/**
 * use '-' and '|' in between '<>' to indicate horizontal or vertical bond
 * ((0,1)<->(0,0)) for horizontal bond
 * ((1,1)<|>(0,1)) for vertical bond
 * @param os
 * @param bond
 * @return
 */
std::ostream& operator<<(std::ostream& os, const Bond& bond)
{
    if(bond.isActive()) {
        // place '-' for horizontal bond and '|' for vertical bong
        if(bond.isHorizontal()) {
            return os << '(' << bond._end1 << "<->" << bond._end2 << ')';
        }
        return os << '(' << bond._end1 << "<|>" << bond._end2 << ')';
    }
    else
        return os << "(**)";
}

/**
 * use '-' and '|' in between '<>' to indicate horizontal or vertical bond
 * ((0,1)<->(0,0)) for horizontal bond
 * ((1,1)<|>(0,1)) for vertical bond
 * @param os
 * @param bond
 * @return
 */
std::ostream& operator<<(std::ostream& os, const Bond_v2& bond) {
    if (bond.isActive()) {
        // place '-' for horizontal bond and '|' for vertical bong
        return os << '(' << bond.getIndex() << ')';
    } else
        return os << "(**)";
}

bool operator==(Bond a, Bond b)
{
    if(a.isHorizontal() && b.isHorizontal())
    {
        return (a.id().row_ == b.id().row_) && (a.id().column_ == b.id().column_);
    }
    if(a.isVertical() && b.isVertical()){
        return (a.id().row_ == b.id().row_) && (a.id().column_ == b.id().column_);
    }

    return false;
}

bool  operator<(const Bond& bond1, const Bond& bond2){
    if(bond1.isHorizontal() && bond2.isHorizontal()){
        return bond1._end1.column_ < bond2._end1.column_;
    }
    if(bond1.isVertical() && bond2.isVertical()){
        return bond1._end1.row_ < bond2._end1.row_;
    }
    return bond1.isHorizontal();
}


bool  operator>(const Bond& bond1, const Bond& bond2){
    if(bond1.isHorizontal() && bond2.isHorizontal()){
        return bond1._end1.column_ > bond2._end1.column_;
    }
    if(bond1.isVertical() && bond2.isVertical()){
        return bond1._end1.row_ > bond2._end1.row_;
    }
    return bond1.isVertical();
}

/**********************************
 * Version 2 of class Bond
 */
Index Bond_v2::getSiteA() const {
    return Index(_id.row_, _id.column_);
}

Index Bond_v2::getSiteB(size_t length) const {
    auto index = Index(_id.row_, _id.column_);
    if(_id.bondType == BondType::Horizontal){
        index.column_ = (index.column_ + 1) % length;
    }else if(_id.bondType == BondType::Vertical){
        index.row_ = (index.row_ + 1) % length;
    }
    return index;
}

std::vector<Index> Bond_v2::getSiteIndex(size_t length) const {
    auto index = Index(_id.row_, _id.column_);
    Index index2 = index;
    if(_id.bondType == BondType::Horizontal){
        index2.column_ = (index.column_ + 1) % length;
    }else if(_id.bondType == BondType::Vertical){
        index2.row_ = (index.row_ + 1) % length;
    }

    return {index, index2};
}
