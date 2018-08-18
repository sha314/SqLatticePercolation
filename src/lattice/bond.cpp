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
