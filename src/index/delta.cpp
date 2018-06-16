//
// Created by shahnoor on 29-Apr-18.
//

#include "delta.h"

using namespace std;

ostream& operator<<(ostream& os, const Delta& delta)
{
    return os <<'(' << delta.dx << ',' << delta.dy << ')';
//    return os << '(' << std::setw(2) << delta.dx << ',' << std::setw(2) << delta.dy << ')';
}


/**
 * dx = siteA.col() - siteB.col()
 * dy = siteA.row() - siteB.row()
 *
 * @param siteA
 * @param siteB
 * @return Delta{dx, dy}
 */
Delta getDelta(Index siteA, Index siteB){
    int dx = int(siteA.col_) - int(siteB.col_);
    int dy = int(siteA.row_) - int(siteB.row_);
    return Delta{dx, dy};
}

