//
// Created by shahnoor on 29-Apr-18.
//

#ifndef PERCOLATION_DELTA_H
#define PERCOLATION_DELTA_H

#include <ostream>
#include "index.h"

/**
 * Stores the difference in x and y
 */
struct Delta{
    int dx;
    int dy;

    ~Delta() = default;
    Delta(int delta_x, int delta_y):dx{delta_x},dy{delta_y} {}

};

std::ostream& operator<<(std::ostream& os, const Delta& delta);

Delta getDelta(Index siteA, Index siteB);


#endif //PERCOLATION_DELTA_H
