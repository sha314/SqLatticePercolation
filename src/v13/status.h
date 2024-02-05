//
// Created by shahnoor on 6/28/21.
//

#ifndef SQLATTICEPERCOLATION_STATUS_H
#define SQLATTICEPERCOLATION_STATUS_H


enum class P_STATUS{
    SUCESS,
    CURRENT_SITE_NOT_EMPTY, // there are empty sites but the current selected site is not empty
    EMPTY_SITE_LIST,         // no free site available. for site percolation
    EMPTY_BOND_LIST         // no free bond is available. for bond percolation

};

#endif //SQLATTICEPERCOLATION_STATUS_H
