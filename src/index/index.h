//
// Created by shahnoor on 10/2/2017.
//

#ifndef SITEPERCOLATION_INDEX_H
#define SITEPERCOLATION_INDEX_H


#include <ostream>
#include <iostream>
#include <sstream>
#include <vector>

#include "../types.h"
#include "../exception/exceptions.h"
#include "../lattice/bond_type.h"


/**
 *
 */
struct Index{
    value_type row_{};
    value_type col_{};

    ~Index()                      = default;
    Index()                       = default;

//    Index(Index&)               = default;
//    Index(Index&&)              = default;
//    Index& operator=(Index&)    = default;
//    Index& operator=(Index&&)   = default;

    Index(value_type x, value_type y) : row_{x}, col_{y} {}

};


/**
 * Template index
 * @tparam T
 */
class IndexRelative{
public:
    int x_{};
    int y_{};

    ~IndexRelative()                      = default;
    IndexRelative()                       = default;

    IndexRelative(int x, int y) : x_{x}, y_{y} {}

};

/**
 *
 */
struct BondIndex{
    // if hv_ == 0 -> horizontal bond
    // if hv_ == 1 -> vertical bond
    BondType bondType;
//    value_type hv_;
    value_type row_;
    value_type column_;

    ~BondIndex()                        = default;
    BondIndex()                         = default;
//    BondIndex(BondIndex&)               = default;
//    BondIndex(BondIndex&&)              = default;
//    BondIndex& operator=(BondIndex&)    = default;
//    BondIndex& operator=(BondIndex&&)   = default;

//    BondIndex(value_type hv, value_type row, value_type column)
//            : hv_{hv}, x_{row}, column_{column}
//    {
//        if(hv_ != 0 && hv_ != 1){
//            std::cout << "Cannot construct object. throwing exception : line " << __LINE__ << std::endl;
//            throw InvalidIndex{"hv_ != 0 || hv_ != 1 : line " + std::to_string(__LINE__) };
//        }
//    }


    BondIndex(BondType hv, value_type row, value_type column)
            :  row_{row}, column_{column}
    {
        bondType = hv;
    }

    bool horizontal() const { return bondType == BondType::Horizontal;}
    bool vertical() const { return bondType == BondType::Vertical;}
};


std::ostream& operator<<(std::ostream& os, const Index& index);
bool operator==(const Index& index1, const Index& index2);
bool operator<(const Index& index1, const Index& index2);

std::ostream& operator<<(std::ostream& os, const IndexRelative& index);

std::ostream& operator<<(std::ostream& os, const BondIndex& index);
bool operator==(BondIndex index1, BondIndex index2);
bool operator<(BondIndex index1, BondIndex index2);


/**
 *  Get 2nd nearest neightbor / sin the direction of 1st nearest neighbor, while @var center is the center
 */
Index get_2nn_in_1nn_direction(Index center, Index nn_1, value_type length);
std::vector<Index> get_2nn_s_in_1nn_s_direction(Index center, const std::vector<Index> &nn_1, value_type length);

#endif //SITEPERCOLATION_INDEX_H
