//
// Created by shahnoor on 10/2/2017.
//

#include <iomanip>
#include "index.h"

using namespace std;


ostream& operator<<(ostream& os, const Index& index)
{
    return os <<'(' << index.row_ << ',' << index.column_ << ')';
//    return os << '(' << std::setw(2) << index.x_ << ',' << std::setw(2) << index.y_ << ')';
}

ostream& operator<<(ostream& os, const IndexRelative& index)
{
//    return os <<'(' << index.x_ << ',' << index.y_ << ')';
    return os << '(' << std::setw(3) << index.x_ << ',' << std::setw(3) << index.y_ << ')';
}

bool operator==(const Index& index1, const Index& index2){
    return (index1.row_ == index2.row_) && (index1.column_ == index2.column_);
}

bool operator<(const Index& index1, const Index& index2){
    if(index1.row_ < index2.row_)
            return true;
    if(index1.row_ == index2.row_){
        return index1.column_ < index2.column_;
    }
    return false;
}


ostream& operator<<(ostream& os, const BondIndex& index){
    if(index.horizontal()){
        // horizontal bond
        os << "<->" ;
    }
    if (index.vertical()){
        // vertical bond
        os << "<|>" ;
    }
    return os << '(' << index.row_ << ',' << index.column_ << ')';
}

bool operator==(BondIndex index1, BondIndex index2){
    if(index1.horizontal() ==  index2.horizontal() || index1.vertical() ==  index2.vertical()){
        // horizontal or vertical
        return index1.row_ == index2.row_ && index1.column_ == index2.column_;
    }
    return false;
}


bool operator<(BondIndex index1, BondIndex index2){
    cout << "not yet defined : line " << __LINE__ << endl;
    return false;
}


/**
 * Get the 2nd nearest nearest neighbor in the direction of 1st nearest neighbor.
 * Periodicity is not considered here.
 * @param center
 * @param nn_1
 * @param length
 * @return
 */
Index get_2nn_in_1nn_direction(Index center, Index nn_1, value_type length){
    int delta_c = int(nn_1.column_) - int(center.column_);
    int delta_r = int(nn_1.row_) - int(center.row_);
    if (delta_c == 0 && delta_r == 0){
        cout << "Both indices are same : line " << __LINE__ << endl;
    }
    else if(delta_c > 1 || delta_r > 1){
        // meaning, the site_index_sequence are on the opposite edges
//        cout << "2nd index is not the First nearest neighbor : line " << __LINE__ << " : file " << __FILE__ << endl;
    }

    return Index{(nn_1.row_ + delta_r + length) % length, (nn_1.column_ + delta_c + length) % length};
}


/**
 * Get all second nearest neighbors based on the first nearest neighbors.
 * Periodicity is not considered here
 * @param center
 * @param nn_1
 * @param length
 * @return
 */
vector<Index> get_all_2nn_in_1nn_s_direction(Index center, const vector<Index> &nn_1, value_type length){
    vector<Index> nn_2(nn_1.size());

    for(size_t i{}; i != nn_1.size() ; ++i){
        int delta_c = int(nn_1[i].column_) - int(center.column_);
        int delta_r = int(nn_1[i].row_) - int(center.row_);
        if (delta_c == 0 && delta_r == 0){
            cout << "Both indices are same : line " << __LINE__ << endl;
        }
        else if(delta_c > 1 || delta_r > 1){
            // meaning, the site_index_sequence are on the opposite edges
//            cout << "center " << center << " nn " << nn_1 << endl;
//            cout << "2nd index is not the First nearest neighbor : line " << __LINE__ << " : file " << __FILE__ << endl;
        }

        nn_2[i] =  Index{(nn_1[i].row_ + delta_r + length) % length, (nn_1[i].column_ + delta_c + length) % length};
    }
    return nn_2;
}

Index IndexTranslator::translate1DToSite(size_t length, index_type index) {
    size_t c = index % length;
    size_t r = index / length;
    return Index(r, c);
}

index_type IndexTranslator::translateSiteTo1D(size_t length, Index index) {
    index_type i  = length * index.row_ + index.column_;
    return i;
}

BondIndex IndexTranslator::translate1DToBond(size_t length, index_type index) {
    size_t  length_squared = length * length;
    size_t t = index / length_squared; // determine horizontal or vertical
    size_t half = index - t*length_squared;
    size_t r = half / length;
    size_t c = half % length;
    auto bt = (t == 0)? BondType::Horizontal : BondType::Vertical;
    return BondIndex(bt, r, c);
}

index_type IndexTranslator::translateBondTo1D(size_t length, BondIndex bondIndex) {
    index_type  length_squared = length * length;
    index_type t = (bondIndex.bondType == BondType ::Horizontal) ? 0 : 1;
    index_type index  = length_squared * t + length * bondIndex.row_ + bondIndex.column_;

    return index;
}
