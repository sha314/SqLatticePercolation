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
