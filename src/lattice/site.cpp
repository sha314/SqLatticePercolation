//
// Created by shahnoor on 10/2/2017.
//

#include <iomanip>
#include "site.h"

using namespace std;


std::ostream& operator<<(std::ostream& os, const Site& site)
{
    if(site.isActive())
        return os << site._id;
    else
        return os << "(*)";
}



bool operator==(Site& site1, Site& site2){
    return (site1.ID().row_ == site2.ID().row_) && (site1.ID().column_ == site2.ID().column_);
}


std::ostream& operator<<(std::ostream& os, const Site_v12& site)
{
    return os << site.get_id() << "," << site.get_index() << "," << site.get_groupID();
}

void Site_v12::addToRelativeIndex(int dx, int dy) {
    cout << "old " << _relative_index << endl;
    _relative_index.x_ += dx;
    _relative_index.y_ += dy;
    cout << "new " << _relative_index << endl;
}

void Site_v12::addToRelativeIndex(IndexRelative dx_dy) {
    cout << "old " << _relative_index << endl;
    _relative_index.x_ += dx_dy.x_;
    _relative_index.y_ += dx_dy.y_;
    cout << "new " << _relative_index << endl;
}
