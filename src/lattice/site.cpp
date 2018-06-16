//
// Created by shahnoor on 10/2/2017.
//

#include <iomanip>
#include "site.h"


std::ostream& operator<<(std::ostream& os, const Site& site)
{
    if(site.isActive())
        return os << site._id;
    else
        return os << "(*)";
}



bool operator==(Site& site1, Site& site2){
    return (site1.ID().row_ == site2.ID().row_) && (site1.ID().col_ == site2.ID().col_);
}
