//
// Created by shahnoor on 6/28/19.
//

#include <iomanip>
#include "lattice.h"
#include "../util/printer.h"

using namespace std;


/**
 *
 * @param length         -> length of the lattice
 * @param activate_bonds -> if true all bonds are activated by default
 *                          and will not be deactivated as long as the object exists,
 *                          even if SLattice::reset function is called.
 * @param activate_sites -> if true all sites are activated by default
 *                          and will not be deactivated as long as the object exists,
 *                          even if SLattice::reset function is called.
 */
SqLattice_v2::SqLattice_v2(value_type length) {
//    cout << "Constructing Lattice object : line " << __LINE__ << endl;
    _length = length;
    _length_squared = _length * _length;
    _sites = std::vector<Site>(_length_squared);
    _bonds = std::vector<Bond_v2>(2 * _length_squared);
}

/********************************************************************
 * Activation and Deactivation
 ********************************************************************/
void SqLattice_v2::activateAllSite()
{
    for(value_type i{} ; i != _length_squared ; ++i) {
        _sites[i].activate();
    }
}

/**
 *
 */
void SqLattice_v2::activateAllBond()
{
    for(value_type i{} ; i != _length_squared ; ++i) {
        _bonds[i].activate();
        _bonds[i  + _length_squared].activate();
    }
}

void SqLattice_v2::activate_site(Index index) {
//    cout << "activating site " << index << endl;
    size_t i = translateSiteTo1D(index);
    _sites[i].activate();
}


void SqLattice_v2::activateBond(BondIndex bond) {
    // check if the bond is vertical or horizontal
    // then call appropriate function to activate _h_bond or _v_bond
    size_t i = translateBondTo1D(bond);
    if(_bonds[i].isActive()){
        cout << "Bond is already activated : line " << __LINE__ << endl;
    }
    _bonds[i].activate();
}


void SqLattice_v2::deactivate_site(Index index){
    size_t i = translateSiteTo1D(index);
    _sites[i].deactivate();
}


void SqLattice_v2::deactivate_bond(BondIndex bond) {
    // check if the bond is vertical or horizontal
    // then call appropriate function to activate _h_bond or _v_bond
    size_t i = translateBondTo1D(bond);
    _bonds[i].deactivate();
}


/*********************************************************
 * Viewing methods
 *********************************************************/
/**
 *  View the sites of the lattice
 *  place (*) if the site is not active
 */
void SqLattice_v2::view_sites()
{
    std::cout << "view site_index_sequence" << std::endl;
    std::cout << '{';
    for(value_type r{} ; r != _length ; ++r) {
        if(r!=0) std::cout << "  ";
        else std::cout << '{';
        for (value_type c{}; c != _length; ++c) {
            value_type i_site = r * _length + c;
            if(_sites[i_site].isActive()){
                std::cout << _sites[i_site] ;
            }
            else{
                std::cout << "(*)";
            }
            if(c != _length-1)
                std::cout << ',';
        }
        std::cout << '}';
        if(r != _length-1)
            std::cout << std::endl;
    }
    std::cout << '}';
    std::cout << std::endl;
}

/**
 *  View the sites of the lattice
 *  place (*) if the site is not active
 *  Shows the group_id along with sites
 *
 *  Very good output format. Up to lattice size < 100
 */
void SqLattice_v2::view_sites_extended()
{
    std::cout << "view site_index_sequence" << std::endl;
    std::cout << '{';
    for(value_type r{} ; r != _length ; ++r) {
        if(r!=0) std::cout << "  ";
        else std::cout << '{';
        for (value_type c{}; c != _length; ++c) {
            value_type i_site = r * _length + c;
            std::cout << std::setw(3) << _sites[i_site].get_groupID() << ":";

            if(_sites[i_site].isActive()) {
                cout << '(' << std::setw(2) << _sites[i_site]._id.row_ << ','
                     << std::setw(2) << _sites[i_site]._id.column_ << ')';
            }
            else{
                cout << std::setw(7) << "(*)";
            }

            if(c != _length-1)
                std::cout << ',';
        }
        std::cout << '}';
        if(r != _length-1)
            std::cout << std::endl;
    }
    std::cout << '}';
    std::cout << std::endl;
}



/**
 * Displays group ids of sites in a matrix form
 */
void SqLattice_v2::view_sites_by_id() {
    std::cout << "Sites by id : line " << __LINE__ << endl;
    cout << "   ";
    for(value_type j{}; j != _length; ++ j){
        cout << " " << setw(3) << j;
    }
    cout << endl << " __|";
    for(value_type j{}; j != _length; ++ j){
        cout << "_ _ ";
    }
    cout << endl;
    for(value_type r{} ; r != _length; ++r){
        cout << setw(3) << r << "|";
        for(value_type c{} ; c != _length ; ++ c){
            value_type i_site = r * _length + c;
            cout << setw(3) << _sites[i_site].get_groupID() << ' ';
        }
        cout << endl;
    }
}


/**
 *
 */
void SqLattice_v2::view_sites_by_relative_index(){
    std::cout << "Relative index : line " << __LINE__ << endl;
    cout << "Format: \"id(x,y)\"" << endl;
    cout << "   |";
    for(value_type j{}; j != _length; ++ j){
        cout << setw(4) << j << "         |";
    }
    cout << endl;
    print_h_barrier(_length, "___|__", "___________|__");
    for(value_type r{} ; r != _length; ++r){
        cout << setw(3) << r << "|";
        for(value_type c{} ; c != _length ; ++ c){
            value_type i_site = r * _length + c;
            if(_sites[i_site].get_groupID() == -1){
                // left blank
                cout << setw(4) << _sites[i_site].get_groupID() << "         |";
                continue;
            }
            cout << setw(4) << _sites[i_site].get_groupID() << _sites[i_site].relativeIndex() << "|";
        }
        cout << endl;
        print_h_barrier(_length, "___|__", "___________|__");
    }
//    print_h_barrier(_length, "___|__", "___________|_");
}

/**
 * View bonds in the lattice by relative index.
 * format : id(relative_index)
 */
void SqLattice_v2::view_bonds_by_relative_index() {

    std::cout << "Bonds by id : line " << __LINE__ << endl;
//    // printing indices for columns
//    std::cout << "    | ";
//    for(value_type i{}; i != _length; ++i){
//        std::cout << "   " << setw(4) << i << "         | ";
//    }
//    std::cout << std::endl;
//
//    // pringing H,V label
//
//    print_h_barrier(_length, "    |  ", "V           H  |  ");
//    print_h_barrier(_length, "____|__", "_______________|__");
//
//    // for each row there will be two columns
//    for(value_type i{}; i != _length; ++i){
//        std::cout << i << ' ';
//        std::cout << "H |";
//        for(value_type j1{}; j1 != _length; ++j1){
//            std::cout << "     " << std::setw(3) << _h_bonds[i][j1].get_groupID()
//                      << _h_bonds[i][j1].relativeIndex() << "|";
//        }
//        std::cout << std::endl;
//        print_h_barrier(_length, "    |  ", "               |  "); // just for better viewing
//        std::cout << "  " << "V |";
//        for(value_type j2{}; j2 != _length; ++j2){
//            std::cout << std::setw(3) << _v_bonds[i][j2].get_groupID()
//                      << _v_bonds[i][j2].relativeIndex() << "     |";
//        }
//        std::cout << std::endl;
//
//        // printing horizontal separator
//        print_h_barrier(_length, "____|__", "_______________|__");
//    }
//    std::cout << std::endl;
}


/**
 * View bonds in the lattice by relative index. id of the site is showed
 * format : id for site or id(relative_index) for bond
 */
void SqLattice_v2::view_bonds_by_relative_index_v2() {

    std::cout << "Bonds by id : line " << __LINE__ << endl;
//    cout << "site id -1 means isolated site and 0 means connected site in bond percolation(definition)" << endl;
//    print_h_barrier(15, "_", "___", "_\n");
//    cout << "|(site id) (horizontal bond id(relative index))|" << endl;
//    cout << "|(vertical bond id(relative index))            |" << endl;
//    print_h_barrier(15, "-", "---", "-\n");
//    // printing indices for columns
//    std::cout << "    | ";
//    for(value_type i{}; i != _length; ++i){
//        std::cout << i << "                  | ";
//    }
//    std::cout << std::endl;
//
//    // pringing H,V label
//
//    print_h_barrier(_length, "    |  ", " V            H   |  ");
//    print_h_barrier(_length, "____|__", "__________________|__");
//
//    // for each row there will be two columns
//    for(value_type i{}; i != _length; ++i){
//        std::cout << i << ' ';
//        std::cout << "H |";
//        for(value_type j1{}; j1 != _length; ++j1){
//            std::cout << std::setw(3) << _sites[i][j1].get_groupID() ;
//            std::cout << "     " << std::setw(3) << _h_bonds[i][j1].get_groupID()
//                      << _h_bonds[i][j1].relativeIndex() << "|";
//        }
//        std::cout << std::endl;
//        print_h_barrier(_length, "    |  ", "                  |  "); // just to see a better view
//        std::cout << "  " << "V |";
//        for(value_type j2{}; j2 != _length; ++j2){
//            std::cout << std::setw(3) << _v_bonds[i][j2].get_groupID()
//                      << _v_bonds[i][j2].relativeIndex() << "        |";
//        }
//        std::cout << std::endl;
//
//        // printing horizontal separator
//        print_h_barrier(_length, "____|__", "__________________|__");
//    }
//
//    std::cout << std::endl;
}


/**
 * View bonds in the lattice by relative index. id of the site is showed
 * format : id(relative index) for site and only id for bond
 */
void SqLattice_v2::view_bonds_by_relative_index_v3() {

    std::cout << "Bonds by id : line " << __LINE__ << endl;
//    //cout << "site id -1 means isolated site and 0 means connected site in bond percolation(definition)" << endl;
//    print_h_barrier(15, "_", "___", "_\n");
//    cout << "|(site id) (horizontal bond id(relative index))|" << endl;
//    cout << "|(vertical bond id(relative index))            |" << endl;
//    print_h_barrier(15, "-", "---", "-\n");
//    // printing indices for columns
//    std::cout << "    | ";
//    for(value_type i{}; i != _length; ++i){
//        std::cout << i << "                | ";
//    }
//    std::cout << std::endl;
//
//    // pringing H,V label
//
//    print_h_barrier(_length, "    |  ", " V           H  |  ");
//    print_h_barrier(_length, "____|__", "________________|__");
//
//    // for each row there will be two columns
//    for(value_type i{}; i != _length; ++i){
//        std::cout << i << ' ';
//        std::cout << "H |";
//        for(value_type j1{}; j1 != _length; ++j1){
//            std::cout << std::setw(3) << _sites[i][j1].get_groupID() << _sites[i][j1].relativeIndex() ;
//            std::cout << "   " << std::setw(3) << _h_bonds[i][j1].get_groupID() << "|";
//        }
//        std::cout << std::endl;
//        print_h_barrier(_length, "    |  ", "                |  "); // just to see a better view
//        std::cout << "  " << "V |";
//        for(value_type j2{}; j2 != _length; ++j2){
//            std::cout << std::setw(3) << _v_bonds[i][j2].get_groupID() << "               |";
//        }
//        std::cout << std::endl;
//
//        // printing horizontal separator
//        print_h_barrier(_length, "____|__", "________________|__");
//    }
//
//    std::cout << std::endl;
}

/**
 * View bonds in the lattice by relative index. id of the site is showed
 * format : id(relative index) for site and only id for bond
 * if any site is isolated relative index is not shown
 */
void SqLattice_v2::view_bonds_by_relative_index_v4() {

    std::cout << "Bonds by id : line " << __LINE__ << endl;
//    //cout << "site id -1 means isolated site and 0 means connected site in bond percolation(definition)" << endl;
//    print_h_barrier(15, "_", "___", "_\n");
//    cout << "|(site id) (horizontal bond id(relative index))|" << endl;
//    cout << "|(vertical bond id(relative index))            |" << endl;
//    print_h_barrier(15, "-", "---", "-\n");
//    // printing indices for columns
//    std::cout << "    | ";
//    for(value_type i{}; i != _length; ++i){
//        std::cout << i << "                | ";
//    }
//    std::cout << std::endl;
//
//    // pringing H,V label
//
//    print_h_barrier(_length, "    |  ", " V           H  |  ");
//    print_h_barrier(_length, "____|__", "________________|__");
//
//    // for each row there will be two columns
//    for(value_type i{}; i != _length; ++i){
//        std::cout << i << ' ';
//        std::cout << "H |";
//        for(value_type j1{}; j1 != _length; ++j1){
//            int id = _sites[i][j1].get_groupID();
//            std::cout << std::setw(3) << id  ;
//            if(id != -1){
//                cout << _sites[i][j1].relativeIndex();
//            }else{
//                cout << "(-,-)    ";
//            }
//            std::cout << "   " << std::setw(3) << _h_bonds[i][j1].get_groupID() << "|";
//        }
//        std::cout << std::endl;
//        print_h_barrier(_length, "    |  ", "                |  "); // just to see a better view
//        std::cout << "  " << "V |";
//        for(value_type j2{}; j2 != _length; ++j2){
//            std::cout << std::setw(3) << _v_bonds[i][j2].get_groupID() << "               |";
//        }
//        std::cout << std::endl;
//
//        // printing horizontal separator
//        print_h_barrier(_length, "____|__", "________________|__");
//    }
//
//    std::cout << std::endl;
}



/**
 * View lattice (sites and bonds) by relative index.
 * format : id(relative_index)
 */
void SqLattice_v2::view_by_relative_index() {

    std::cout << "Bonds by id : line " << __LINE__ << endl;
//    // printing indices for columns
//    std::cout << "    | ";
//    for(value_type i{}; i != _length; ++i){
//        std::cout << i << "                         | ";
//    }
//    std::cout << std::endl;
//
//    // pringing H,V label
//
//    print_h_barrier(_length, "    |  ", " V                   H   |  ");
//    print_h_barrier(_length, "____|__", "_________________________|__");
//
//    // for each row there will be two columns
//    for(value_type i{}; i != _length; ++i){
//        std::cout << i << ' ';
//        std::cout << "H |";
//        for(value_type j1{}; j1 != _length; ++j1){
//            std::cout << std::setw(3) << _sites[i][j1].get_groupID() << _sites[i][j1].relativeIndex();
//            std::cout << "   " << std::setw(3) << _h_bonds[i][j1].get_groupID()
//                      << _h_bonds[i][j1].relativeIndex() << "|";
//        }
//        std::cout << std::endl;
//        print_h_barrier(_length, "    |  ", "                         |  "); // just to see a better view
//        std::cout << "  " << "V |";
//        for(value_type j2{}; j2 != _length; ++j2){
//            std::cout << std::setw(3) << _v_bonds[i][j2].get_groupID()
//                      << _v_bonds[i][j2].relativeIndex() << "               |";
//        }
//        std::cout << std::endl;
//
//        // printing horizontal separator
//        print_h_barrier(_length, "____|__", "_________________________|__");
//    }
//
////        view_h_bonds_extended();
////        view_v_bonds_extended();
//    std::cout << std::endl;

}

/**
 * View lattice (sites and bonds) by relative index.
 * format : id(relative_index)
 */
void SqLattice_v2::view() {

    std::cout << "Bonds by id : line " << __LINE__ << endl;
    cout << "Structure " << endl;
    print_h_barrier(10, "_", "___", "_\n");
    cout << "|(site id) (horizontal bond id)|" << endl;
    cout << "|(vertical bond id)            |" << endl;
    print_h_barrier(10, "-", "---", "-\n");
    // printing indices for columns
    std::cout << "    | ";
    for(value_type i{}; i != _length; ++i){
        std::cout << i << "      | ";
    }
    std::cout << std::endl;

    // pringing H,V label

    print_h_barrier(_length, "    |  ", "V    H|  ");
    print_h_barrier(_length, "____|__", "______|__");

    // for each row there will be two columns
    value_type i_site{};
    for(value_type i{}; i != _length; ++i){
        std::cout << i << ' ';
        std::cout << "H |";
        for(value_type j1{}; j1 != _length; ++j1){
            i_site = i * _length + j1;
            std::cout << std::setw(3) << _sites[i_site].get_groupID() ;
            std::cout << "  " << std::setw(3) << _bonds[i_site].get_groupID() << "|";
        }
        std::cout << std::endl;
        print_h_barrier(_length, "    |  ", "      |  "); // just to see a better view
        std::cout << "  " << "V |";
        for(value_type j2{}; j2 != _length; ++j2){
            i_site = i * _length + j2;
            std::cout << std::setw(3) << _bonds[i_site].get_groupID() << "     |";
        }
        std::cout << std::endl;

        // printing horizontal separator
        print_h_barrier(_length, "____|__", "______|__");
    }

//        view_h_bonds_extended();
//        view_v_bonds_extended();
    std::cout << std::endl;

}


/**
 *
 */
void SqLattice_v2::view_h_bonds()
{
    std::cout << "view horizontal bonds" << std::endl;
    std::cout << '{';
    for(value_type r{} ; r != _length ; ++r) {
        if(r!=0) std::cout << "  ";
        else std::cout << '{';
        for (value_type c{}; c != _length; ++c) {
            value_type i_site = r * _length + c;
            std::cout << _bonds[i_site] ;
            if(c != _length-1)
                std::cout << ',';
        }
        std::cout << '}';
        if(r != _length-1)
            std::cout << std::endl;
    }
    std::cout << '}';
    std::cout << std::endl;
}

/**
 *
 */
void SqLattice_v2::view_v_bonds()
{
    std::cout << "view vertical bonds" << std::endl;
    std::cout << '{';
    for(value_type r{} ; r != _length ; ++r) {
        if(r!=0) std::cout << "  ";
        else std::cout << '{';
        for (value_type c{}; c != _length; ++c) {
            value_type i_site = r * _length + c;
            std::cout << _bonds[i_site + _length_squared] ;
            if(c != _length-1)
                std::cout << ',';
        }
        std::cout << '}';
        if(r != _length-1)
            std::cout << std::endl;
    }
    std::cout << '}';
    std::cout << std::endl;
}


/**
 *
 */
void SqLattice_v2::view_h_bonds_extended(){
    std::cout << "view horizontal bonds" << std::endl;
    std::cout << '{';
    for(value_type r{} ; r != _length ; ++r) {
        if(r!=0) std::cout << "  ";
        else std::cout << '{';
        for (value_type c{}; c != _length; ++c) {
            value_type i_site = r * _length + c;
            std::cout << "(" << _bonds[i_site].get_groupID() << ":" << _bonds[i_site] << ")" ;
            if(c != _length-1)
                std::cout << ',';
        }
        std::cout << '}';
        if(r != _length-1)
            std::cout << std::endl;
    }
    std::cout << '}';
    std::cout << std::endl;
}

/**
 *
 */
void SqLattice_v2::view_v_bonds_extended(){
    std::cout << "view vertical bonds" << std::endl;
    std::cout << '{';
    for(value_type r{} ; r != _length ; ++r) {
        if(r!=0) std::cout << "  ";
        else std::cout << '{';
        for (value_type c{}; c != _length; ++c) {
            value_type i_site = r * _length + c;
            std::cout << "(" << _bonds[i_site + _length_squared].get_groupID() << ":" << _bonds[i_site] << ")" ;
            if(c != _length-1)
                std::cout << ',';
        }
        std::cout << '}';
        if(r != _length-1)
            std::cout << std::endl;
    }
    std::cout << '}';
    std::cout << std::endl;
}


/**
 *
 */
void SqLattice_v2::view_bonds_by_id(){
    std::cout << "Bonds by id : line " << __LINE__ << endl;
    cout << "Structure " << endl;
//    print_h_barrier(8, "__", "___", "_\n");
//    cout << "|      (horizontal bond id)|" << endl;
//    cout << "|(vertical bond id)        |" << endl;
//    print_h_barrier(8, "--", "---", "-\n");
//    // printing indices for columns
//    std::cout << "    | ";
//    for(value_type i{}; i != _length; ++i){
//        std::cout << i << "    | ";
//    }
//    std::cout << std::endl;
//
//    // pringing H,V label
//
//    print_h_barrier(_length, "    |  ", "V  H|  ");
//    print_h_barrier(_length, "____|__", "____|__");
//
//    // for each row there will be two columns
//    for(value_type i{}; i != _length; ++i){
//        std::cout << i << ' ';
//        std::cout << "H |";
//        for(value_type j1{}; j1 != _length; ++j1){
//            std::cout << "   " << std::setw(3) << _h_bonds[i][j1].get_groupID() << "|";
//        }
//        std::cout << std::endl;
//        std::cout << "  " << "V |";
//        for(value_type j2{}; j2 != _length; ++j2){
//            std::cout << std::setw(3) << _v_bonds[i][j2].get_groupID() << "   |";
//        }
//        std::cout << std::endl;
//
//        // printing horizontal separator
//        print_h_barrier(_length, "____|__", "____|__");
//    }
//
////        view_h_bonds_extended();
////        view_v_bonds_extended();
//    std::cout << std::endl;
}


Site& SqLattice_v2::getSite(Index index) {
    size_t i = translateSiteTo1D(index);
    return _sites[i];
}

//Site&& SqLattice_v2::getSiteR(Index index) {
//    Site a = _sites[index.x_][index.y_];
//    return std::move(a);
//}

const Site& SqLattice_v2::getSite(Index index) const {
    size_t i = translateSiteTo1D(index);
    return _sites[i];
}

void SqLattice_v2::setGroupID(Index index, int group_id){
    size_t i = translateSiteTo1D(index);
    _sites[i].set_groupID(group_id);
}

void SqLattice_v2::setGroupID(BondIndex index, int group_id){
    size_t i = translateBondTo1D(index);
    _bonds[i].set_groupID(group_id);

}

int SqLattice_v2::getGroupID(Index index){
    size_t i = translateSiteTo1D(index);
    return _sites[i].get_groupID();
}

int SqLattice_v2::getGroupID(BondIndex index){
    size_t i = translateBondTo1D(index);
    return _bonds[i].get_groupID();
}

void SqLattice_v2::setGroupIDSite(size_t index, int group_id){
    _sites[index].set_groupID(group_id);
}

void SqLattice_v2::setGroupIDBond(size_t index, int group_id){
    _bonds[index].set_groupID(group_id);
}

int SqLattice_v2::getGroupIDSite(size_t index){
    return _sites[index].get_groupID();
}

int SqLattice_v2::getGroupIDBond(size_t index){
    return _bonds[index].get_groupID();
}


Bond_v2& SqLattice_v2::getBond(BondIndex index) {
    size_t i = translateBondTo1D(index);
    return _bonds[i];

}

void SqLattice_v2::reset(bool reset_all) {
    if(reset_all){
        reset_sites();
        reset_bonds();
        return;
    }
    // setting all group id to -1
    if(_site_resetting_flag) {
        reset_sites();
    }
//    cout << "Bond resetting is disabled : line " << __LINE__ << endl;
    if(_bond_resetting_flag) {
        reset_bonds();
    }

}

/**
 *
 */
void SqLattice_v2::reset_bonds() {
    for(value_type i{}; i != _bonds.size(); ++i){
        _bonds[i].deactivate();
    }
}

/**
 *
 */
void SqLattice_v2::reset_sites() {
    for(value_type i{}; i != _sites.size(); ++i){
        _sites[i].deactivate();
    }
}



/******************************************************************************
 * Get Neighbor from given index
 ******************************************************************************/
/**
 * Periodic case only.
 * Each site has four neighbor sites.
 * @param site
 * @return
 */
std::vector<Index> SqLattice_v2::get_neighbor_sites(Index site){
    std::vector<Index> sites(4);
    sites[0] = {(site.row_ + 1) % _length, site.column_};
    sites[1] = {(site.row_ - 1 + _length) % _length, site.column_};
    sites[2] = {site.row_, (site.column_ + 1) % _length};
    sites[3] = {site.row_, (site.column_ - 1 + _length) % _length};
    return sites;
}

/**
 * Periodic case only.
 * Each bond has six neibhbor bonds.
 * @param site
 * @return
 */
std::vector<BondIndex> SqLattice_v2::get_neighbor_bonds(BondIndex bond) {
    value_type next_column = (bond.column_ + 1) % _length;
    value_type prev_column = (bond.column_ - 1 + _length) % _length;
    value_type prev_row = (bond.row_ - 1 + _length) % _length;
    value_type next_row = (bond.row_ + 1) % _length;

    vector<BondIndex> bonds(6);

    // horizontal bond case
    if (bond.horizontal()) {
        // increase column index for the right neighbor

        // left end of bond
        bonds[0] = {BondType::Vertical, bond.row_, bond.column_};
        bonds[1] = {BondType::Vertical, prev_row, bond.column_};
        bonds[2] = {BondType::Horizontal, bond.row_, prev_column};

        // right end bond
        bonds[3] = {BondType::Vertical, prev_row, next_column};
        bonds[4] = {BondType::Vertical, bond.row_, next_column};
        bonds[5] = {BondType::Horizontal, bond.row_, next_column};

    }
        // vertical bond case
    else if (bond.vertical()) {
        // increase row index

        // top end of bond
        bonds[0] = {BondType::Horizontal, bond.row_, bond.column_};
        bonds[1] = {BondType::Horizontal, bond.row_, prev_column};
        bonds[2] = {BondType::Vertical, prev_row, bond.column_};

        // bottom end of bond
        bonds[3] = {BondType::Horizontal, next_row, bond.column_};
        bonds[4] = {BondType::Horizontal, next_row, prev_column};
        bonds[5] = {BondType::Vertical, next_row, bond.column_};

    }


    return bonds;
}

std::vector<Index> SqLattice_v2::get_neighbor_sites(BondIndex bond) {
    value_type r = bond.row_;
    value_type c = bond.column_;
    vector<Index> sites(2);
    sites[0] = {r, c};
    if(bond.horizontal()){
        sites[1] = {r, (c+1) % _length};
    }else{
        sites[1] = {(r+1) % _length, c};
    }
    return sites;
}

///************************************************************
// * Static methods
// */
//std::vector<Index> SqLattice_v2::get_neighbor_sites(size_t length, Index site){
//    std::vector<Index> sites(4);
//    sites[0] = {(site.row_ + 1) % length, site.column_};
//    sites[1] = {(site.row_ - 1 + length) % length, site.column_};
//    sites[2] = {site.row_, (site.column_ + 1) % length};
//    sites[3] = {site.row_, (site.column_ - 1 + length) % length};
//    return sites;
//}
//
///**
// * Periodic case only.
// * Each bond has six neibhbor bonds.
// * @param site
// * @return
// */
//std::vector<BondIndex> SqLattice_v2::get_neighbor_bond_indices(size_t length, BondIndex bond) {
//    value_type next_column = (bond.column_ + 1) % length;
//    value_type prev_column = (bond.column_ - 1 + length) % length;
//    value_type prev_row = (bond.row_ - 1 + length) % length;
//    value_type next_row = (bond.row_ + 1) % length;
//
//    vector<BondIndex> bonds(6);
//
//    // horizontal bond case
//    if (bond.horizontal()) {
//        // increase column index for the right neighbor
//
//        // left end of bond
//        bonds[0] = {BondType::Vertical, bond.row_, bond.column_};
//        bonds[1] = {BondType::Vertical, prev_row, bond.column_};
//        bonds[2] = {BondType::Horizontal, bond.row_, prev_column};
//
//        // right end bond
//        bonds[3] = {BondType::Vertical, prev_row, next_column};
//        bonds[4] = {BondType::Vertical, bond.row_, next_column};
//        bonds[5] = {BondType::Horizontal, bond.row_, next_column};
//
//    }
//        // vertical bond case
//    else if (bond.vertical()) {
//        // increase row index
//
//        // top end of bond
//        bonds[0] = {BondType::Horizontal, bond.row_, bond.column_};
//        bonds[1] = {BondType::Horizontal, bond.row_, prev_column};
//        bonds[2] = {BondType::Vertical, prev_row, bond.column_};
//
//        // bottom end of bond
//        bonds[3] = {BondType::Horizontal, next_row, bond.column_};
//        bonds[4] = {BondType::Horizontal, next_row, prev_column};
//        bonds[5] = {BondType::Vertical, next_row, bond.column_};
//
//    }
//
//
//    return bonds;
//}
//
//std::vector<Index> SqLattice_v2::get_neighbor_indices(size_t length, BondIndex bond) {
//    value_type r = bond.row_;
//    value_type c = bond.column_;
//    vector<Index> sites(2);
//    sites[0] = {r, c};
//    if(bond.horizontal()){
//        sites[1] = {r, (c+1) % length};
//    }else{
//        sites[1] = {(r+1) % length, c};
//    }
//    return sites;
//}

std::vector<Index> SqLattice_v2::getSites() {
    std::vector<Index> indices;
    for(size_t i{}; i < _sites.size(); ++i){
        indices.push_back(_sites[i]._id);
    }
    return indices;
}

std::vector<BondIndex> SqLattice_v2::getBonds() {
    std::vector<BondIndex> indices;
    // number of horizontal bonds and number of vertical bonds are equal
    for(size_t i{}; i < _bonds.size(); ++i){
        indices.push_back(_bonds[i]._id);
    }
    return indices;
}

std::vector<BondIndex> SqLattice_v2::get_neighbor_bonds(Index site) {
    value_type next_column = (site.column_ + 1) % _length;
    value_type prev_column = (site.column_ - 1 + _length) % _length;
    value_type prev_row = (site.row_ - 1 + _length) % _length;
    value_type next_row = (site.row_ + 1) % _length;

    vector<BondIndex> bond_neighbor(4);
    bond_neighbor[0] = {BondType::Horizontal, site.row_, site.column_};
    bond_neighbor[1] = {BondType::Horizontal, site.row_, prev_column};
    bond_neighbor[2] = {BondType::Vertical,    site.row_, site.column_};
    bond_neighbor[3] = {BondType::Vertical, prev_row, site.column_};

    return bond_neighbor;
}


/**
 * Take a bond index only if the corresponding site is active
 * takes longer? time than version 1?, i.e.,  connection()
 * @param site
 * @param site_neighbor
 * @param bond_neighbor
 */
void SqLattice_v2::get_neighbors(Index site, vector<Index> &site_neighbor, vector<BondIndex> &bond_neighbor)
{
    site_neighbor.clear();
    bond_neighbor.clear();

    value_type prev_column  = (site.column_ + _length - 1) % _length;
    value_type prev_row     = (site.row_ + _length - 1) % _length;
    value_type next_row     = (site.row_ + 1) % _length;
    value_type next_column  = (site.column_ + 1) % _length;

    // 1 level inside the lattice
    // not in any the boundary
    site_neighbor.resize(4);
    site_neighbor[0] = {site.row_, next_column};
    site_neighbor[1] = {site.row_, prev_column};
    site_neighbor[2] = {next_row, site.column_};
    site_neighbor[3] = {prev_row, site.column_};

    bond_neighbor.reserve(4);

//    if(!_lattice.getSite(site_neighbor[0]).isActive()) {
//        bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
//    }
//    if(!_lattice.getSite(site_neighbor[1]).isActive()){
//        bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
//    }
//    if(!_lattice.getSite(site_neighbor[2]).isActive()){
//        bond_neighbor.push_back({BondType::Vertical,    site.row_, site.column_});
//    }
//    if(!_lattice.getSite(site_neighbor[3]).isActive()) {
//        bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
//    }
    bond_neighbor.resize(4);
    bond_neighbor[0] = {BondType::Horizontal, site.row_, site.column_};
    bond_neighbor[1] = {BondType::Horizontal, site.row_, prev_column};
    bond_neighbor[2] = {BondType::Vertical,    site.row_, site.column_};
    bond_neighbor[3] = {BondType::Vertical, prev_row, site.column_};

}

void SqLattice_v2::setRelativeIndex(Index index, IndexRelative ir) {
    size_t i = translateSiteTo1D(index);
    _sites[i].relativeIndex(ir);
}

IndexRelative SqLattice_v2::getRelativeIndex(Index index) {
    size_t i = translateSiteTo1D(index);
    return _sites[i].relativeIndex();
}

void SqLattice_v2::setRelativeIndex(size_t index, IndexRelative ir) {
    _sites[index].relativeIndex(ir);
}

IndexRelative SqLattice_v2::getRelativeIndex(size_t index) {
    return _sites[index].relativeIndex();
}

/**
 *  (b) 1D index to 2D BondIndex and vice versa.
 *          (t,r,c) = t*L^2 + r*L + c = index
 *      here, r -> row index
 *            c -> column index
 *            t -> bond type. 0 if horizontal and 1 if vertical
 * @param bondIndex
 * @return
 */
uint SqLattice_v2::translateBondTo1D(BondIndex bondIndex) const{
    size_t t = (bondIndex.bondType == BondType ::Horizontal) ? 0 : 1;
    size_t index  = _length_squared * t + _length * bondIndex.row_ + bondIndex.column_;

    return index;
}

BondIndex SqLattice_v2::translate1DToBond(uint index) const{
    size_t t = index / _length_squared; // determine horizontal or vertical
    size_t half = index - t*_length_squared;
    size_t r = half / _length;
    size_t c = half % _length;
    auto bt = (t == 0)? BondType::Horizontal : BondType::Vertical;
    return BondIndex(bt, r, c);
}

/**
 *  (a) 1D index to 2D and vice versa for Index.
 *          (r,c) = r*L + c = index
 *
 *      here, r -> row index
 *            c -> column index
 * @param index
 * @return
 */
uint SqLattice_v2::translateSiteTo1D(Index index) const{
    size_t i  = _length * index.row_ + index.column_;
    return i;
}

Index SqLattice_v2::translate1DToSite(uint index) const{
    size_t c = index % _length;
    size_t r = index / _length;
    return Index(r, c);
}


void SqLattice_v2::init(bool activate_bonds, bool activate_sites, bool bond_reset, bool site_reset) {
    _bond_resetting_flag = bond_reset;
    _site_resetting_flag = site_reset;
    for(size_t i{}; i < _length_squared; ++i){
        _sites[i] = Site(translate1DToSite(i), _length);
        _bonds[i] = Bond_v2(translate1DToBond(i)); // horizontal
        _bonds[i+_length_squared] = Bond_v2(translate1DToBond(i + _length_squared)); // vertical
        if(activate_bonds) {
            _bonds[i].activate(); // horizontal bond
            _bonds[i+_length_squared].activate(); // vertical bond
        }
        if(activate_sites){
            _sites[i].activate();
        }
    }
}

std::vector<uint> SqLattice_v2::getSiteIndices() {
    std::vector<uint> a(_sites.size());
    for(uint i{}; i < _sites.size(); ++i){a[i]=i;}
    return a;
}

std::vector<uint> SqLattice_v2::getBondIndices() {
    std::vector<uint> a(_bonds.size());
    for(uint i{}; i < _bonds.size(); ++i){a[i]=i;}
    return a;
}

void SqLattice_v2::view_bonds() {
    std::cout << "view all bonds" << std::endl;
    value_type  step = _length_squared - 1;
    std::cout << '{';
    for(value_type r{} ; r != 2*_length_squared ; ++r) {
        if(r % _length_squared == 0) cout << endl;
        if(r!=0) std::cout << "  ";
        else std::cout << '{';
        std::cout << _bonds[r] << ",";
        if((r+1) % _length == 0)  std::cout << '}' << endl;

    }
    std::cout << '}';
    std::cout << std::endl;
}






