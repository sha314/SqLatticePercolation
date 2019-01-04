//
// Created by shahnoor on 10/2/2017.
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
SqLattice::SqLattice(
        value_type length,
        bool activate_bonds, bool activate_sites,
        bool bond_reset, bool site_reset)
        : _length{length}, _bond_resetting_flag{bond_reset}, _site_resetting_flag{site_reset}
{
    cout << "Constructing Lattice object : line " << __LINE__ << endl;
    _sites = std::vector<std::vector<Site>>(_length);
    _h_bonds = std::vector<std::vector<Bond>>(_length);
    _v_bonds = std::vector<std::vector<Bond>>(_length);
    if(!activate_bonds && !activate_sites) {    // both are deactivated by default
        for (value_type i{}; i != _length; ++i) {
            _sites[i] = std::vector<Site>(_length);
            _h_bonds[i] = std::vector<Bond>(_length);
            _v_bonds[i] = std::vector<Bond>(_length);
            for (value_type j{}; j != _length; ++j) {
                _sites[i][j] = Site(Index(i, j), _length);
                _h_bonds[i][j] = {Index(i, j), Index(i, (j + 1) % _length), _length};
                _v_bonds[i][j] = {Index(i, j), Index((i + 1) % _length, j), _length};
            }
        }
    }
    else if(activate_bonds && !activate_sites) {    // all bonds are activated by default
        for (value_type i{}; i != _length; ++i) {
            _sites[i] = std::vector<Site>(_length);
            _h_bonds[i] = std::vector<Bond>(_length);
            _v_bonds[i] = std::vector<Bond>(_length);
            for (value_type j{}; j != _length; ++j) {
                _sites[i][j] = Site(Index(i, j), _length);
                _h_bonds[i][j] = {Index(i, j), Index(i, (j + 1) % _length), _length};
                _v_bonds[i][j] = {Index(i, j), Index((i + 1) % _length, j), _length};
                _h_bonds[i][j].activate();
                _v_bonds[i][j].activate();
            }
        }
    }
    else if(!activate_bonds && activate_sites) {    // all sites are activated by default
        for (value_type i{}; i != _length; ++i) {
            _sites[i] = std::vector<Site>(_length);
            _h_bonds[i] = std::vector<Bond>(_length);
            _v_bonds[i] = std::vector<Bond>(_length);
            for (value_type j{}; j != _length; ++j) {
                _sites[i][j] = Site(Index(i, j), _length);
                _sites[i][j].activate();
                _h_bonds[i][j] = {Index(i, j), Index(i, (j + 1) % _length), _length};
                _v_bonds[i][j] = {Index(i, j), Index((i + 1) % _length, j), _length};
            }
        }
    }
    else {
        for (value_type i{}; i != _length; ++i) {   // all bonds and sites are activated by default
            _sites[i] = std::vector<Site>(_length);
            _h_bonds[i] = std::vector<Bond>(_length);
            _v_bonds[i] = std::vector<Bond>(_length);
            for (value_type j{}; j != _length; ++j) {
                _sites[i][j] = Site(Index(i, j), _length);
                _sites[i][j].activate();
                _h_bonds[i][j] = {Index(i, j), Index(i, (j + 1) % _length), _length};
                _v_bonds[i][j] = {Index(i, j), Index((i + 1) % _length, j), _length};
                _h_bonds[i][j].activate();
                _v_bonds[i][j].activate();
            }
        }
    }

}

/********************************************************************
 * Activation and Deactivation
 ********************************************************************/
void SqLattice::activateAllSite()
{
    for(value_type i{} ; i != _length ; ++i) {
        for (value_type j{}; j != _length; ++j) {
            _sites[i][j].activate();
        }
    }
}

/**
 *
 */
void SqLattice::activateAllBond()
{
    for(value_type i{} ; i != _length ; ++i) {
        for (value_type j{}; j != _length; ++j) {
            _h_bonds[i][j].activate();
            _v_bonds[i][j].activate();
        }
    }
}

void SqLattice::activate_site(Index index) {
//    cout << "activating site " << index << endl;
    _sites[index.row_][index.column_].activate();
}


void SqLattice::activateBond(BondIndex bond) {
    // check if the bond is vertical or horizontal
    // then call appropriate function to activate _h_bond or _v_bond
    if(bond.horizontal()){ // horizontal
        if(_h_bonds[bond.row_][bond.column_].isActive()){
            cout << "Bond is already activated : line " << __LINE__ << endl;
        }
        _h_bonds[bond.row_][bond.column_].activate();
    }
    else if(bond.vertical()) // vertical
    {
        if(_v_bonds[bond.row_][bond.column_].isActive()){
            cout << "Bond is already activated : line " << __LINE__ << endl;
        }
        _v_bonds[bond.row_][bond.column_].activate();
    }
    else{
        cout << bond << " is not a valid bond : line " << __LINE__ << endl;
    }

}


void SqLattice::deactivate_site(Index index){
    _sites[index.row_][index.column_].deactivate();
}


void SqLattice::deactivate_bond(Bond bond) {
    // check if the bond is vertical or horizontal
    // then call appropriate function to activate _h_bond or _v_bond
    if(bond.isHorizontal()){
        if(_h_bonds[bond.id().row_][bond.id().column_].isActive()){
            cout << "Bond is already activated : line " << __LINE__ << endl;
        }
        _h_bonds[bond.id().row_][bond.id().column_].deactivate();
    }
    else if(bond.isVertical())
    {
        if(_v_bonds[bond.id().row_][bond.id().column_].isActive()){
            cout << "Bond is already activated : line " << __LINE__ << endl;
        }
        _v_bonds[bond.id().row_][bond.id().column_].deactivate();
    }
    else{
        bond.activate();
        cout << bond << " is not a valid bond : line " << __LINE__ << endl;
    }
}


/*********************************************************
 * Viewing methods
 *********************************************************/
/**
 *  View the sites of the lattice
 *  place (*) if the site is not active
 */
void SqLattice::view_sites()
{
    std::cout << "view sites" << std::endl;
    std::cout << '{';
    for(value_type i{} ; i != _length ; ++i) {
        if(i!=0) std::cout << "  ";
        else std::cout << '{';
        for (value_type j{}; j != _length; ++j) {
            if(_sites[i][j].isActive()){
                std::cout << _sites[i][j] ;
            }
            else{
                std::cout << "(*)";
            }
            if(j != _length-1)
                std::cout << ',';
        }
        std::cout << '}';
        if(i != _length-1)
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
void SqLattice::view_sites_extended()
{
    std::cout << "view sites" << std::endl;
    std::cout << '{';
    for(value_type i{} ; i != _length ; ++i) {
        if(i!=0) std::cout << "  ";
        else std::cout << '{';
        for (value_type j{}; j != _length; ++j) {
            std::cout << std::setw(3) << _sites[i][j].get_groupID() << ":";

            if(_sites[i][j].isActive()) {
                cout << '(' << std::setw(2) << _sites[i][j]._id.row_ << ','
                     << std::setw(2) << _sites[i][j]._id.column_ << ')';
            }
            else{
                cout << std::setw(7) << "(*)";
            }

            if(j != _length-1)
                std::cout << ',';
        }
        std::cout << '}';
        if(i != _length-1)
            std::cout << std::endl;
    }
    std::cout << '}';
    std::cout << std::endl;
}



/**
 * Displays group ids of sites in a matrix form
 */
void SqLattice::view_sites_by_id() {
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
    for(value_type i{} ; i != _length; ++i){
        cout << setw(3) << i << "|";
        for(value_type j{} ; j != _length ; ++ j){
            cout << setw(3) << _sites[i][j].get_groupID() << ' ';
        }
        cout << endl;
    }
}


/**
 *
 */
void SqLattice::view_sites_by_relative_index(){
    std::cout << "Relative index : line " << __LINE__ << endl;
    cout << "Format: \"id(x,y)\"" << endl;
    cout << "   |";
    for(value_type j{}; j != _length; ++ j){
        cout << setw(4) << j << "         |";
    }
    cout << endl;
    print_h_barrier(_length, "___|__", "___________|__");
    for(value_type i{} ; i != _length; ++i){
        cout << setw(3) << i << "|";
        for(value_type j{} ; j != _length ; ++ j){
            if(_sites[i][j].get_groupID() == -1){
                // left blank
                cout << setw(4) << _sites[i][j].get_groupID() << "         |";
                continue;
            }
            cout << setw(4) << _sites[i][j].get_groupID() << _sites[i][j].relativeIndex() << "|";
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
void SqLattice::view_bonds_by_relative_index() {

    std::cout << "Bonds by id : line " << __LINE__ << endl;
    // printing indices for columns
    std::cout << "    | ";
    for(value_type i{}; i != _length; ++i){
        std::cout << "   " << setw(4) << i << "         | ";
    }
    std::cout << std::endl;

    // pringing H,V label

    print_h_barrier(_length, "    |  ", "V           H  |  ");
    print_h_barrier(_length, "____|__", "_______________|__");

    // for each row there will be two columns
    for(value_type i{}; i != _length; ++i){
        std::cout << i << ' ';
        std::cout << "H |";
        for(value_type j1{}; j1 != _length; ++j1){
            std::cout << "     " << std::setw(3) << _h_bonds[i][j1].get_groupID()
                      << _h_bonds[i][j1].relativeIndex() << "|";
        }
        std::cout << std::endl;
        print_h_barrier(_length, "    |  ", "               |  "); // just for better viewing
        std::cout << "  " << "V |";
        for(value_type j2{}; j2 != _length; ++j2){
            std::cout << std::setw(3) << _v_bonds[i][j2].get_groupID()
                      << _v_bonds[i][j2].relativeIndex() << "     |";
        }
        std::cout << std::endl;

        // printing horizontal separator
        print_h_barrier(_length, "____|__", "_______________|__");
    }
    std::cout << std::endl;
}


/**
 * View bonds in the lattice by relative index. id of the site is showed
 * format : id for site or id(relative_index) for bond
 */
void SqLattice::view_bonds_by_relative_index_v2() {

    std::cout << "Bonds by id : line " << __LINE__ << endl;
    cout << "site id -1 means isolated site and 0 means connected site in bond percolation(definition)" << endl;
    print_h_barrier(15, "_", "___", "_\n");
    cout << "|(site id) (horizontal bond id(relative index))|" << endl;
    cout << "|(vertical bond id(relative index))            |" << endl;
    print_h_barrier(15, "-", "---", "-\n");
    // printing indices for columns
    std::cout << "    | ";
    for(value_type i{}; i != _length; ++i){
        std::cout << i << "                  | ";
    }
    std::cout << std::endl;

    // pringing H,V label

    print_h_barrier(_length, "    |  ", " V            H   |  ");
    print_h_barrier(_length, "____|__", "__________________|__");

    // for each row there will be two columns
    for(value_type i{}; i != _length; ++i){
        std::cout << i << ' ';
        std::cout << "H |";
        for(value_type j1{}; j1 != _length; ++j1){
            std::cout << std::setw(3) << _sites[i][j1].get_groupID() ;
            std::cout << "     " << std::setw(3) << _h_bonds[i][j1].get_groupID()
                      << _h_bonds[i][j1].relativeIndex() << "|";
        }
        std::cout << std::endl;
        print_h_barrier(_length, "    |  ", "                  |  "); // just to see a better view
        std::cout << "  " << "V |";
        for(value_type j2{}; j2 != _length; ++j2){
            std::cout << std::setw(3) << _v_bonds[i][j2].get_groupID()
                      << _v_bonds[i][j2].relativeIndex() << "        |";
        }
        std::cout << std::endl;

        // printing horizontal separator
        print_h_barrier(_length, "____|__", "__________________|__");
    }

    std::cout << std::endl;
}


/**
 * View bonds in the lattice by relative index. id of the site is showed
 * format : id(relative index) for site and only id for bond
 */
void SqLattice::view_bonds_by_relative_index_v3() {

    std::cout << "Bonds by id : line " << __LINE__ << endl;
    //cout << "site id -1 means isolated site and 0 means connected site in bond percolation(definition)" << endl;
    print_h_barrier(15, "_", "___", "_\n");
    cout << "|(site id) (horizontal bond id(relative index))|" << endl;
    cout << "|(vertical bond id(relative index))            |" << endl;
    print_h_barrier(15, "-", "---", "-\n");
    // printing indices for columns
    std::cout << "    | ";
    for(value_type i{}; i != _length; ++i){
        std::cout << i << "                | ";
    }
    std::cout << std::endl;

    // pringing H,V label

    print_h_barrier(_length, "    |  ", " V           H  |  ");
    print_h_barrier(_length, "____|__", "________________|__");

    // for each row there will be two columns
    for(value_type i{}; i != _length; ++i){
        std::cout << i << ' ';
        std::cout << "H |";
        for(value_type j1{}; j1 != _length; ++j1){
            std::cout << std::setw(3) << _sites[i][j1].get_groupID() << _sites[i][j1].relativeIndex() ;
            std::cout << "   " << std::setw(3) << _h_bonds[i][j1].get_groupID() << "|";
        }
        std::cout << std::endl;
        print_h_barrier(_length, "    |  ", "                |  "); // just to see a better view
        std::cout << "  " << "V |";
        for(value_type j2{}; j2 != _length; ++j2){
            std::cout << std::setw(3) << _v_bonds[i][j2].get_groupID() << "               |";
        }
        std::cout << std::endl;

        // printing horizontal separator
        print_h_barrier(_length, "____|__", "________________|__");
    }

    std::cout << std::endl;
}

/**
 * View bonds in the lattice by relative index. id of the site is showed
 * format : id(relative index) for site and only id for bond
 * if any site is isolated relative index is not shown
 */
void SqLattice::view_bonds_by_relative_index_v4() {

    std::cout << "Bonds by id : line " << __LINE__ << endl;
    //cout << "site id -1 means isolated site and 0 means connected site in bond percolation(definition)" << endl;
    print_h_barrier(15, "_", "___", "_\n");
    cout << "|(site id) (horizontal bond id(relative index))|" << endl;
    cout << "|(vertical bond id(relative index))            |" << endl;
    print_h_barrier(15, "-", "---", "-\n");
    // printing indices for columns
    std::cout << "    | ";
    for(value_type i{}; i != _length; ++i){
        std::cout << i << "                | ";
    }
    std::cout << std::endl;

    // pringing H,V label

    print_h_barrier(_length, "    |  ", " V           H  |  ");
    print_h_barrier(_length, "____|__", "________________|__");

    // for each row there will be two columns
    for(value_type i{}; i != _length; ++i){
        std::cout << i << ' ';
        std::cout << "H |";
        for(value_type j1{}; j1 != _length; ++j1){
            int id = _sites[i][j1].get_groupID();
            std::cout << std::setw(3) << id  ;
            if(id != -1){
                cout << _sites[i][j1].relativeIndex();
            }else{
                cout << "(-,-)    ";
            }
            std::cout << "   " << std::setw(3) << _h_bonds[i][j1].get_groupID() << "|";
        }
        std::cout << std::endl;
        print_h_barrier(_length, "    |  ", "                |  "); // just to see a better view
        std::cout << "  " << "V |";
        for(value_type j2{}; j2 != _length; ++j2){
            std::cout << std::setw(3) << _v_bonds[i][j2].get_groupID() << "               |";
        }
        std::cout << std::endl;

        // printing horizontal separator
        print_h_barrier(_length, "____|__", "________________|__");
    }

    std::cout << std::endl;
}



/**
 * View lattice (sites and bonds) by relative index.
 * format : id(relative_index)
 */
void SqLattice::view_by_relative_index() {

    std::cout << "Bonds by id : line " << __LINE__ << endl;
    // printing indices for columns
    std::cout << "    | ";
    for(value_type i{}; i != _length; ++i){
        std::cout << i << "                         | ";
    }
    std::cout << std::endl;

    // pringing H,V label

    print_h_barrier(_length, "    |  ", " V                   H   |  ");
    print_h_barrier(_length, "____|__", "_________________________|__");

    // for each row there will be two columns
    for(value_type i{}; i != _length; ++i){
        std::cout << i << ' ';
        std::cout << "H |";
        for(value_type j1{}; j1 != _length; ++j1){
            std::cout << std::setw(3) << _sites[i][j1].get_groupID() << _sites[i][j1].relativeIndex();
            std::cout << "   " << std::setw(3) << _h_bonds[i][j1].get_groupID()
                      << _h_bonds[i][j1].relativeIndex() << "|";
        }
        std::cout << std::endl;
        print_h_barrier(_length, "    |  ", "                         |  "); // just to see a better view
        std::cout << "  " << "V |";
        for(value_type j2{}; j2 != _length; ++j2){
            std::cout << std::setw(3) << _v_bonds[i][j2].get_groupID()
                      << _v_bonds[i][j2].relativeIndex() << "               |";
        }
        std::cout << std::endl;

        // printing horizontal separator
        print_h_barrier(_length, "____|__", "_________________________|__");
    }

//        view_h_bonds_extended();
//        view_v_bonds_extended();
    std::cout << std::endl;

}

/**
 * View lattice (sites and bonds) by relative index.
 * format : id(relative_index)
 */
void SqLattice::view() {

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
    for(value_type i{}; i != _length; ++i){
        std::cout << i << ' ';
        std::cout << "H |";
        for(value_type j1{}; j1 != _length; ++j1){
            std::cout << std::setw(3) << _sites[i][j1].get_groupID() ;
            std::cout << "  " << std::setw(3) << _h_bonds[i][j1].get_groupID() << "|";
        }
        std::cout << std::endl;
        print_h_barrier(_length, "    |  ", "      |  "); // just to see a better view
        std::cout << "  " << "V |";
        for(value_type j2{}; j2 != _length; ++j2){
            std::cout << std::setw(3) << _v_bonds[i][j2].get_groupID() << "     |";
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
void SqLattice::view_h_bonds()
{
    std::cout << "view horizontal bonds" << std::endl;
    std::cout << '{';
    for(value_type i{} ; i != _length ; ++i) {
        if(i!=0) std::cout << "  ";
        else std::cout << '{';
        for (value_type j{}; j != _length; ++j) {
            std::cout << _h_bonds[i][j] ;
            if(j != _length-1)
                std::cout << ',';
        }
        std::cout << '}';
        if(i != _length-1)
            std::cout << std::endl;
    }
    std::cout << '}';
    std::cout << std::endl;
}

/**
 *
 */
void SqLattice::view_v_bonds()
{
    std::cout << "view vertical bonds" << std::endl;
    std::cout << '{';
    for(value_type i{} ; i != _length ; ++i) {
        if(i!=0) std::cout << "  ";
        else std::cout << '{';
        for (value_type j{}; j != _length; ++j) {
            std::cout << _v_bonds[i][j] ;
            if(j != _length-1)
                std::cout << ',';
        }
        std::cout << '}';
        if(i != _length-1)
            std::cout << std::endl;
    }
    std::cout << '}';
    std::cout << std::endl;
}


/**
 *
 */
void SqLattice::view_h_bonds_extended(){
    std::cout << "view horizontal bonds" << std::endl;
    std::cout << '{';
    for(value_type i{} ; i != _length ; ++i) {
        if(i!=0) std::cout << "  ";
        else std::cout << '{';
        for (value_type j{}; j != _length; ++j) {
            std::cout << "(" << _h_bonds[i][j].get_groupID() << ":" << _h_bonds[i][j] << ")" ;
            if(j != _length-1)
                std::cout << ',';
        }
        std::cout << '}';
        if(i != _length-1)
            std::cout << std::endl;
    }
    std::cout << '}';
    std::cout << std::endl;
}

/**
 *
 */
void SqLattice::view_v_bonds_extended(){
    std::cout << "view vertical bonds" << std::endl;
    std::cout << '{';
    for(value_type i{} ; i != _length ; ++i) {
        if(i!=0) std::cout << "  ";
        else std::cout << '{';
        for (value_type j{}; j != _length; ++j) {
            std::cout << "(" << _v_bonds[i][j].get_groupID() << ":" << _v_bonds[i][j] << ")" ;
            if(j != _length-1)
                std::cout << ',';
        }
        std::cout << '}';
        if(i != _length-1)
            std::cout << std::endl;
    }
    std::cout << '}';
    std::cout << std::endl;
}


/**
 *
 */
void SqLattice::view_bonds_by_id(){
    std::cout << "Bonds by id : line " << __LINE__ << endl;
    cout << "Structure " << endl;
    print_h_barrier(8, "__", "___", "_\n");
    cout << "|      (horizontal bond id)|" << endl;
    cout << "|(vertical bond id)        |" << endl;
    print_h_barrier(8, "--", "---", "-\n");
    // printing indices for columns
    std::cout << "    | ";
    for(value_type i{}; i != _length; ++i){
        std::cout << i << "    | ";
    }
    std::cout << std::endl;

    // pringing H,V label

    print_h_barrier(_length, "    |  ", "V  H|  ");
    print_h_barrier(_length, "____|__", "____|__");

    // for each row there will be two columns
    for(value_type i{}; i != _length; ++i){
        std::cout << i << ' ';
        std::cout << "H |";
        for(value_type j1{}; j1 != _length; ++j1){
            std::cout << "   " << std::setw(3) << _h_bonds[i][j1].get_groupID() << "|";
        }
        std::cout << std::endl;
        std::cout << "  " << "V |";
        for(value_type j2{}; j2 != _length; ++j2){
            std::cout << std::setw(3) << _v_bonds[i][j2].get_groupID() << "   |";
        }
        std::cout << std::endl;

        // printing horizontal separator
        print_h_barrier(_length, "____|__", "____|__");
    }

//        view_h_bonds_extended();
//        view_v_bonds_extended();
    std::cout << std::endl;
}


Site& SqLattice::getSite(Index index) {
    return _sites[index.row_][index.column_];
}

//Site&& SqLattice::getSiteR(Index index) {
//    Site a = _sites[index.x_][index.y_];
//    return std::move(a);
//}

const Site& SqLattice::getSite(Index index) const {
    return _sites[index.row_][index.column_];
}

void SqLattice::setGroupID(Index index, int group_id){
    _sites[index.row_][index.column_].set_groupID(group_id);
}

void SqLattice::setGroupID(BondIndex index, int group_id){
    if(index.horizontal()){
        _h_bonds[index.row_][index.column_].set_groupID(group_id);
    }
    if(index.vertical()){
        _v_bonds[index.row_][index.column_].set_groupID(group_id);
    }
}

int SqLattice::getGroupID(Index index){
    return _sites[index.row_][index.column_].get_groupID();
}

int SqLattice::getGroupID(BondIndex index){
    if(index.horizontal()){
        return _h_bonds[index.row_][index.column_].get_groupID();
    }
    if(index.vertical()){
        return  _v_bonds[index.row_][index.column_].get_groupID();
    }
    return -1;
}

//Bond Lattice::get_h_bond(Index set_ID) {
//    return _h_bonds[id.x_][set_ID.y_];
//}
//
//Bond Lattice::get_v_bond(Index set_ID) {
//    return _v_bonds[id.x_][set_ID.y_];
//}


Bond& SqLattice::get_h_bond(Index id) {
    return _h_bonds[id.row_][id.column_];
}

Bond& SqLattice::get_v_bond(Index id) {
    return _v_bonds[id.row_][id.column_];
}


Bond& SqLattice::getBond(BondIndex index) {
    if(index.horizontal())
        return _h_bonds[index.row_][index.column_];
    if(index.vertical())
        return _v_bonds[index.row_][index.column_];
    // todo throw exception
    throw InvalidBond{"Invalid bond : line " + to_string(__LINE__)};
}

void SqLattice::reset(bool reset_all) {
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
void SqLattice::reset_bonds() {
    for(value_type i{}; i != _h_bonds.size(); ++i){
        for (int j{}; j != _h_bonds[i].size(); ++j) {
            // deactivating. automatically set group id == - and relative index == (0,0)
            // setting group id = -1 and deactivating the bond
            _h_bonds[i][j].deactivate();
            _v_bonds[i][j].deactivate();


        }
    }
}

/**
 *
 */
void SqLattice::reset_sites() {
    for(value_type i{}; i != _sites.size(); ++i){
        for(value_type j{}; j != _sites[i].size(); ++j) {
            // deactivating. automatically set group id == - and relative index == (0,0)
            // setting group id = -1 and deactivating the site
            _sites[i][j].deactivate();
        }
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
std::vector<Index> SqLattice::get_neighbor_site_indices(Index site){
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
std::vector<BondIndex> SqLattice::get_neighbor_bond_indices(BondIndex bond) {
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

std::vector<Index> SqLattice::get_neighbor_indices(BondIndex bond) {
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

/************************************************************
 * Static methods
 */
std::vector<Index> SqLattice::get_neighbor_site_indices(size_t length, Index site){
    std::vector<Index> sites(4);
    sites[0] = {(site.row_ + 1) % length, site.column_};
    sites[1] = {(site.row_ - 1 + length) % length, site.column_};
    sites[2] = {site.row_, (site.column_ + 1) % length};
    sites[3] = {site.row_, (site.column_ - 1 + length) % length};
    return sites;
}

/**
 * Periodic case only.
 * Each bond has six neibhbor bonds.
 * @param site
 * @return
 */
std::vector<BondIndex> SqLattice::get_neighbor_bond_indices(size_t length, BondIndex bond) {
    value_type next_column = (bond.column_ + 1) % length;
    value_type prev_column = (bond.column_ - 1 + length) % length;
    value_type prev_row = (bond.row_ - 1 + length) % length;
    value_type next_row = (bond.row_ + 1) % length;

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

std::vector<Index> SqLattice::get_neighbor_indices(size_t length, BondIndex bond) {
    value_type r = bond.row_;
    value_type c = bond.column_;
    vector<Index> sites(2);
    sites[0] = {r, c};
    if(bond.horizontal()){
        sites[1] = {r, (c+1) % length};
    }else{
        sites[1] = {(r+1) % length, c};
    }
    return sites;
}

std::vector<Site> SqLattice::getSites() {
    cout << "undefined : line " << __LINE__ << endl;
//    return _sites;
    return {};
}




