//
// Created by shahnoor on 10/2/2017.
//

#include <iomanip>
#include "lattice.h"

using namespace std;

/**
 *
 * @param length   -> length of the lattice
 */
SqLattice::SqLattice(value_type length)
        : _length{length}
{
    cout << "Constructing Lattice object : line " << __LINE__ << endl;
    cout << "Bond resetting is disabled : line " << __LINE__ << endl;
    _sites = std::vector<std::vector<Site>>(_length);
    _h_bonds = std::vector<std::vector<Bond>>(_length);
    _v_bonds = std::vector<std::vector<Bond>>(_length);
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
SqLattice::SqLattice(value_type length, bool activate_bonds, bool activate_sites)
        : _length{length}, _bond_resetting_flag{!activate_bonds}, _site_resetting_flag{!activate_sites}
{
    cout << "Constructing Lattice object : line " << __LINE__ << endl;
    cout << "Bond resetting is disabled : line " << __LINE__ << endl;
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

/**
 *
 */
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
                     << std::setw(2) << _sites[i][j]._id.col_ << ')';
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
    cout << "   ";
    for(value_type j{}; j != _length; ++ j){
        cout << " " << setw(4) << j << setw(8) << ' ';
    }
    cout << endl << " __|";
    for(value_type j{}; j != _length; ++ j){
        cout << "_ _ _ _ _ _ _";
    }
    cout << endl;
    for(value_type i{} ; i != _length; ++i){
        cout << setw(3) << i << "|";
        for(value_type j{} ; j != _length ; ++ j){
            if(_sites[i][j].get_groupID() == -1){
                // left blank
                cout << setw(4) << _sites[i][j].get_groupID() << setw(9) << ' ';
                continue;
            }
            cout << setw(4) << _sites[i][j].get_groupID() << _sites[i][j].relativeIndex();
        }
        cout << endl;
    }
    cout << " __|";
    for(value_type j{}; j != _length; ++ j){
        cout << "_ _ _ _ _ _ _";
    }
    cout << endl;
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
            std::cout << "(" << _h_bonds[i][j].groupID() << ":" << _h_bonds[i][j] << ")" ;
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
            std::cout << "(" << _v_bonds[i][j].groupID() << ":" << _v_bonds[i][j] << ")" ;
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


void SqLattice::activate_site(Index index) {
    if(debug_activate_site){
        cout << "activating site " << index << endl;
    }
    _sites[index.row_][index.col_].activate();
}

void SqLattice::activate_h_bond(Bond bond){
    if(bond.isHorizontal()) {// horizontal bond
        _h_bonds[bond.id().row_][bond.id().col_].activate();
    }
    else{
        bond.activate();
        cout << bond << " is not a valid horizontal bond : line " << __LINE__ << endl;
    }
}

void SqLattice::activate_v_bond(Bond bond){
    if(bond.isVertical()){// vertical bond
        _v_bonds[bond.id().row_][bond.id().col_].activate();
    }
    else{
        bond.activate();
        cout << bond << " is not a valid vertical bond : line " << __LINE__ << endl;
    }
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
    _sites[index.row_][index.col_].deactivate();
}

void SqLattice::deactivate_h_bond(Bond bond){
    if(bond.isHorizontal()) {// horizontal bond
        _h_bonds[bond.id().row_][bond.id().col_].deactivate();
    }
    else{
        bond.activate();
        cout << bond << " is not a valid horizontal bond : line " << __LINE__ << endl;
    }
}

void SqLattice::deactivate_v_bond(Bond bond){
    if(bond.isVertical()){// vertical bond
        _v_bonds[bond.id().row_][bond.id().col_].deactivate();
    }
    else{
        bond.activate();
        cout << bond << " is not a valid vertical bond : line " << __LINE__ << endl;
    }
}


void SqLattice::deactivate_bond(Bond bond) {
    // check if the bond is vertical or horizontal
    // then call appropriate function to activate _h_bond or _v_bond
    if(bond.isHorizontal()){
        if(_h_bonds[bond.id().row_][bond.id().col_].isActive()){
            cout << "Bond is already activated : line " << __LINE__ << endl;
        }
        _h_bonds[bond.id().row_][bond.id().col_].deactivate();
    }
    else if(bond.isVertical())
    {
        if(_v_bonds[bond.id().row_][bond.id().col_].isActive()){
            cout << "Bond is already activated : line " << __LINE__ << endl;
        }
        _v_bonds[bond.id().row_][bond.id().col_].deactivate();
    }
    else{
        bond.activate();
        cout << bond << " is not a valid bond : line " << __LINE__ << endl;
    }
}

Site& SqLattice::getSite(Index index) {
    return _sites[index.row_][index.col_];
}

//Site&& SqLattice::getSiteR(Index index) {
//    Site a = _sites[index.x_][index.y_];
//    return std::move(a);
//}

const Site& SqLattice::getSite(Index index) const {
    return _sites[index.row_][index.col_];
}

void SqLattice::GroupID(Index index, int group_id){
    _sites[index.row_][index.col_].set_groupID(group_id);
}

int SqLattice::GroupID(Index index){
    return _sites[index.row_][index.col_].get_groupID();
}

//Bond Lattice::get_h_bond(Index set_ID) {
//    return _h_bonds[id.x_][set_ID.y_];
//}
//
//Bond Lattice::get_v_bond(Index set_ID) {
//    return _v_bonds[id.x_][set_ID.y_];
//}


Bond& SqLattice::get_h_bond(Index id) {
    return _h_bonds[id.row_][id.col_];
}

Bond& SqLattice::get_v_bond(Index id) {
    return _v_bonds[id.row_][id.col_];
}


Bond& SqLattice::getBond(BondIndex index) {
    if(index.horizontal())
        return _h_bonds[index.row_][index.column_];
    if(index.vertical())
        return _v_bonds[index.row_][index.column_];
    // todo throw exception
    throw InvalidBond{"Invalid bond : line " + to_string(__LINE__)};
}

void SqLattice::reset() {
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
//            _h_bonds[i][j].set_groupID(-1);
//            _v_bonds[i][j].set_groupID(-1);

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
//            _sites[i][j].set_groupID(-1);
            // setting group id = -1 and deactivating the site
            _sites[i][j].deactivate();
        }
    }
}

/**
 * get Site ID's, i.e., Cluster ID's (separated by space character) as string
 * which can be used to write ids of the lattice to a file to analyze later to obtain
 *      One time data generation
 * @return
 */
string SqLattice::getLatticeIDs() {
    stringstream ss;
    for(value_type i{} ; i != _length ; ++i) {
        for (value_type j{}; j != _length; ++j) {
            ss << _sites[i][j].get_groupID() << " ";
        }
        if(i != _length-1)
            ss << std::endl;
    }
    return ss.str();
}


/******************************************************************************
 * Get Neighbor from given index
 ******************************************************************************/
std::vector<Index> SqLattice::getNeighbrs(Index site){
    std::vector<Index> sites(4);
    sites[0] = {(site.row_ + 1) % _length, site.col_};
    sites[1] = {(site.row_ - 1 + _length) % _length, site.col_};
    sites[2] = {site.row_, (site.col_ + 1) % _length};
    sites[3] = {site.row_, (site.col_ - 1 + _length) % _length};
    return sites;
}


/**
 *
 */
void SqLattice::view_bonds_by_id(){
    std::cout << "Bonds by id : line " << __LINE__ << endl;
    // printing indices for columns
    std::cout << "    | ";
    for(value_type i{}; i != _length; ++i){
        std::cout << i << "    | ";
    }
    std::cout << std::endl;

    // pringing H,V label
    std::cout << "    |";
    for(value_type i{}; i != _length; ++i){
        std::cout << " V  H |";
    }
    std::cout << std::endl;

    print_horizontal_separator();

    // for each row there will be two columns
    for(value_type i{}; i != _length; ++i){
        std::cout << i << ' ';
        std::cout << "H |";
        for(value_type j1{}; j1 != _length; ++j1){
            std::cout << "   " << std::setw(3) << _h_bonds[i][j1].groupID() << "|";
        }
        std::cout << std::endl;
        std::cout << "  " << "V |";
        for(value_type j2{}; j2 != _length; ++j2){
            std::cout << std::setw(3) << _v_bonds[i][j2].groupID() << "   |";
        }
        std::cout << std::endl;

        // printing horizontal separator
        print_horizontal_separator();
    }

//        view_h_bonds_extended();
//        view_v_bonds_extended();
    std::cout << std::endl;
}

/**
 *
 */
void SqLattice::print_horizontal_separator() const {
    cout << "____|__";
    for(value_type k{}; k != _length; ++k){
            cout << "____|__";
        }
    cout << endl;
}

