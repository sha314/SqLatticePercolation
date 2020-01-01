//
// Created by shahnoor on 12/23/19.
//


#include <iomanip>
#include "lattice.h"


using namespace std;

Lattice_v12::Lattice_v12(int length) {
    if(length < 2) cerr << "length should be > 2" << endl;
    _length = length;
//    cout << "Constructing Lattice object : line " << __LINE__ << endl;
    _sites_2d = std::vector<std::vector<Site_v12>>(size_t(_length));
    _bonds_2d = std::vector<std::vector<Bond_v12>>(size_t(_length*2));

//    init_lattice_v1();
    init_lattice_v2();

}

/**
 * Initialize lattice in a regular fashion. To enable natural view,
 *  for sites row and column is interchanged.
 *  for bonds a little complicated
 *
 * This function is "deprecated". Because it misplaces id for bonds.
 * `_sites` and `_bonds` supposed to be map from id to index. misplacing id is bad.
 *
 *
 */
void Lattice_v12::init_lattice_v1()  {
    int id{};

    for (int i{}; i != _length; ++i) {   // all bonds and site_index_sequence are activated by default
        _sites_2d[i].resize(size_t(_length)); // resize new row
        for (int j{}; j != _length; ++j) {
            auto a = Index(i, j);
            _sites_2d[i][j] = Site_v12(a, id++);
            _sites.emplace_back(a);
//            _sites[i][j].activate();
        }
    }

    id = 0;
    for (int i{}; i != _length; ++i) {   // all bonds and site_index_sequence are activated by default
        _bonds_2d[i].resize(size_t(_length));
        _bonds_2d[i + _length].resize(size_t(_length));

        for (int j{}; j != _length; ++j) {
            auto a = _sites_2d[i][j].get_id();

            auto k = (j+1) % _length;
            auto b1 = _sites_2d[i][k].get_id();

            _sites_2d[i][j].connectBond(id);
            _sites_2d[i][k].connectBond(id);

            _bonds_2d[i][j] = {Link(a, b1), id++}; // horizontal while assigning but vertical when viewing
            _bonds.emplace_back(Index(i, j)); // id misplaced bug
            k = (i+1) % _length;
            auto b2 = _sites_2d[k][j].get_id();
            auto m = i + _length;

            _sites_2d[i][j].connectBond(id);
            _sites_2d[k][j].connectBond(id);

            _bonds_2d[m][j] = {Link(a, b2), id++};
            _bonds.emplace_back(Index(m, j)); // id misplaced bug

        }
    }
}

/**
 * You still have to interchange row and column but bonds are much simpler to view
 * TODO
 */
void Lattice_v12::init_lattice_v2()  {
    int id{};

    _sites.resize(_length*_length);
    for (int i{}; i != _length; ++i) {   // all bonds and site_index_sequence are activated by default
        _sites_2d[i].resize(size_t(_length)); // resize new row
        for (int j{}; j != _length; ++j) {
            auto a = Index(i, j);
            _sites_2d[i][j] = Site_v12(a, id);
            _sites[id] = a;
//            _sites[i][j].activate();
            ++id;
        }
    }

    id = 0;
    _bonds.resize(2*_length*_length);
    for (int i{}; i != _length; ++i) {   // all bonds and site_index_sequence are activated by default
        _bonds_2d[i].resize(size_t(_length));
        _bonds_2d[i + _length].resize(size_t(_length));

        for (int j{}; j != _length; ++j) {
            auto a = _sites_2d[i][j].get_id();

//            auto k = (j+1) % _length; // next column
//            auto  k2 = (i+_length-1) % _length; // previous row

            auto k1 = (j+_length-1) % _length; // previous column
            auto k2 = (i+1) % _length; // next row

            auto b1 = _sites_2d[i][k1].get_id();

            _sites_2d[i][j].connectBond(id);
            _sites_2d[i][k1].connectBond(id);

            _bonds_2d[i][j] = {Link(a, b1), id}; // horizontal while assigning but vertical when viewing
            _bonds[id] = Index(i, j);
            ++id;

            auto b2 = _sites_2d[k2][j].get_id();
            auto m = i + _length;

            _sites_2d[i][j].connectBond(id);
            _sites_2d[k2][j].connectBond(id);

            _bonds_2d[m][j] = {Link(a, b2), id};
            _bonds[id] = Index(m, j);
            ++id;

        }
    }
}

void Lattice_v12::view_as_assigned() {
    cout << "Lattice_v12::view_as_assigned" << endl;
    cout << "-----------------------------------------" << endl;
    cout << "| Site                |  Link to right  |" << endl;
    cout << "| Link to bottom one  |                 |" << endl;
    cout << "-----------------------------------------" << endl;
    for(int x{}; x < _length; ++x) { // column
        cout << "--------------------";
    }
    cout << endl << "y  |" << endl;
    for(int r{0}; r < _length; ++r){ // row
        cout << setw(3) << r << "|";
        for(int c{}; c < _length; ++c){ // column
//            cout << "(x,y)=(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[r][c];

            auto bond = _bonds_2d[r][c];

            cout << setw(3) << site.get_id() << "," <<  site.get_index().get_string() << "," << setw(3) << site.get_groupID() << " |";

            cout << setw(3) << bond.get_id()
                 << "," << bond.getIndex().get_string()
                 << "," <<  setw(3) << bond.get_groupID();
            cout << " |";


        }
        cout << endl;
        cout << "   |";
        for(int c{}; c < _length; ++c){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto k = (r+_length);

            auto bond =_bonds_2d[k][c];

            cout << setw(3) << bond.get_id()
                 << "," << bond.getIndex().get_string()
                 << "," <<  setw(3) << bond.get_groupID();

            cout << " |" << setw(17) << "|" ;
        }
        cout << endl;
        cout << "---|";
        for(int x{}; x < _length; ++x) { // column
            cout << "----------------------------------";
        }
        cout << endl;
    }
    cout << "___|";
    for(int x{}; x < _length; ++x) { // column
        cout << "__________________________________";
    }
    cout << endl;
    cout << "x->|";
    for(int x{}; x < _length; ++x) { // column
        cout << "              " << setw(3) << x << "                |";
    }
    cout << endl;

}

void Lattice_v12::view_all_v2() {
    cout << "Lattice_v12::view_all" << endl;
//    print_box();
    cout << "-----------------------------------------" << endl;
    cout << "| Site                |  Link to right  |" << endl;
    cout << "| Link to bottom one  |                 |" << endl;
    cout << "-----------------------------------------" << endl;
    /*
     * Switching row and column makes horizontal bonds vertical and vice versa
     */
    for(int x{}; x < _length; ++x) { // column
        cout << "--------------------";
    }
    cout << endl << "y  |" << endl;
    for(int y{_length-1}; y >= 0; --y){ // row
        cout << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(x,y)=(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[x][y];

            auto k = (x+_length); // horizontally next site

            auto bond =_bonds_2d[k][y]; // bond to horizontally next site


            cout << setw(3) << site.get_id() << ","
                 << site.get_index().get_string() << ","
                 << setw(3) << site.get_groupID() << " |";

            cout << setw(3) << bond.get_id()
                 << "," << bond.getIndex().get_string()
                 << "," <<  setw(3) << bond.get_groupID();
            cout << " |";


        }
        cout << endl;
        cout << "   |";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
//            x = (x + _length ) % _length;
//            auto y2 = (_length - y - 1);
//            cout << "y = " << y << " y2 = " << y2 << endl;
            auto bond = _bonds_2d[x][y]; // bond to vertically next site
            cout << setw(3) << bond.get_id()
                 << "," << bond.getIndex().get_string()
                 << "," <<  setw(3) << bond.get_groupID();

            cout << " |" << setw(17) << "|" ;
        }
        cout << endl;
        cout << "---|";
        for(int x{}; x < _length; ++x) { // column
            cout << "----------------------------------";
        }
        cout << endl;
//        exit(0);
    }
    cout << "___|";
    for(int x{}; x < _length; ++x) { // column
        cout << "__________________________________";
    }
    cout << endl;
    cout << "x->|";
    for(int x{}; x < _length; ++x) { // column
        cout << "              " << setw(3) << x << "                |";
    }
    cout << endl;

}

/**
 * If this method is fllowed then we can change what to be printed to screen before printing it.
 * Just like we can edit strings.
 * This function is an exact copy of Lattice_v12::view_all_v2 except for the fact that this function
 * first stores strings inside an array then prints them.
 *
 * We can do this in a more sophisticated way:
 *      1. 2D array to store sites and bonds information as strings
 *      2. boundaries will not be stored but printed when printing elements in the end
 *      3. This way all the information is accecible as strings at any time
 *      May be in future I will do it.
 */
void Lattice_v12::view_from_str_v1() {
    vector<string> lines;
    stringstream ss;
    cout << "Lattice_v12::view_all" << endl;
    print_box();

    /*
     * Switching row and column makes horizontal bonds vertical and vice versa
     */
    for(int x{}; x < _length; ++x) { // column
        ss << "--------------------";
    }
    lines.emplace_back(ss.str());
    ss.str("");
    ss << "y  |";
    lines.emplace_back(ss.str());
    ss.str("");
    for(int y{_length-1}; y >= 0; --y){ // row
        ss << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(x,y)=(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[x][y];

            auto k = (x+_length);

            auto bond =_bonds_2d[k][y];


            ss << setw(3) << site.get_id() << ","
                 << site.get_index().get_string() << ","
                 << setw(3) << site.get_groupID() << " |";

            ss << setw(3) << bond.get_id()
                 << "," << bond.getIndex().get_string()
                 << "," <<  setw(3) << bond.get_groupID();
            ss << " |";


        }
        lines.emplace_back(ss.str());
        ss.str("");

        ss << "   |";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto y2 = (y + _length - 1) % _length;
//            cout << "y = " << y << " y2 = " << y2 << endl;
            auto bond = _bonds_2d[x][y2];
            ss << setw(3) << bond.get_id()
                 << "," << bond.getIndex().get_string()
                 << "," <<  setw(3) << bond.get_groupID();

            ss << " |" << setw(17) << "|" ;
        }
        lines.emplace_back(ss.str());
        ss.str("");

        ss << "---|";
        for(int x{}; x < _length; ++x) { // column
            ss << "----------------------------------";
        }

        lines.emplace_back(ss.str());
        ss.str("");
//        exit(0);
    }
    ss << "___|";
    for(int x{}; x < _length; ++x) { // column
        ss << "__________________________________";
    }

    lines.emplace_back(ss.str());
    ss.str("");
    ss << "x->|";
    for(int x{}; x < _length; ++x) { // column
        ss << "              " << setw(3) << x << "                |";
    }

    lines.emplace_back(ss.str());
    ss.str("");
    for(const auto &line: lines){
        cout << line << endl;
    }

}

void Lattice_v12::view_all_v1() {
    cout << "Lattice_v12::view_all_v2" << endl;
//    print_box();
    cout << "--------------------------------------" << endl;
    cout << "| Site                | Link to right |" << endl;
    cout << "| Link to bottom one  |               |" << endl;
    cout << "--------------------------------------" << endl;
    /*
     * Switching row and column makes horizontal bonds vertical and vice versa
     */
    for(int x{}; x < _length; ++x) { // column
        cout << "--------------------";
    }
    cout << endl << "y  |" << endl;
    for(int y{_length-1}; y >= 0; --y){ // row
        cout << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(x,y)=(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[x][y];

            auto k = (x+_length);

            auto bond =_bonds_2d[k][y];


            cout << setw(3) << site.get_id() << ","
                 << site.get_index().get_string() << ","
                 << setw(3) << site.get_groupID() << " |";

            cout << setw(3) << bond.get_id()
                 << "," << bond.getIndex().get_string()
                 << "," <<  setw(3) << bond.get_groupID();
            cout << " |";


        }
        cout << endl;
        cout << "   |";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
//            x = (x + _length ) % _length;
            auto y2 = (y + _length - 1) % _length;
//            cout << "y = " << y << " y2 = " << y2 << endl;
            auto bond = _bonds_2d[x][y2];
            cout << setw(3) << bond.get_id()
                 << "," << bond.getIndex().get_string()
                 << "," <<  setw(3) << bond.get_groupID();

            cout << " |" << setw(17) << "|" ;
        }
        cout << endl;
        cout << "---|";
        for(int x{}; x < _length; ++x) { // column
            cout << "----------------------------------";
        }
        cout << endl;
//        exit(0);
    }
    cout << "___|";
    for(int x{}; x < _length; ++x) { // column
        cout << "__________________________________";
    }
    cout << endl;
    cout << "x->|";
    for(int x{}; x < _length; ++x) { // column
        cout << "              " << setw(3) << x << "                |";
    }
    cout << endl;

}

/**
 * Just relative index of sites with group id and id.
 */
void Lattice_v12::view_by_relative_index_v1() {
    cout << "Lattice_v12::view_by_relative_index" << endl;
    /*
     * Switching row and column makes horizontal bonds vertical and vice versa
     */
    for(int x{}; x < _length; ++x) { // column
        cout << "-------------------";
    }
    cout << endl << "y  |" << endl;
    for(int y{_length-1}; y >= 0; --y){ // row
        cout << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(x,y)=(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[x][y];

            auto k = (x+_length);

            auto bond =_bonds_2d[k][y];


            cout << setw(3) << site.get_id() << ","
                 << get_string(site.relativeIndex()) << ","
                 << setw(3) << site.get_groupID() << " |";

//            cout << setw(3) << bond.get_id()
//                 << "," << get_string(bond.getIndex())
//                 << "," <<  setw(3) << bond.get_groupID();
//            cout << " |";


        }
//        cout << endl;
//        cout << "   |";
//        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
//        auto y2 = (y + _length - 1) % _length;
//            auto bond = _bonds_2d[x][y2];
//            cout << setw(3) << bond.get_id()
//                 << "," << get_string(bond.getIndex())
//                 << "," <<  setw(3) << bond.get_groupID();

//            cout << " |" << setw(17) << "|" ;
//        }
        cout << endl;
        cout << "---|";
        for(int x{}; x < _length; ++x) { // column
            cout << "-------------------";
        }
        cout << endl;
    }
    cout << "___|";
    for(int x{}; x < _length; ++x) { // column
        cout << "___________________";
    }
    cout << endl;
    cout << "x->|";
    for(int x{}; x < _length; ++x) { // column
        cout << "       " << setw(3) << x << "        |";
    }
    cout << endl;

}


/**
 * Just relative index of sites with group id and id.
 */
void Lattice_v12::view_by_relative_index_v2() {
    cout << "Lattice_v12::view_by_relative_index_v2" << endl;
    /*
     * Switching row and column makes horizontal bonds vertical and vice versa
     */
    for(int x{}; x < _length; ++x) { // column
        cout << "-------------------";
    }
    cout << endl << "y  |" << endl;
    for(int y{_length-1}; y >= 0; --y){ // row
        cout << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(x,y)=(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[x][y];

            cout << setw(3) << site.get_id() << ","
                 << get_string(site.relativeIndex()) << ","
                 << setw(3) << site.get_groupID() << " |";

//            auto k = (x+_length);
//            auto bond =_bonds_2d[k][y];
//            cout << setw(3) << bond.get_id()
//                 << "," << get_string(bond.getIndex())
//                 << "," <<  setw(3) << bond.get_groupID();
//            cout << " |";


        }
//        cout << endl;
//        cout << "   |";
//        for(int x{}; x < _length; ++x){ // column
////            cout << "(" << x << "," << y << ")" << endl;
//            auto bond = _bonds_2d[x][y];
//            cout << setw(3) << bond.get_id()
//                 << "," << get_string(bond.getIndex())
//                 << "," <<  setw(3) << bond.get_groupID();
//
//            cout << "   |" << setw(17) << "|" ;
//        }
        cout << endl;
        cout << "---|";
        for(int x{}; x < _length; ++x) { // column
            cout << "-------------------";
        }
        cout << endl;
    }
    cout << "___|";
    for(int x{}; x < _length; ++x) { // column
        cout << "___________________";
    }
    cout << endl;
    cout << "x->|";
    for(int x{}; x < _length; ++x) { // column
        cout << "       " << setw(3) << x << "        |";
    }
    cout << endl;

}


void Lattice_v12::view_by_id_v1() {
    cout << "Lattice_v12::view_by_id" << endl;
    print_box();
    /*
     * Switching row and column makes horizontal bonds vertical and vice versa
     */
    cout << endl << "y  |" ;
    for(int x{}; x < _length; ++x) { // column
        cout << "__________";
    }
    cout << endl;
    for(int y{_length-1}; y >= 0; --y){ // row
        cout << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(x,y)=(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[x][y];

            auto k = (x+_length);

            auto bond =_bonds_2d[k][y];


            cout << setw(3) << site.get_id()  << " |";

            cout << setw(3) << bond.get_id() << " |";


        }
        cout << endl;
        cout << "   |";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto y2 = (y + _length - 1) % _length;
            auto bond = _bonds_2d[x][y2];
            cout << setw(3) << bond.get_id() << " |" << setw(5) << "|" ;
        }
        cout << endl;
        cout << "---|";
        for(int x{}; x < _length; ++x) { // column
            cout << "----------";
        }
        cout << endl;
    }
    cout << "___|";
    for(int x{}; x < _length; ++x) { // column
        cout << "__________";
    }
    cout << endl;
    cout << "x->|";
    for(int x{}; x < _length; ++x) { // column
        cout << "  " << setw(3) << x << "    |";
    }
    cout << endl;

}

void Lattice_v12::view_by_id_v2() {
    cout << "Lattice_v12::view_by_id" << endl;
    print_box();
    /*
     * Switching row and column makes horizontal bonds vertical and vice versa
     */
    cout << endl << "y  |" ;
    for(int x{}; x < _length; ++x) { // column
        cout << "__________";
    }
    cout << endl;
    for(int y{_length-1}; y >= 0; --y){ // row
        cout << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(x,y)=(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[x][y];

            auto k = (x+_length);

            auto bond =_bonds_2d[k][y];


            cout << setw(3) << site.get_id()  << " |";

            cout << setw(3) << bond.get_id() << " |";


        }
        cout << endl;
        cout << "   |";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto bond = _bonds_2d[x][y];
            cout << setw(3) << bond.get_id() << " |" << setw(5) << "|" ;
        }
        cout << endl;
        cout << "---|";
        for(int x{}; x < _length; ++x) { // column
            cout << "----------";
        }
        cout << endl;
    }
    cout << "___|";
    for(int x{}; x < _length; ++x) { // column
        cout << "__________";
    }
    cout << endl;
    cout << "x->|";
    for(int x{}; x < _length; ++x) { // column
        cout << "  " << setw(3) << x << "    |";
    }
    cout << endl;

}

/**
 * View by group id
 */
void Lattice_v12::view_by_gid_v1() {
    cout << "Lattice_v12::view_by_gid" << endl;
    print_box();
    /*
     * Switching row and column makes horizontal bonds vertical and vice versa
     */

    cout << endl << "y  |" ;
    for(int x{}; x < _length; ++x) { // column
        cout << "__________";
    }
    cout << endl;
    for(int y{_length-1}; y >= 0; --y){ // row
        cout << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(x,y)=(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[x][y];

            auto k = (x+_length);

            auto bond =_bonds_2d[k][y];


            cout << setw(3) << site.get_groupID()  << " |";

            cout << setw(3) << bond.get_groupID() << " |";


        }
        cout << endl;
        cout << "   |";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto y2 = (y + _length - 1) % _length;
            auto bond = _bonds_2d[x][y2];
            cout << setw(3) << bond.get_groupID() << " |" << setw(5) << "|" ;
        }
        cout << endl;
        cout << "---|";
        for(int x{}; x < _length; ++x) { // column
            cout << "----------";
        }
        cout << endl;
    }
    cout << "___|";
    for(int x{}; x < _length; ++x) { // column
        cout << "__________";
    }
    cout << endl;
    cout << "x->|";
    for(int x{}; x < _length; ++x) { // column
        cout << "  " << setw(3) << x << "    |";
    }
    cout << endl;

}

/**
 * View by group id
 */
void Lattice_v12::view_by_gid_v2() {
    cout << "Lattice_v12::view_by_gid" << endl;
    print_box();
    /*
     * Switching row and column makes horizontal bonds vertical and vice versa
     */

    cout << endl << "y  |" ;
    for(int x{}; x < _length; ++x) { // column
        cout << "__________";
    }
    cout << endl;
    for(int y{_length-1}; y >= 0; --y){ // row
        cout << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(x,y)=(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[x][y];

            auto k = (x+_length);

            auto bond =_bonds_2d[k][y];


            cout << setw(3) << site.get_groupID()  << " |";

            cout << setw(3) << bond.get_groupID() << " |";


        }
        cout << endl;
        cout << "   |";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto bond = _bonds_2d[x][y];
            cout << setw(3) << bond.get_groupID() << " |" << setw(5) << "|" ;
        }
        cout << endl;
        cout << "---|";
        for(int x{}; x < _length; ++x) { // column
            cout << "----------";
        }
        cout << endl;
    }
    cout << "___|";
    for(int x{}; x < _length; ++x) { // column
        cout << "__________";
    }
    cout << endl;
    cout << "x->|";
    for(int x{}; x < _length; ++x) { // column
        cout << "  " << setw(3) << x << "    |";
    }
    cout << endl;

}



void Lattice_v12::view_sites_v1() {
    cout << "Lattice_v12::view_sites" << endl;
    for(int x{}; x < _length; ++x) { // column
        cout << "--------------------";
    }
    cout << endl << "y  |" << endl;
    // generating row backwards
    for(int y{_length-1}; y >= 0; --y){ // row
        cout << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[x][y]; // first column then row index

            cout << setw(3) << site.get_id() << "," <<  setw(4) << site.get_index() << "," << setw(3) << site.get_groupID() << "  |";

        }
        cout << endl;
        cout << "---|";
        for(int x{}; x < _length; ++x) { // column
            cout << "-------------------";
        }
        cout << endl;
    }

    cout << "___|";
    for(int x{}; x < _length; ++x) { // column
        cout << "___________________";
    }
    cout << endl;
    cout << "x->|";
    for(int x{}; x < _length; ++x) { // column
        cout << "     " << setw(3) << x << "          |";
    }
    cout << endl;

}

void Lattice_v12::view_sites_by_relative_index_v1() {
    cout << "Lattice_v12::view_sites_by_relative_index" << endl;

    cout << endl << "y  |" ;
    for(int x{}; x < _length; ++x) { // column
        cout << "______________";
    }
    cout << endl;
    // generating row backwards
    for(int y{_length-1}; y >= 0; --y){ // row
        cout << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[x][y]; // first column then row index

            cout << setw(3) << site.get_groupID() << "," <<  get_string(site.relativeIndex())  << "  |";

        }
        cout << endl;
        cout << "---|";
        for(int x{}; x < _length; ++x) { // column
            cout << "--------------";
        }
        cout << endl;
    }

    cout << "___|";
    for(int x{}; x < _length; ++x) { // column
        cout << "______________";
    }
    cout << endl;
    cout << "x->|";
    for(int x{}; x < _length; ++x) { // column
        cout << "  " << setw(3) << x << "        |";
    }
    cout << endl;

}

void Lattice_v12::view_sites_by_id_gid_v1() {
    cout << "Lattice_v12::view_sites_by_id" << endl;
    for(int x{}; x < _length; ++x) { // column
        cout << "-----------";
    }


    cout << endl << "  y|" << endl;
    for(int y{_length-1}; y >= 0; --y){ // row
        cout << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[x][y];

            cout << setw(3) << site.get_id() << ","  << setw(3) << site.get_groupID() << "  |";

        }
        cout << endl;
        cout << "---|";
        for(int x{}; x < _length; ++x) { // column
            cout << "----------";
        }
        cout << endl;
    }

    cout << "x->|";
    for(int x{}; x < _length; ++x) { // column
        cout << "  " << setw(3) << x << "    |";
    }
    cout << endl;

}

void Lattice_v12::view_sites_by_id_v1() {
    cout << "Lattice_v12::view_sites_by_id" << endl;
    for(int x{}; x < _length; ++x) { // column
        cout << "-------";
    }


    cout << endl << "  y|" << endl;
    for(int y{_length-1}; y >= 0; --y){ // row
        cout << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[x][y];

            cout << setw(3) << site.get_id()  << "  |";

        }
        cout << endl;
        cout << "---|";
        for(int x{}; x < _length; ++x) { // column
            cout << "------";
        }
        cout << endl;
    }

    cout << "___|";
    for(int x{}; x < _length; ++x) { // column
        cout << "______";
    }
    cout << endl;

    cout << "x->|";
    for(int x{}; x < _length; ++x) { // column
        cout << setw(3) << x << "  |";
    }
    cout << endl;

}



void Lattice_v12::view_sites_by_gid_v1() {
    cout << "Lattice_v12::view_sites_by_id" << endl;
    cout << endl << "y  |" ;
    for(int x{}; x < _length; ++x) { // column
        cout << "_____";
    }
    cout << endl;
    for(int y{_length-1}; y >= 0; --y){ // row
        cout << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto site =  _sites_2d[x][y];

            cout << setw(3) << site.get_groupID() << " |";

        }
        cout << endl;
        cout << "---|";
        for(int x{}; x < _length; ++x) { // column
            cout << "-----";
        }
        cout << endl;
    }
    cout << "___|";
    for(int x{}; x < _length; ++x) { // column
        cout << "_____";
    }
    cout << endl;
    cout << "x->|";
    for(int x{}; x < _length; ++x) { // column
        cout << setw(3) << x << " |";
    }
    cout << endl;

}


/**
 * print two bonds for one site. one horizontal and one vertical
 * @param i
 * @param j
 */
void Lattice_v12::print_box() {
//    cout << "-----------------------------------------------------------" << endl;
//    cout << "Site id, Site index, site gid  |  Bond id, Link, Bond gid  |" << endl;
//    cout << "Bond id, Link, Bond gid        |                           |" << endl;
//    cout << "-----------------------------------------------------------" << endl;

//    cout << "--------------------------------------" << endl;
//    cout << "| Site               |  Link to right  |" << endl;
//    cout << "| Link to above one  |                 |" << endl;
//    cout << "--------------------------------------" << endl;

    cout << "--------------------------------------" << endl;
    cout << "| Site               |  Link to right  |" << endl;
    cout << "| Link to bottom one |                 |" << endl;
    cout << "--------------------------------------" << endl;
}

std::string Lattice_v12::get_string(const Link &lnk) const {
    stringstream ss;
    ss << "(" << setw(2) << lnk.nodeA() << "&" << setw(2) << lnk.nodeB() << ")";
    return ss.str();
}

std::string Lattice_v12::get_string(const IndexRelative& index) const{
    stringstream ss;
    ss << '(' << std::setw(3) << index.x_ << ',' << std::setw(3) << index.y_ << ')';
    return ss.str();
}

std::string Lattice_v12::get_string(const Index& index) const{
    stringstream ss;
    ss << '(' << std::setw(2) << index.row_ << ',' << std::setw(2) << index.column_ << ')';
    return ss.str();
}

Bond_v12& Lattice_v12::getBond(int id) {
    auto i = _bonds[id];
    return  _bonds_2d[i.row_][i.column_];
}

Site_v12& Lattice_v12::getSite(int id) {
    auto i = _sites[id];
    return  _sites_2d[i.row_][i.column_];
}

int Lattice_v12::getGroupIDBond(int id) {
    auto i = _bonds[id];
    return  _bonds_2d[i.row_][i.column_].get_groupID();
}

int Lattice_v12::getGroupIDSite(int id) {
    auto i = _sites[id];
    return  _sites_2d[i.row_][i.column_].get_groupID();
}

int Lattice_v12::getGroupIDBond(Index i) {
    return  _bonds_2d[i.row_][i.column_].get_groupID();
}

int Lattice_v12::getGroupIDSite(Index i) {
    return  _sites_2d[i.row_][i.column_].get_groupID();
}

void Lattice_v12::setGroupIDBond(int id, int group_id) {
    auto i = _bonds[id];
    _bonds_2d[i.row_][i.column_].set_groupID(group_id);
}

void Lattice_v12::setGroupIDSite(int id, int group_id) {
    auto i = _sites[id];
    _sites_2d[i.row_][i.column_].set_groupID(group_id);
}

void Lattice_v12::setGroupIDBond(Index index, int group_id) {
    _bonds_2d[index.row_][index.column_].set_groupID(group_id);
}

void Lattice_v12::setGroupIDSite(Index index, int group_id) {
    _sites_2d[index.row_][index.column_].set_groupID(group_id);
}

//void Lattice_v12::get_neighbors(Index site, std::vector<Index> &site_neighbor, std::vector<Index> &bond_neighbor) {
//    site_neighbor.clear();
//    bond_neighbor.clear();
//
//    value_type prev_column  = (site.column_ + _length - 1) % _length;
//    value_type prev_row     = (site.row_    + _length - 1) % _length;
//    value_type next_row     = (site.row_    + 1)           % _length;
//    value_type next_column  = (site.column_ + 1)           % _length;
//
//    // 1 level inside the lattice
//    // not in any the boundary
//    site_neighbor.resize(4);
//    site_neighbor[0] = {site.row_, next_column};
//    site_neighbor[1] = {site.row_, prev_column};
//    site_neighbor[2] = {next_row, site.column_};
//    site_neighbor[3] = {prev_row, site.column_};
//
//    bond_neighbor.reserve(4);
//
////    if(!_lattice.getSite(site_neighbor[0]).isActive()) {
////        bond_neighbor.push_back({BondType::Horizontal, site.row_, site.column_});
////    }
////    if(!_lattice.getSite(site_neighbor[1]).isActive()){
////        bond_neighbor.push_back({BondType::Horizontal, site.row_, prev_column});
////    }
////    if(!_lattice.getSite(site_neighbor[2]).isActive()){
////        bond_neighbor.push_back({BondType::Vertical,    site.row_, site.column_});
////    }
////    if(!_lattice.getSite(site_neighbor[3]).isActive()) {
////        bond_neighbor.push_back({BondType::Vertical, prev_row, site.column_});
////    }
//    bond_neighbor.resize(4);
//
//    bond_neighbor[0] = {site.row_, site.column_};
//    bond_neighbor[1] = {site.row_, prev_column};
//    bond_neighbor[2] = {site.row_ + _length, site.column_};
//    bond_neighbor[3] = {prev_row + _length, site.column_};
//
//    for(size_t i{}; i < bond_neighbor.size();++i){
//        auto a = getBondLink(bond_neighbor[i]);
//        cout << a << endl;
//    }
//}


void Lattice_v12::view_sites_list() {
    cout << "Lattice_v12::view_sites_list" << endl;
    for(size_t i{}; i < _sites.size(); ++i){
        auto tmp = _sites_2d[_sites[i].row_][_sites[i].column_];
        auto bs = tmp.connectedBondIDs();
        cout << "[" << tmp.get_id() << ", " << _sites[i] << ", " << tmp.get_groupID() << " ]{";
        for(auto b: bs){
            cout << b << ",";
        }
        cout << "}" << endl;
    }
}

void Lattice_v12::view_bonds_list() {
    cout << "Lattice_v12::view_bonds_list" << endl;
    for(size_t i{}; i < _bonds.size(); ++i){
        auto tmp = _bonds_2d[_bonds[i].row_][_bonds[i].column_];
        auto ss = tmp.connectedSites();
        cout << "[" << tmp.get_id() << ", " << _bonds[i] << ", " << tmp.get_groupID() << " ]{";
        for(auto s: ss){
            cout << s << ",";
        }
        cout << "}" << endl;
    }
}

std::vector<Index> Lattice_v12::get_neighbor_sites_of_site(Index site) {
    vector<Index> site_neighbor(4);


    value_type prev_column  = (site.column_ + _length - 1) % _length;
    value_type prev_row     = (site.row_    + _length - 1) % _length;
    value_type next_row     = (site.row_    + 1)           % _length;
    value_type next_column  = (site.column_ + 1)           % _length;

    // 1 level inside the lattice
    // not in any the boundary

    site_neighbor[0] = {site.row_, next_column};
    site_neighbor[1] = {site.row_, prev_column};
    site_neighbor[2] = {next_row, site.column_};
    site_neighbor[3] = {prev_row, site.column_};

    return site_neighbor;
}

std::vector<Index> Lattice_v12::get_neighbor_sites_of_site(int id) {
    auto site = _sites[id];
    return get_neighbor_sites_of_site(site);
}

IndexRelative& Lattice_v12::getRelativeIndex(Index index) {
    return _sites_2d[index.row_][index.column_].relativeIndex();
}

void Lattice_v12::setRelativeIndex(Index index, IndexRelative ir) {
    _sites_2d[index.row_][index.column_].relativeIndex(ir);
}


IndexRelative& Lattice_v12::getRelativeIndex(int id) {
    Index index = _sites[id];
    return getRelativeIndex(index);
}

void Lattice_v12::setRelativeIndex(int id, IndexRelative ir) {
    Index index = _sites[id];
    setRelativeIndex(index, ir);
}

void Lattice_v12::reset(bool reset_all) {
    reset_sites();
    reset_bonds();
}

void Lattice_v12::reset_sites() {
    for(size_t i{}; i < _sites_2d.size(); ++i){
        for(size_t j{}; j < _sites_2d[i].size(); ++j){
            _sites_2d[i][j].reset();
        }
    }
}

void Lattice_v12::reset_bonds() {
    for(size_t i{}; i < _bonds_2d.size(); ++i){
        for(size_t j{}; j < _bonds_2d[i].size(); ++j){
            _bonds_2d[i][j].reset();
        }
    }
}
