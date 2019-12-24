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

    int id{};

    for (int i{}; i != _length; ++i) {   // all bonds and site_index_sequence are activated by default
        _sites_2d[i].resize(size_t(length)); // resize new row
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
            _bonds_2d[i][j] = {Link(a, b1), id++}; // horizontal while assigning but vertical when viewing
            _bonds.emplace_back(Index(i, j));
            k = (i+1) % _length;
            auto b2 = _sites_2d[k][j].get_id();
            auto m = i + length;
            _bonds_2d[m][j] = {Link(a, b2), id++};
            _bonds.emplace_back(Index(m, j));

        }
    }

}

void Lattice_v12::view_as_assigned() {
    cout << "Lattice_v12::view_as_assigned" << endl;
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

void Lattice_v12::view_all() {
    cout << "Lattice_v12::view_all" << endl;
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
            auto bond = _bonds_2d[x][y];
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


void Lattice_v12::view_by_relative_index() {
    cout << "Lattice_v12::view_all" << endl;
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
                 << get_string(site.relativeIndex()) << ","
                 << setw(3) << site.get_groupID() << " |";

            cout << setw(3) << bond.get_id()
                 << "," << get_string(bond.getIndex())
                 << "," <<  setw(3) << bond.get_groupID();
            cout << " |";


        }
        cout << endl;
        cout << "   |";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto bond = _bonds_2d[x][y];
            cout << setw(3) << bond.get_id()
                 << "," << get_string(bond.getIndex())
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


void Lattice_v12::view_by_id() {
    cout << "Lattice_v12::view_by_id" << endl;
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
void Lattice_v12::view_by_gid() {
    cout << "Lattice_v12::view_by_id" << endl;
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



void Lattice_v12::view_sites() {
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

void Lattice_v12::view_sites_by_relative_index() {
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

void Lattice_v12::view_sites_by_id_gid() {
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

void Lattice_v12::view_sites_by_id() {
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



void Lattice_v12::view_sites_by_gid() {
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
void Lattice_v12::print_box(int i, int j) {

}

std::string Lattice_v12::get_string(const Link &lnk) const {
    stringstream ss;
    ss << "(" << setw(2) << lnk.nodeA() << "&" << setw(2) << lnk.nodeB() << ")";
    return ss.str();
}

std::string Lattice_v12::get_string(const IndexRelative& index) const{
    stringstream ss;
    ss << '(' << std::setw(2) << index.x_ << ',' << std::setw(2) << index.y_ << ')';
    return ss.str();
}

std::string Lattice_v12::get_string(const Index& index) const{
    stringstream ss;
    ss << '(' << std::setw(2) << index.row_ << ',' << std::setw(2) << index.column_ << ')';
    return ss.str();
}

Bond_v12 Lattice_v12::getBond(int id) {
    auto i = _bonds[id];
    return  _bonds_2d[i.row_][i.column_];
}

Site_v12 Lattice_v12::getSite(int id) {
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
