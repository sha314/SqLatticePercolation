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
    _sites = std::vector<std::vector<Site_v12>>(size_t(_length));
    _bonds = std::vector<std::vector<Bond_v12>>(size_t(_length*2));

    int id{};

    for (int i{}; i != _length; ++i) {   // all bonds and site_index_sequence are activated by default
        _sites[i].resize(size_t(length)); // resize new row
        for (int j{}; j != _length; ++j) {
            _sites[i][j] = Site_v12(Index(i, j), ++id);
//            _sites[i][j].activate();
        }
    }

//    id = 0;
    for (int i{}; i != _length; ++i) {   // all bonds and site_index_sequence are activated by default
        _bonds[i].resize(size_t(_length));
        _bonds[i + _length].resize(size_t(_length));

        for (int j{}; j != _length; ++j) {
            auto a = _sites[i][j].get_id();

            auto k = (j+1) % _length;
            auto b1 = _sites[i][k].get_id();
            _bonds[i][j] = {Link(a, b1), ++id}; // horizontal while assigning but vertical when viewing

            k = (i+1) % _length;
            auto b2 = _sites[k][j].get_id();
            _bonds[i + _length][j] = {Link(a, b2), ++id};

        }
    }

}

void Lattice_v12::view_as_assigned() {
    cout << "Lattice_v12::view_as_assigned" << endl;

    for(int r{0}; r < _length; ++r){ // row
        for(int c{}; c < _length; ++c){ // column
//            cout << "(x,y)=(" << x << "," << y << ")" << endl;
            auto site =  _sites[r][c];

            auto bond = _bonds[r][c];

            cout << setw(3) << site.get_id() << "," <<  site.get_index().get_string() << "," << setw(3) << site.get_groupID() << "  |";

            cout << setw(3) << bond.get_id()
                 << "," << bond.getIndex().get_string()
                 << "," <<  setw(3) << bond.get_groupID();
            cout << " |";


        }
        cout << endl;
        for(int c{}; c < _length; ++c){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto k = (r+_length);

            auto bond =_bonds[k][c];

            cout << setw(3) << bond.get_id()
                 << "," << bond.getIndex().get_string()
                 << "," <<  setw(3) << bond.get_groupID();

            cout << "|" << setw(17) << "|" ;
        }
        cout << endl;
        cout << "--------------------------------------------------";
        cout << endl;
    }

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
            auto site =  _sites[x][y];

            auto k = (x+_length);

            auto bond =_bonds[k][y];


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
            auto bond = _bonds[x][y];
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
            auto site =  _sites[x][y]; // first column then row index

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
            auto site =  _sites[x][y];

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
            auto site =  _sites[x][y];

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
    for(int x{}; x < _length; ++x) { // column
        cout << "------";
    }


    cout << endl << "  y|" << endl;
    for(int y{_length-1}; y >= 0; --y){ // row
        cout << setw(3) << y << "|";
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto site =  _sites[x][y];

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
        cout << "______";
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
