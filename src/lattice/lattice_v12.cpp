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
        _sites[i].resize(size_t(length));
        for (int j{}; j != _length; ++j) {
            _sites[i][j] = Site_v12(Index(i, j), ++id);
//            _sites[i][j].activate();
        }
    }

    for (int i{}; i != _length; ++i) {   // all bonds and site_index_sequence are activated by default
        _bonds[i] = std::vector<Bond_v12>(size_t(_length));
        _bonds[i + _length] = std::vector<Bond_v12>(size_t(_length));

        for (int j{}; j != _length; ++j) {
            auto a = _sites[i][j].get_id();
            auto k = (j+1) % _length;
            auto b1 = _sites[i][k].get_id();
            k = (i+1) % _length;
            auto b2 = _sites[k][j].get_id();
            _bonds[i][j] = {Link(a, b1), ++id};
            _bonds[i + _length][j] = {Link(a, b2), ++id};

        }
    }

}


void Lattice_v12::view_all() {

    for(int y{_length-1}; y >= 0; --y){ // row
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto site =  _sites[y][x];
            auto bond = _bonds[y][x];

            cout << setw(3) << site.get_id() << "," <<  setw(7) << site.get_index() << "," << setw(3) << site.get_groupID() << "  |";

            cout << setw(3) << bond.get_id()
                 << "," << setw(7) << bond.getIndex()
                 << "," <<  setw(3) << bond.get_groupID();
            cout << " |";


        }
        cout << endl;
        for(int x{}; x < _length; ++x){ // column
            auto k = (y+_length) % _length;

            auto bond =_bonds[k][x];
            cout << setw(3) << bond.get_id()
                 << "," << setw(7) << bond.getIndex()
                 << "," <<  setw(3) << bond.get_groupID();

            cout << "|" << setw(23) << "|" ;
        }
        cout << endl;
        cout << "--------------------------------------------------";
        cout << endl;
    }

}


void Lattice_v12::view_sites() {

    for(int y{_length-1}; y >= 0; --y){ // row
        for(int x{}; x < _length; ++x){ // column
//            cout << "(" << x << "," << y << ")" << endl;
            auto site =  _sites[y][x];
            auto bond = _bonds[y][x];

            cout << setw(3) << site.get_id() << "," <<  setw(7) << site.get_index() << "," << setw(3) << site.get_groupID() << "  |";

        }
        cout << "--------------------------------------------------";
        cout << endl;
    }

    for(int x{}; x < _length; ++x) { // column
        cout << setw(5) << x << " | ";
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
