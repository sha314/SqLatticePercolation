//
// Created by shahnoor on 2/15/21.
//

#include <set>
#include "lattice_v14.h"
#include "bond_v13.h"


using namespace std;

/**
 * horizontal axis->y axis. y increases as you move right. columns
 * vertical axis->x axis. x increases as you move down. rows
 * L = length of the lattice
 *
 * we can get id from (row,col)=(x,y)
 * id of site = col + L * row
 *
 * or we can get row and col from id
 * row = id%L
 * col = id - L*row
 *
 *
 *
 * There can be two types of bonds. Horizontal bonds (0) and vertical bonds(1). If we associate only right and bottom bond
 * to a site then left and top bond can be found if we look for right and bottom bond of left and top sites.
 * Or,
 * right bond of left site = left bond of current site.
 * bottom bond of top site = top bond of current site.
 *
 * Horizontal bonds will have the same id as the site associated with it. So, there will be L^2 number of horizontal bonds.
 * To ensure unique ids, we need to add L^2 to all the ids of vertical bonds.
 * site id        = s
 * right bond id  = s
 * bottom bond id = s + L^2
 *
 */
void Lattice_v14::init_lattice() {
    int length_squared = length*length;
    for (int rr=0; rr < length; ++rr){
        for(int cc=0; cc < length; ++ cc){
            int s_index = rr*length + cc;
            int v_bond_index = s_index+length_squared;
            site_matrix[s_index] = Site_v14(rr, cc);
            site_matrix[s_index].set_id(s_index);

            bond_matrix[s_index] = Bond_v14(rr, cc, 0);
            bond_matrix[s_index].set_id(s_index);
            bond_matrix[v_bond_index] = Bond_v14(rr, cc, 1);
            bond_matrix[v_bond_index].set_id(v_bond_index)   ;

            site_ids.push_back(s_index);
            bond_ids.push_back(s_index);
            bond_ids.push_back(v_bond_index);
        }
    }


//# print("self.bond_matrix ", self.bond_matrix)

}

///**
// * Which site is neighbor to which bond and vice versa.
// *
// * Each bond connects two sites. And each site is connected to four bonds.
// *
// * Each site will contain id of four bonds it's connected to.
// * Each bond will contain id of the two sites it is connecting.
// */
//void Lattice_v14::init_neighbor_ids() {
//
//    for (int rr=0; rr < length; ++rr) {
//        for (int cc = 0; cc < length; ++cc) {
//            int s0_index = rr*length + cc;
//            std::vector<int> bonds = find_neighbor_bonds(s0_index);
//            site_matrix[s0_index].add_connecting_bonds(bonds);
//            for(auto bb: bonds){
//                bond_matrix[bb].add_connected_site(s0_index);
//            }
//        }
//    }
//
//}

void Lattice_v14::reset() {
    for(Site_v14& site: site_matrix){
        site.reset();
    }
    for(Bond_v14& bond: bond_matrix){
        bond.reset();
    }
}

int Lattice_v14::bottom_bond_of_site(int s0_index) {
    return s0_index + site_count;
}

/**
 * top index of a site 'm' is the bottom index of some other site 'n'.
            n = (m + 2*L*L - L) (mod L*L)

 * @param s0_index
 * @return
 */
int Lattice_v14::top_bond_of_site(int s0_index) {

    int of_bottom_site = (s0_index+bond_count-length)%site_count;
    return bottom_bond_of_site(of_bottom_site);
}

int Lattice_v14::right_bond_of_site(int s0_index) {
    return s0_index;
}

/***
left index of a site 'm' is the right index of some other site 'n'.
n =
*/
int Lattice_v14::left_bond_of_site(int s0_index) {

    int row = get_row_from_id(s0_index);
    int col = get_col_from_id(s0_index);

    int left_site = row*length + (col + length - 1) % length;
//# print("left of ", s0_index, " is the right of ", left_site)
    return right_bond_of_site(left_site);
}

int Lattice_v14::get_row_from_id(int s0_index) {
    int row = s0_index / length;
    int col = s0_index % length;
    return row;
}

int Lattice_v14::get_col_from_id(int s0_index) {
    int row = s0_index / length;
    int col = s0_index % length;
    return col;
}

/**
 * There can be two types of bonds. Horizontal bonds (0) and vertical bonds(1). If we associate only right and bottom bond
 * to a site then left and top bond can be found if we look for right and bottom bond of left and top sites.
 * Or,
 * right bond of left site = left bond of current site.
 * bottom bond of top site = top bond of current site.
 *
 * Horizontal bonds will have the same id as the site associated with it. So, there will be L^2 number of horizontal bonds.
 * To ensure unique ids, we need to add L^2 to all the ids of vertical bonds.
 * site id        = s
 * right bond id  = s
 * bottom bond id = s + L^2
 *
 * @param s0_index : current site id
 * @return
 */
std::vector<int> Lattice_v14::find_neighbor_bonds(int s0_index) {
    int right_bond = right_bond_of_site(s0_index);
//# bottom_bond = s0_index + self.site_count
    int bottom_bond = bottom_bond_of_site(s0_index);
//# left_bond = (s0_index + self.length - 1) % self.length
    int left_bond = left_bond_of_site(s0_index);
//# top_bond = (s0_index+self.bond_count-self.length)%self.site_count + self.site_count
    int top_bond = top_bond_of_site(s0_index);
    return {right_bond, bottom_bond, left_bond, top_bond};
}

std::vector<int> Lattice_v14::get_neighbor_bonds(int s0_index) {
    return site_matrix[s0_index].get_connecting_bonds(length);
}

std::vector<int> Lattice_v14::get_neighbor_sites(int b0_index) {
    return bond_matrix[b0_index].get_connected_sites();
}

std::vector<int> Lattice_v14::get_site_neighbor_of_site(int s0_index) {
//# print("get_site_neighbor_of_site : ", s0_index)
    std::vector<int> bonds = site_matrix[s0_index].get_connecting_bonds(length);
//# print("bonds ", bonds)
    std::vector<int>  out_list, nn;
    for (auto bb : bonds){
        nn = get_neighbor_sites(bb);
        for(int n=0;n < nn.size(); ++n){
            if (nn[n] == s0_index){
                nn.erase(nn.begin() + n);
                break;
            }
        }
        out_list.insert(out_list.end(), nn.begin(), nn.end());
        out_list.push_back(nn[0]);
    }

    return out_list;
}

std::vector<int> Lattice_v14::get_sites_for_wrapping_test(int s0_index) {
//# print("get_site_neighbor_of_site : ", s0_index)
    auto central_site = site_matrix[s0_index];
    int gid_central = central_site.get_gid();
    auto bonds = central_site.get_connecting_bonds(length);
//# print("bonds ", bonds)
    std::vector<int>  out_list, nn;
    for (auto bb : bonds) {
        nn = get_neighbor_sites(bb);
        for(int n=0;n < nn.size(); ++n){
            if (nn[n] == s0_index){
                nn.erase(nn.begin() + n);
                break;
            }
        }
        int gid = site_matrix[nn[0]].get_gid();
        if (gid == gid_central) {
            out_list.push_back(nn[0]);
        }
    }

    return out_list;
}

void Lattice_v14::print_bonds() {
    cout << "print_bonds" << endl;
    for(auto i: bond_ids) {
        auto bbonds = bond_matrix[i];
        auto c_sites = bbonds.get_connected_sites();
        if (c_sites.size() != 2) {
            cout << "warning : number of connected sites must be 2" << endl;
        }
        cout << "[" << i << "] gid=" << bbonds.get_gid() << " id="
             << bbonds.get_id() << " sites= {";
        for (auto c : c_sites) {
          cout << c << ",";
        }
        cout << "}" << endl;
    }
}

void Lattice_v14::print_sites(){
    cout << "print_sites" << endl;
    for(auto i: site_ids) {
        auto ssite = site_matrix[i];
        auto c_bonds = ssite.get_connecting_bonds(length);
        if (c_bonds.size() != 2) {
            cout << "warning : number of connected sites must be 2" << endl;
        }
        cout << "[" << i << "] gid=" << ssite.get_gid() << " id="
             << ssite.get_id() << " sites= {";
        for (auto c : c_bonds) {
            cout << c << ",";
        }
        cout << "}" << endl;
    }
}



Site_v14& Lattice_v14::get_site_by_id(int  id) {
#ifdef UNIT_TEST
    if (id >= site_matrix.size()){
        cout << "index out of range [" << id  << "] : " << __FILE__<< " : " << __LINE__  << endl;
        exit(-1);
    }
#endif
//    cout << "return site[" << id << "] => " << site_matrix[id].get_str(0) << endl;
    return site_matrix[id];
}

Bond_v14& Lattice_v14::get_bond_by_id(int id) {
    return bond_matrix[id];
}

void Lattice_v14::set_site_gid_by_id(int id, int gid){
    site_matrix[id].set_gid(gid);
}

void Lattice_v14::set_bond_gid_by_id(int id, int gid) {
//# print("id ", id)
    bond_matrix[id].set_gid(gid);
}

int Lattice_v14::get_site_gid_by_id(int id) {
#ifdef UNIT_TEST
    if(id >= site_matrix.size()){
        cerr << "Error. site id index out of range " << endl;
        cout << "id " << id << " length " << length << " site_matrix size " << site_matrix.size();
        exit(-1);
    }
#endif
    return site_matrix[id].get_gid();
}

int Lattice_v14::get_bond_gid_by_id(int id){
//# print("id ", id)
    return bond_matrix[id].get_gid();
}

Site_v14 Lattice_v14::get_site_by_index(size_t row, size_t col) {
    size_t s0_index = row * length + col;
    return site_matrix[s0_index];
}

/**
 * hv_flag:    0    means    horizontal and 1    means    vertical
 */
Bond_v14 Lattice_v14::get_bond_by_index(int row, int col, int hv_flag) {

    if (abs(hv_flag) > 1) {
        cout << "invalid hv flag" << endl;
    }
    int s0_index = row * length + col;
    s0_index += hv_flag * length * length;
    return bond_matrix[s0_index];
}


void Lattice_v14::init_relative_index(int id) {
    int s_index = site_ids[id];
    check_if_it_is_not_equal(id, s_index, __LINE__);
    site_matrix[s_index].init_relative_index();
}

void Lattice_v14::check_if_it_is_not_equal(int id, int s_index, int line) const {
    if (id != s_index) {
        cout << id << " != " << s_index << " : line " << line << endl;
    }
}

void Lattice_v14::set_relative_index(int id, RelativeIndex_v13 relative_index){
    int s_index = site_ids[id];
    check_if_it_is_not_equal(id, s_index, __LINE__);
    site_matrix[s_index].set_relative_index(relative_index);
}

vector<string> Lattice_v14::get_row_str(int row, int format) {
    vector<string> r_string;
    for(int cc =0; cc < length; ++cc){
        int s_index = row * length + cc;
        string str1 = site_matrix[s_index].get_str(format);
        string str2 = bond_matrix[s_index].get_str(format);
        stringstream str0;
        int tmp_l = str2.size() - str1.size();
        for (int i=0; i < tmp_l ; ++i){
            str0 << " ";
        }

        stringstream ss_main;
        ss_main << str1 << str0.str() << " " << str2 << " |";
        r_string.push_back(ss_main.str());
    }
    r_string.push_back("\n");
    return r_string;
}

vector<string> Lattice_v14::get_row_v_str(int row, int format) {
    vector<string> r_string;
    for(int cc =0; cc < length; ++cc){
        int s_index = row * length + cc;
        int v_bond_index = s_index + length * length;
        string str1 = site_matrix[s_index].get_str(format);
        string str3 = bond_matrix[v_bond_index].get_str(format);
        stringstream str4;
        for(auto s: str1){
            str4 << " ";
        }


        stringstream ss_main;
//        ss_main << str3 << "  " << str4.str() << " |";
        ss_main << str4.str()<< " " << str3  << " |";
        r_string.push_back(ss_main.str());
    }
    r_string.push_back("\n");
    return r_string;
}

void Lattice_v14::print_row_separator(int str_sz) {
    stringstream str_str;
    for(int i=0; i < str_sz; ++i){
        str_str << "-";
    }

    for(int cc =0; cc < length; ++ cc){
        cout << str_str.str();
    }

    cout << endl;
}

void Lattice_v14::view(int formatt) {
    if (formatt == 1){
        cout << "[gid, id]" << endl;
    } else if(formatt == 2){
        cout << "gid(row, col)" << endl;
    }

    cout << "format : " << formatt << endl;
    cout << "{site}           {horizontal bond}"<< endl;
    cout << "{vertical bond}  {               }"<< endl;
    cout << "The Lattice_v14 : " << endl;
    cout << "<--VIEW BEGIN-->"<< endl;
    print_row_separator(34);
    for(int rr=0; rr < length; ++rr) {
        auto aa = get_row_str(rr, formatt);
        for(auto a: aa){
            cout << a;
        }

        auto bb = get_row_v_str(rr, formatt);
        for(auto b: bb){
            cout << b;
        }
    }

    print_row_separator(34);

    cout << "<--VIEW END-->" << endl;

}


void Lattice_v14::view_relative_index() {

    cout << "Upper Left corner is <0,0>" << endl;
    cout << "<x,y> means relative index" << endl;
    cout << "90 degree Clockwise rotated Coordinate system" << endl;
    cout << "y value increases as we go rightward. Like columns" << endl;
    cout << "x value increases as we go downward . Like rows" << endl;
    cout << "Format 'gid<x,y>'" << endl;
    cout << "<--Relative index - VIEW BEGIN-->" << endl;
    int row_unit_str = 14;
    print_row_separator(row_unit_str);
    cout << setw(5) << right << "|";
    for(int cc=0; cc < length; ++cc) {

        cout << setw(6) << cc << setw(6) << right << "|";
    }
    cout << endl;
    print_row_separator(row_unit_str);
    for(int rr=0; rr < length; ++rr) {
        cout << setw(3) << rr << " |";
        for(int cc=0; cc < length; ++cc) {
            int s_index = rr * length + cc;
            auto site_s = site_matrix[s_index];
            auto a = site_s.get_relative_index();
            cout << setw(2) << site_s.get_gid() << a.get_str() << "|";
//# print("{:7}".foramt(a), end=' |')

        }
        cout << endl;
    }
    print_row_separator(row_unit_str);
    cout << "<--Relative index - VIEW END-->" <<endl;
}

void Lattice_v14::view_site_gids() {
    cout << "Upper Left corner is <0,0>"<< endl;
    cout << "<x,y> means relative index"<< endl;
    cout << "90 degree Clockwise rotated Coordinate system"<< endl;
    cout << "y value increases as we go rightward. Like columns"<< endl;
    cout << "x value increases as we go downward . Like rows"<< endl;
    cout << "Format 'gid<x,y>'"<< endl;
    cout << "<--Relative index - VIEW BEGIN-->"<< endl;

    int row_unit_str = 5;
    print_row_separator(row_unit_str);
    cout << setw(5) << right << "|";
    for(int cc=0; cc < length; ++cc) {

        cout << setw(3) << left << cc << right << "|";
    }
    cout << endl;
    print_row_separator(row_unit_str);

    for(int rr=0; rr < length; ++rr) {
        cout << setw(3) << rr << " |";
        for (int cc = 0; cc < length; ++cc) {
            int s_index = rr * length + cc;
            auto site_s = site_matrix[s_index];
//            RelativeIndex_v13 a = site_s.get_relative_index();
            cout << setw(3) << site_s.get_gid() << "|";
//# print("{:7}".foramt(a), end=' |')

        }
        cout << endl;
    }

    print_row_separator(row_unit_str);
    cout << "<--Relative index - VIEW END-->" << endl;

}

vector<int> Lattice_v14::get_all_neighbor_sites(int central_site_id) {
    auto bonds = get_neighbor_bonds(central_site_id);
    vector<int> four_neighbors;
    for(auto bid: bonds){
        auto connected_sites = get_neighbor_sites(bid);
        vector<int> tmp(connected_sites.begin(), connected_sites.end());
        if (tmp.size() != 2){
            cout << "func:get_all_neighbor_sites -> len(connected_sites) != 1" << endl;
//            cout << get_bond_by_id(bid) << endl;
            exit(-1);
        }
        for(auto aa: tmp){
            if(aa == central_site_id) continue;
            four_neighbors.push_back(aa);
        }

    }

    if(four_neighbors.size() != 4){
        cout << "Number of neighbors must be 4" << endl;
        exit(-1);
    }


    return four_neighbors;
}

void Lattice_v14::test_id_index_conversion() {
    int sq_length = length*length;
    for(int i=0; i < length; ++i){

        size_t rnd = rand() % sq_length;
        auto id1 = site_ids[rnd];
        check_if_it_is_not_equal(rnd, id1, __LINE__);
        auto id2 = site_matrix[id1].get_id();
        auto index = site_matrix[id1].get_index();
        auto id3 = index.row() * length + index.col();
        if(id1 != id2 || id1 != id3){
            cout << "Error. Site  " << endl;
            exit(-1);
        }

    }

    for(int i=0; i < length; ++i){

        size_t rnd = rand() % sq_length;
        auto id1 = bond_ids[rnd];
        auto id2 = bond_matrix[id1].get_id();
        auto rr = bond_matrix[id1].get_row();
        auto cc = bond_matrix[id1].get_col();
        auto tt = bond_matrix[id1].get_type();
        auto id3 = tt * sq_length + rr* length + cc;
        if(id1 != id2 || id1 != id3){
            cout << "Error. Bond  " << endl;
            exit(-1);
        }

    }

}

void Lattice_v14::viewLattice_clsuter_k(int cluster_gid) {
    cout << "Upper Left corner is <0,0>"<< endl;
    cout << "<x,y> means relative index"<< endl;
    cout << "90 degree Clockwise rotated Coordinate system"<< endl;
    cout << "y value increases as we go rightward. Like columns"<< endl;
    cout << "x value increases as we go downward . Like rows"<< endl;
    cout << "Format 'gid<x,y>'"<< endl;
    cout << "<--Relative index - VIEW BEGIN-->"<< endl;

    int row_unit_str = 5;
    print_row_separator(row_unit_str);
    cout << setw(5) << right << "|";
    for(int cc=0; cc < length; ++cc) {

        cout << setw(3) << left << cc << right << "|";
    }
    cout << endl;
    print_row_separator(row_unit_str);

    for(int rr=0; rr < length; ++rr) {
        cout << setw(3) << rr << " |";
        for (int cc = 0; cc < length; ++cc) {
            int s_index = rr * length + cc;
            auto site_s = site_matrix[s_index];
//            RelativeIndex_v13 a = site_s.get_relative_index();
            int gidvalue = site_s.get_gid();
            if (gidvalue != -1 && gidvalue != cluster_gid){
                gidvalue = -2;
            }
            cout << setw(3) << gidvalue << "|";
//# print("{:7}".foramt(a), end=' |')

        }
        cout << endl;
    }

    print_row_separator(row_unit_str);
    cout << "<--Relative index - VIEW END-->" << endl;


}

void Lattice_v14::print_site_bond_list() {
    for (auto ss=0; ss < site_matrix.size(); ss++){
        cout << "[" << setw(3) << ss << "]" << site_matrix[ss].get_gid() << site_matrix[ss].get_index().get_str()
             << ": {";
        for(auto i: site_matrix[ss].get_connecting_bonds(length)){
            cout << setw(3) << i << bond_matrix[i].get_str() << ",";
        }
        cout << "}" << endl;
    }

}
