//
// Created by shahnoor on 11/26/2017.
//

#include "fractal_tool.h"

using namespace std;

/**
 *
 * @param length
 * @return
 */
vector<value_type> get_delta(value_type length){
    vector<value_type> v;
    for(value_type i{1}; i != length; ++i){
        if( length % i == 0){
            v.push_back(i);
        }
    }
    return v;
}


/**
 *
 * @param length
 */
Fractal_2d::Fractal_2d(value_type length)
        : _length{length}
{
    if(_length % 2 != 0){
        cout << "Odd length error. Cannot divide_in_four_parts odd length" << endl;
    }


    middle.reserve(15);
    value_type x{_length};
    while (x % 2 != 0){
        middle.push_back(x);
        x = _length/2;
    }
    middle.shrink_to_fit();
    rectangles.reserve(30);
    rectangles.push_back({0, 0, length, length});
    _delta = length;
    grid = SqLattice(_length);
    grid.activateAllSite(); // initially
}



/**
 *
 * x0 \ x0   x1  x2
 * .        .   .
 * x1       .   .
 * .        .   .
 * x2       .   .
 *
 * @param min
 * @param max
 */
void Fractal_2d::divide_in_four_parts(value_type min, value_type max) {
    _delta /= 2;
    value_type dx = max / 2;
    value_type x0 = min;
    value_type x1 = min + dx;
    value_type x2 = max;
    // (x0,x0)  to (x1, x1) -> deactivation, id = -1
    for(value_type i{x0}; i != x1; ++i){
        for(value_type j{x0}; j != x1; ++j){
            // upper left
            grid.getSite({i, j}).deactivate();
        }
    }
    // (x0, x1) to (x1, x2) -> activation, id = new id
    ++id;
    for(value_type i{x0}; i != x1; ++i){
        for(value_type j{x1}; j != x2; ++j){
            // upper left
            grid.getSite({i, j}).set_groupID(id);
        }
    }

    // (x1, x0) to (x2, x1) -> activation, id = new id
    ++id;
    for(value_type i{x1}; i != x2; ++i){
        for(value_type j{x0}; j != x1; ++j){
            // upper left
            grid.getSite({i, j}).set_groupID(id);
        }
    }
    // (x1, x1) to (x2, x2) -> activation, id = new id
    ++id;
    for(value_type i{x1}; i != x2; ++i){
        for(value_type j{x1}; j != x2; ++j){
            // upper left
            grid.getSite({i, j}).set_groupID(id);
        }
    }

    rectangles.push_back({x0, x1, x1, x2});   // store reusable or active rectangles
    rectangles.push_back({x1, x0, x2, x1});
    rectangles.push_back({x1, x1, x2, x2});

}



/**
 *
 * r0 \ c0   c1  c2
 * .        .   .
 * r1       .   .
 * .        .   .
 * r2       .   .
 *
 * @param min
 * @param max
 */
void Fractal_2d::divide_in_four_parts(const Rectangle& rectangle) {
    _delta /= 2;
    value_type dr = rectangle.delta_row() / 2;
    value_type dc = rectangle.delta_col() / 2;

    value_type c0 = rectangle.col_i;
    value_type c1 = rectangle.col_i + dc;
    value_type c2 = rectangle.col_j;

    value_type r0 = rectangle.row_i;
    value_type r1 = rectangle.row_i + dr;
    value_type r2 = rectangle.row_j;

    // (x0,x0)  to (x1, x1) -> deactivation, id = -1
    for(value_type i{r0}; i != r1; ++i){
        for(value_type j{c0}; j != c1; ++j){
            // upper left
            grid.getSite({i, j}).deactivate();
        }
    }
    // (x0, x1) to (x1, x2) -> activation, id = new id
    ++id;
    for(value_type i{r0}; i != r1; ++i){
        for(value_type j{c1}; j != c2; ++j){
            // upper left
            grid.getSite({i, j}).set_groupID(id);
        }
    }

    // (x1, x0) to (x2, x1) -> activation, id = new id
    ++id;
    for(value_type i{r1}; i != r2; ++i){
        for(value_type j{c0}; j != c1; ++j){
            // upper left
            grid.getSite({i, j}).set_groupID(id);
        }
    }
    // (x1, x1) to (x2, x2) -> activation, id = new id
    ++id;
    for(value_type i{r1}; i != r2; ++i){
        for(value_type j{c1}; j != c2; ++j){
            // upper left
            grid.getSite({i, j}).set_groupID(id);
        }
    }

    rectangles.push_back({r0, c1, r1, c2});   // store reusable or active rectangles
    rectangles.push_back({r1, c0, r2, c1});
    rectangles.push_back({r1, c1, r2, c2});

}




/**
 * this method will call #divide_in_four_parts() method with suitable argument from previous run.
 *
 */
void Fractal_2d::divide() {
//    if(middle.empty()){
//        cout << "Cannot divide further" << endl;
//        return;
//    }
//    auto it = middle.begin();

    vector<Rectangle> tmp(rectangles);
    rectangles.clear();
    rectangles.reserve(tmp.size()*3);
    for(auto a: tmp){
        // divide all the  rectangles and remove the current rectangle.
        // new three reusable rectangle will be added automatically in divide_in_four_parts() function
        divide_in_four_parts(a);
    }



//    middle.erase(it);   // erase used value
}


/**
 *
 */
void Fractal_2d::view()  {
    cout << "delta " << _delta << endl;
    grid.view_sites_extended();
}


/********************************
 *
 *
 *******************************/

void boxCounting::setValueToLookFor(int val){
    _first_spanning_cluster_id = val;
}

/**
* Box counting method for fractal dimension calculation
* @param delta
* @return
*/
value_type boxCounting::box_counting(const SqLattice& lattice, value_type delta) {
    auto _length = lattice.length();
    value_type counter{};
    if(_length % delta == 0){
        for(value_type r{} ; r < _length ; r += delta){
            for(value_type c{}; c < _length ; c += delta){
                if(anyActiveSite(lattice, r, c, delta)){
                    ++ counter;
                }
            }
        }
    }
    else{
        cout << "Delta size is not corrent" << endl;
    }

    return counter;
}


/**
 * Return true if there is at least one active site in the region r, r + delta and c, c + delta
 * @param r
 * @param c
 * @param delta
 * @return
 */
bool boxCounting::anyActiveSite(const SqLattice& lattice, value_type row, value_type col, value_type delta) {
    for(value_type r{} ; r < delta ; ++r) {
        for (value_type c{}; c < delta; ++c) {
            if(lattice.getSite({r+row, c + col}).isActive()){
                return true;
            }
        }
    }
    return false;
}


/**
 * Return true if there is at least one active site in the region r, r + delta and c, c + delta
 * Also this function required that, the site is in the spanning cluster
 * @param row
 * @param col
 * @param delta
 * @return
 */
bool boxCounting::anyActiveSpanningSite(const SqLattice& lattice, value_type row, value_type col, value_type delta) {
    for(value_type r{} ; r < delta ; ++r) {
        for (value_type c{}; c < delta; ++c) {
            if(lattice.getSite({r + row, c + col}).get_groupID() == _first_spanning_cluster_id){
                return true;
            }
        }
    }
    return false;
}

