//
// Created by shahnoor on 11/26/2017.
//

#ifndef PERCOLATION_FRACTAL_TOOL_H
#define PERCOLATION_FRACTAL_TOOL_H


#include "../lattice/lattice.h"

#include <iostream>


std::vector<value_type> get_delta(value_type length);



/**
 * Contains index of initial and final point of a rectangle
 *
 * (i,i)    ...     (i,j)
 * .
 * .
 * .
 * (j,i)    ...     (j,j)
 *
 *
 */
struct Rectangle{
    value_type row_i, row_j;
    value_type col_i, col_j;
public:
    Rectangle(value_type r1, value_type c1, value_type r2, value_type c2)
            : row_i{r1}, col_i{c1}, row_j{r2}, col_j{c2} {}

    value_type delta_row() const {
        return row_j > row_i ? row_j - row_i : row_i - row_j;
    }

    value_type delta_col() const {
        return col_j > col_i ? col_j - col_i : col_i - col_j;
    }

};


/***
 * todo rename to a better name
 * Sierpinski Carpet
 */
class Fractal_2d{
    SqLattice grid;
    value_type _length{};
    int id{};
    std::vector<value_type> middle;
    std::vector<Rectangle> rectangles;
    value_type _delta;


public:
    ~Fractal_2d() = default;
//    Fractal_2d() = default;
    Fractal_2d(value_type length);
    value_type length() const { return _length;}
    void view();

    void suitable_lengths(){
        std::cout << "Program can continue dividing as long as length is even number" << std::endl;
        std::cout << "Length -> generators" << std::endl;
        std::cout << "80 -> 40, 20, 10, 5" << std::endl;
        std::cout << "88 -> 44, 22, 11" << std::endl;
        std::cout << "96 -> 48, 24, 12, 6, 3" << std::endl;
        std::cout << "160 -> 80, 40, 20, 10, 5" << std::endl;
        std::cout << "176 -> 88, 44, 22, 11" << std::endl;
        std::cout << "192 -> 96, 48, 24, 12, 6, 3" << std::endl;
    }

    void divide();
    void divide_in_four_parts(value_type min, value_type max);  // divide_in_four_parts into four equal parts and deletes the upper left corner
    void divide_in_four_parts(const Rectangle& rectangle);    // todo on test_spanning


};


/**
 *
 */
struct boxCounting{
    int _first_spanning_cluster_id;

    void setValueToLookFor(int val);

    value_type box_counting(const SqLattice& lattice, value_type delta);

    bool anyActiveSite(const SqLattice& lattice, value_type row, value_type col, value_type delta);
    bool anyActiveSpanningSite(const SqLattice& lattice, value_type row, value_type col, value_type delta);

};



#endif //PERCOLATION_FRACTAL_TOOL_H
