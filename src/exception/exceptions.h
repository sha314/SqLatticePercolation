//
// Created by shahnoor on 10/12/2017.
//

#ifndef PERCOLATION_EXCEPTIONS_H
#define PERCOLATION_EXCEPTIONS_H


#include <string>
#include <iostream>

/**
 *
 */
struct Mismatch{
    std::string msg_;
    size_t line_;
    Mismatch(size_t line, std::string msg="")
            :line_{line}, msg_{msg}  {}

    void what() const {
        std::cerr << msg_ << "\nId and index mismatch at line " << line_ << std::endl;
    }
};

/**
 *
 */
struct InvalidIndex{
    std::string msg_;
    InvalidIndex(std::string msg)  :msg_{msg}  {}

    void what() const {
        std::cerr << msg_ << std::endl;
    }
};

/**
 *
 */
struct InvalidBond{
    std::string msg_;
    InvalidBond(std::string msg)  :msg_{msg}  {}

    void what() const {
        std::cout << msg_ << std::endl;
    }
};


/**
 * When any neighor is occupied and no suitable neighbor is found, throw this exception
 */
struct OccupiedNeighbor{
    std::string msg_;
    OccupiedNeighbor(std::string msg): msg_{msg}{}

    void what() const {
        std::cout << msg_ << std::endl;
    }
};




#endif //PERCOLATION_EXCEPTIONS_H
