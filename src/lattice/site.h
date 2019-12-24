//
// Created by shahnoor on 10/2/2017.
//

#ifndef SITEPERCOLATION_SITE_H
#define SITEPERCOLATION_SITE_H

#include <array>
#include <set>
#include <vector>
#include <iostream>
#include <memory>
#include <strings.h>

#include "../index/index.h"
#include "../types.h"
#include "../flags.h"



/**
 * single Site of a lattice
 */
struct Site{

    /**
     * if true -> site is placed.
     * if false -> the (empty) position is there but the site is not (required for site percolation)
     */
    bool _status{false};
    int _group_id{-1};
    Index _id{};
    value_type _weight{};    // when one site is occupied multiple times only weight increases by occupation

    //relative distance from the root site. {0,0} if it is the root site
    //very useful for detecting wrapping
    IndexRelative _relative_index{0,0};

//    std::vector<Index> _neighbor{};
    // _connection[0] -> _connection between "set_ID" and "_neighbor[0]"

    /**
     * using set gives following advantages
     * 1. Keeps elements sorted
     * 2. automatically manage repeated value
     */
//    std::set<int> _connection{};

public:

    ~Site()                 = default;
    Site()                  = default;
    Site(const Site&)             = default;
    Site(Site&&)            = default;
    Site& operator=(const Site&)  = default;
    Site& operator=(Site&&) = default;

    Site(Index id, value_type length){
        // I have handle _neighbor or corner points and edge points carefully
        if(id.row_ >= length || id.column_ >= length){
            std::cout << "out of range : line " << __LINE__ << std::endl;
        }
        _id.row_ = id.row_;
        _id.column_ = id.column_;
//        _neighbor = Lattice_Helper::get_neighbor(id, length);
    }


    bool isActive() const { return _status;}
    void activate(){ _status = true;}
    void deactivate() {
        _relative_index = {0,0};
        _group_id = -1;
        _status = false;
    }
    Index ID() const { return  _id;}
    /*
     * Group get_ID is the set_ID of the cluster they are in
     */
    int     get_groupID() const {return _group_id;}
    void    set_groupID(int g_id) {_group_id = g_id;}

    std::stringstream getSite() const {
        std::stringstream ss;
        if(isActive())
            ss << _id;
        else
            ss << "(*)";
        return ss;
    }

    value_type weight() const {return _weight;}
    void increaseWeight(){_weight += 1;}

    void relativeIndex(IndexRelative r){
        _relative_index = r;
    }

    void relativeIndex(int x, int y){
        _relative_index = {x,y};
    }

    IndexRelative relativeIndex() const {return _relative_index;}

};

std::ostream& operator<<(std::ostream& os, const Site& site);
bool operator==(Site& site1, Site& site2);

/*************************
 * Date : 2019.12.23
 * Author : Shahnoor
 *
 * Specially build for v12
 */
class Site_v12{
    /*
     * properties that can be assigned during construction only.
     */
    int _id{-1};
    Index _index{};


    /*
     * Properties that can be assigned later
     */
    int _group_id{-1}; // id of the cluster of which it is a member
    bool _status{false};

    //relative distance from the root site. {0,0} if it is the root site
    //very useful for detecting wrapping
    IndexRelative _relative_index{0,0}; // relative index from root Site of the cluster that it belongs


public:

    ~Site_v12() = default;
    Site_v12() = default;
    Site_v12(Index index, int id){
        _index.row_ = index.row_;
        _index.column_ = index.column_;
        _id = id;
    }


    bool check(int length){
        // I have handle _neighbor or corner points and edge points carefully
        if(_index.row_ >= length || _index.column_ >= length){
            std::cout << "out of range : line " << __LINE__ << std::endl;
            return false;
        }
        return true;
    }
    bool isActive() const { return _status;}
    void activate(){ _status = true;}
    void deactivate() {
        _relative_index = {0,0};
        _group_id = -1;
        _status = false;
    }
    int get_id() const { return  _id;}
    Index get_index() const {return _index;}

    int     get_groupID() const {return _group_id;}
    void    set_groupID(int g_id) {_group_id = g_id;}


    void relativeIndex(IndexRelative r){
        _relative_index = r;
    }

    void relativeIndex(int x, int y){
        _relative_index = {x,y};
    }

    IndexRelative relativeIndex() const {return _relative_index;}

};

std::ostream& operator<<(std::ostream& os, const Site_v12& site);

#endif //SITEPERCOLATION_SITE_H
