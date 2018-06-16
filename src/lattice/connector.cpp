//
// Created by shahnoor on 11/9/2017.
//

#include "connector.h"

using namespace std;

// todo delete comment after February 2018

//
//Connector::Connector(Index site, value_type length, bool periodicity)
//        : _length{length}, _periodicity{periodicity}
//{
//    _site = site;
//    max_index = _length - 1;
//    calculate();
//}
//
//std::vector<BondIndex> Connector::connectingBonds(bool peridicity) {
//    vector<BondIndex> bi;
//    value_type prev_column = (_site.y_ + _length - 1) % _length;
//    value_type prev_row = (_site.x_ + _length - 1) % _length;
//    if(peridicity){
//        bi.resize(4);
//        bi[0] = {BondType::Horizontal,  _site.x_, _site.y_};
//        bi[1] = {BondType::Vertical,    _site.x_, _site.y_};
//        bi[2] = {BondType::Horizontal, _site.x_, prev_column};
//        bi[3] = {BondType::Vertical, prev_row, _site.y_};
//    }
//    else{
//        // without periodicity
//
//        if (_site.x_ == min_index && _site.y_ == min_index) { // upper left corner
//            bi.resize(2);
//
//            bi[0] = {BondType::Horizontal, _site.x_, _site.y_};
//            bi[1] = {BondType::Vertical, _site.x_, _site.y_};
//        } else if (_site.x_ == min_index && _site.y_ == max_index) {// upper right corner
//            bi.resize(2);
//
//            bi[0] = {BondType::Horizontal, _site.x_, prev_column};
//            bi[1] = {BondType::Vertical, _site.x_, _site.y_};
//        } else if (_site.x_ == max_index && _site.y_ == max_index) { // lower right
//            bi.resize(2);
//
//            // lower right corner
//            bi[0] = {BondType::Horizontal, _site.x_, prev_column};
//            bi[1] = {BondType::Vertical, prev_row, _site.y_ };
//        } else if (_site.x_ == max_index && _site.y_ == min_index) { // lower left
//            bi.resize(2);
//
//            // lower left corner
//            bi[0] = {BondType::Horizontal, _site.x_, _site.y_};
//            bi[1] = {BondType::Vertical, prev_row, _site.y_ };
//        }
//            // now the edges
//        else if (_site.x_ == min_index && _site.y_ > min_index && _site.y_ < max_index) {
//            // upper line but not in the upper corners. top edge
//            bi.resize(3);
//
//            bi[0] = {BondType::Horizontal, _site.x_, _site.y_};
//            bi[1] = {BondType::Vertical, _site.x_, _site.y_};
//            bi[2] = {BondType::Horizontal, _site.x_, prev_column};
//        } else if (_site.x_ == max_index && _site.y_ > min_index && _site.y_ < max_index) { // ok
//            //  bottom edge
//            bi.resize(3);
//
//            bi[0] = {BondType::Horizontal, _site.x_, _site.y_};
//            bi[1] = {BondType::Vertical, prev_row, _site.y_ };
//            bi[2] = {BondType::Horizontal, _site.x_, prev_column};
//        } else if (_site.x_ > min_index && _site.x_ < max_index && _site.y_ == min_index) { // left edge
//            bi.resize(3);
//
//            // lower line not in the corners
//            bi[0] = {BondType::Horizontal, _site.x_, _site.y_};
//            bi[1] = {BondType::Vertical, _site.x_, _site.y_};
//            bi[2] = {BondType::Vertical, prev_row, _site.y_ };
//
//        } else if (_site.x_ > min_index && _site.x_ < max_index && _site.y_ == max_index) {
//            // right edge
//            bi.resize(3);
//
//            // left line without corners
//            bi[0] = {BondType::Horizontal, _site.x_, prev_column};
//            bi[1] = {BondType::Vertical, _site.x_, _site.y_};
//            bi[2] = {BondType::Vertical, prev_row, _site.y_ };
//        }
//        else {
//            // 1 level inside the lattice
//            // not in any the boundary
//            bi.resize(4);
//            bi[0] = {BondType::Horizontal,  _site.x_, _site.y_};
//            bi[0] = {BondType::Vertical,    _site.x_, _site.y_};
//            bi[0] = {BondType::Horizontal, _site.x_, prev_column};
//            bi[0] = {BondType::Vertical, prev_row, _site.y_};
//        }
//    }
//    return bi;
//}
//
//
//std::vector<Index> Connector::nearestNeighbors(bool peridicity) {
//    std::vector<Index> neighbors;
//    if(peridicity){
//        neighbors.resize(4);
//        neighbors[0] = {_site.x_, (_site.y_ + 1) % _length};
//        neighbors[1] = {_site.x_, (_site.y_ + _length - 1) % _length};
//        neighbors[2] = {(_site.x_ + 1) % _length, _site.y_};
//        neighbors[3] = {(_site.x_ + _length - 1) % _length, _site.y_};
//    }
//    else{
//        // without periodicity
//        if (_site.x_ == min_index && _site.y_ == min_index) { // upper left corner
//            neighbors.resize(2);
//            neighbors[0] = {_site.x_, (_site.y_ + 1) % _length};
//            neighbors[1] = {(_site.x_ + 1) % _length, _site.y_};
//
//        } else if (_site.x_ == min_index && _site.y_ == max_index) {// upper right corner
//            neighbors.resize(2);
//            neighbors[0] = {_site.x_, (_site.y_ + _length - 1) % _length};
//            neighbors[1] = {(_site.x_ + 1) % _length, _site.y_};
//
//        } else if (_site.x_ == max_index && _site.y_ == max_index) { // lower right
//            neighbors.resize(2);
//            // lower right corner
//            neighbors[0] = {_site.x_, (_site.y_ + _length - 1) % _length};
//            neighbors[1] = {(_site.x_ + _length - 1) % _length, _site.y_};
//
//        } else if (_site.x_ == max_index && _site.y_ == min_index) { // lower left
//            neighbors.resize(2);
//            // lower left corner
//            neighbors[0] = {_site.x_, (_site.y_ + 1) % _length};
//            neighbors[1] = {(_site.x_ + 1) % _length, _site.y_};
//
//        }
//            // now the edges
//        else if (_site.x_ == min_index && _site.y_ > min_index && _site.y_ < max_index) {
//            // upper line but not in the upper corners. top edge
//            neighbors.resize(3);
//            neighbors[0] = {_site.x_, (_site.y_ + 1) % _length};
//            neighbors[1] = {_site.x_, (_site.y_ + _length - 1) % _length};
//            neighbors[2] = {(_site.x_ + 1) % _length, _site.y_};
//
//        } else if (_site.x_ == max_index && _site.y_ > min_index && _site.y_ < max_index) { // ok
//            //  bottom edge
//            neighbors.resize(3);
//            neighbors[0] = {_site.x_, (_site.y_ + 1) % _length};
//            neighbors[1] = {_site.x_, (_site.y_ + _length - 1) % _length};
//            neighbors[2] = {(_site.x_ + _length - 1) % _length, _site.y_};
//
//        } else if (_site.x_ > min_index && _site.x_ < max_index && _site.y_ == min_index) {
//            neighbors.resize(3);
//            // left edge
//            // lower line not in the corners
//            neighbors[0] = {_site.x_, (_site.y_ + 1) % _length};
//            neighbors[1] = {(_site.x_ + 1) % _length, _site.y_};
//            neighbors[2] = {(_site.x_ + _length - 1) % _length, _site.y_};
//
//
//        } else if (_site.x_ > min_index && _site.x_ < max_index && _site.y_ == max_index) {
//            // right edge
//            neighbors.resize(3);
//            // left line without corners
//            // left edge
//            neighbors[0] = {_site.x_, (_site.y_ + _length - 1) % _length};
//            neighbors[1] = {(_site.x_ + 1) % _length, _site.y_};
//            neighbors[2] = {(_site.x_ + _length - 1) % _length, _site.y_};
//        }
//        else{
//            neighbors.resize(4);
//            neighbors[0] = {_site.x_, (_site.y_ + 1) % _length};
//            neighbors[1] = {_site.x_, (_site.y_ + _length - 1) % _length};
//            neighbors[2] = {(_site.x_ + 1) % _length, _site.y_};
//            neighbors[3] = {(_site.x_ + _length - 1) % _length, _site.y_};
//        }
//    }
//    return neighbors;
//}
//
//
//void Connector::calculate() {
//    value_type prev_column = (_site.y_ + _length - 1) % _length;
//    value_type prev_row = (_site.x_ + _length - 1) % _length;
//    value_type next_row = (_site.x_ + 1) % _length;
//    value_type next_column = (_site.y_ + 1) % _length;
//    if(_periodicity){
//        _bond_index.resize(4);
//        _bond_index[0] = {BondType::Horizontal,  _site.x_, _site.y_};
//        _bond_index[1] = {BondType::Vertical,    _site.x_, _site.y_};
//        _bond_index[2] = {BondType::Horizontal, _site.x_, prev_column};
//        _bond_index[3] = {BondType::Vertical, prev_row, _site.y_};
//
//        _neighbor_index.resize(4);
//        _neighbor_index[0] = {_site.x_, next_column};
//        _neighbor_index[1] = {_site.x_, prev_column};
//        _neighbor_index[2] = {next_row, _site.y_};
//        _neighbor_index[3] = {prev_row, _site.y_};
//    }
//    else{
//        // without periodicity
//
//        if (_site.x_ == min_index && _site.y_ == min_index) { // upper left corner
//            _bond_index.resize(2);
//
//            _bond_index[0] = {BondType::Horizontal, _site.x_, _site.y_};
//            _bond_index[1] = {BondType::Vertical, _site.x_, _site.y_};
//
//            _neighbor_index.resize(2);
//            _neighbor_index[0] = {_site.x_, next_column};
//            _neighbor_index[1] = {next_row, _site.y_};
//
//        } else if (_site.x_ == min_index && _site.y_ == max_index) {// upper right corner
//            _bond_index.resize(2);
//
//            _bond_index[0] = {BondType::Horizontal, _site.x_, prev_column};
//            _bond_index[1] = {BondType::Vertical, _site.x_, _site.y_};
//
//            _neighbor_index.resize(2);
//            _neighbor_index[0] = {_site.x_, prev_column};
//            _neighbor_index[1] = {next_row, _site.y_};
//        } else if (_site.x_ == max_index && _site.y_ == max_index) { // lower right
//            _bond_index.resize(2);
//
//            // lower right corner
//            _bond_index[0] = {BondType::Horizontal, _site.x_, prev_column};
//            _bond_index[1] = {BondType::Vertical, prev_row, _site.y_ };
//
//            _neighbor_index.resize(2);
//            // lower right corner
//            _neighbor_index[0] = {_site.x_, prev_column};
//            _neighbor_index[1] = {prev_row, _site.y_};
//
//        } else if (_site.x_ == max_index && _site.y_ == min_index) { // lower left
//            _bond_index.resize(2);
//
//            // lower left corner
//            _bond_index[0] = {BondType::Horizontal, _site.x_, _site.y_};
//            _bond_index[1] = {BondType::Vertical, prev_row, _site.y_ };
//
//            _neighbor_index.resize(2);
//            // lower left corner
//            _neighbor_index[0] = {_site.x_, next_column};
//            _neighbor_index[1] = {next_row, _site.y_};
//        }
//            // now the edges
//        else if (_site.x_ == min_index && _site.y_ > min_index && _site.y_ < max_index) {
//            // upper line but not in the upper corners. top edge
//            _bond_index.resize(3);
//
//            _bond_index[0] = {BondType::Horizontal, _site.x_, _site.y_};
//            _bond_index[1] = {BondType::Vertical, _site.x_, _site.y_};
//            _bond_index[2] = {BondType::Horizontal, _site.x_, prev_column};
//
//            _neighbor_index.resize(3);
//            _neighbor_index[0] = {_site.x_, next_column};
//            _neighbor_index[1] = {_site.x_, prev_column};
//            _neighbor_index[2] = {next_row, _site.y_};
//        } else if (_site.x_ == max_index && _site.y_ > min_index && _site.y_ < max_index) { // ok
//            //  bottom edge
//            _bond_index.resize(3);
//
//            _bond_index[0] = {BondType::Horizontal, _site.x_, _site.y_};
//            _bond_index[1] = {BondType::Vertical, prev_row, _site.y_ };
//            _bond_index[2] = {BondType::Horizontal, _site.x_, prev_column};
//
//            _neighbor_index.resize(3);
//            _neighbor_index[0] = {_site.x_, next_column};
//            _neighbor_index[1] = {_site.x_, prev_column};
//            _neighbor_index[2] = {prev_row, _site.y_};
//        } else if (_site.x_ > min_index && _site.x_ < max_index && _site.y_ == min_index) { // left edge
//            _bond_index.resize(3);
//
//            // lower line not in the corners
//            _bond_index[0] = {BondType::Horizontal, _site.x_, _site.y_};
//            _bond_index[1] = {BondType::Vertical, _site.x_, _site.y_};
//            _bond_index[2] = {BondType::Vertical, prev_row, _site.y_ };
//
//            _neighbor_index.resize(3);
//            // left edge
//            // lower line not in the corners
//            _neighbor_index[0] = {_site.x_, next_column};
//            _neighbor_index[1] = {next_row, _site.y_};
//            _neighbor_index[2] = {prev_row, _site.y_};
//
//        } else if (_site.x_ > min_index && _site.x_ < max_index && _site.y_ == max_index) {
//            // right edge
//            _bond_index.resize(3);
//
//            // left line without corners
//            _bond_index[0] = {BondType::Horizontal, _site.x_, prev_column};
//            _bond_index[1] = {BondType::Vertical, _site.x_, _site.y_};
//            _bond_index[2] = {BondType::Vertical, prev_row, _site.y_ };
//
//            _neighbor_index.resize(3);
//            // left line without corners
//            // left edge
//            _neighbor_index[0] = {_site.x_, prev_column};
//            _neighbor_index[1] = {next_row, _site.y_};
//            _neighbor_index[2] = {prev_row, _site.y_};
//        }
//        else {
//            // 1 level inside the lattice
//            // not in any the boundary
//            _bond_index.resize(4);
//            _bond_index[0] = {BondType::Horizontal,  _site.x_, _site.y_};
//            _bond_index[0] = {BondType::Vertical,    _site.x_, _site.y_};
//            _bond_index[0] = {BondType::Horizontal, _site.x_, prev_column};
//            _bond_index[0] = {BondType::Vertical, prev_row, _site.y_};
//
//            _neighbor_index.resize(4);
//            _neighbor_index[0] = {_site.x_, next_column};
//            _neighbor_index[1] = {_site.x_, prev_column};
//            _neighbor_index[2] = {next_row, _site.y_};
//            _neighbor_index[3] = {prev_row, _site.y_};
//        }
//    }
//}




