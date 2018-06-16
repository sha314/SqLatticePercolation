//
// Created by shahnoor on 10/28/17.
//

#include "time_tracking.h"

#include <iostream>
#include <sstream>

using namespace std;

/**
 *
 * @param t -> Time in seconds
 * @return
 */
std::string getFormattedTime(double t){
    std::string s;

    size_t hr{}, min{};
    double sec{}; // hour, min, sec
    size_t integer_sec = size_t(t);

    if(integer_sec >= 3600)
        hr = integer_sec / 3600;
    if(integer_sec >= 60)
        min = (integer_sec - hr * 3600) / 60;
    sec = t -  min * 60 - hr * 3600;

    s += to_string(hr);
    s += " hr ";
    s += to_string(min);
    s += " min ";
    s += to_string(sec);
    s += " sec ";
    return s;
}

std::string currentTime() {
    time_t t = time(0);   // get birthTime now
    struct tm * now = localtime( & t );
    stringstream ss;
    ss << (now->tm_year + 1900) << '.'
         << (now->tm_mon + 1) << '.'
         <<  now->tm_mday << '_'
         << now->tm_hour << '.' << now->tm_min << '.' << now->tm_sec;

    return ss.str();
}


