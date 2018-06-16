//
// Created by shahnoor on 10/25/17.
//

#include "latticeConfiguration.h"

using namespace std;

std::vector<Index> LatticeSiteConfiguration::configuration_6x6(int x){
    switch (x){
        case 0:
            return configuration_6x6_0();
        case 1:
            return configuration_6x6_1();
        case 2:
            return configuration_6x6_2();
        case 3:
            return configuration_6x6_3();
        case 4:
            return configuration_6x6_4();
        case 5:
            return configuration_6x6_5();
        default:
            cout << "function argument mismatched. using default" << endl;
            return configuration_6x6_4();
    }
}

vector<Index> LatticeSiteConfiguration::configuration_6x6_1(){
    cout << "configuration_6x6_1() : line " << __LINE__ << endl;
    vector<Index> custome_indices;
    custome_indices.push_back({0,2});
//    custome_indices.push_back({0,3});

    custome_indices.push_back({1,1});
    custome_indices.push_back({1,2});
    custome_indices.push_back({1,3});

    custome_indices.push_back({2,1});

    custome_indices.push_back({3,0});
    custome_indices.push_back({3,1});
    custome_indices.push_back({3,2});
    custome_indices.push_back({3,4});
    custome_indices.push_back({3,5});///

    custome_indices.push_back({4,2});
    custome_indices.push_back({4,3});
    custome_indices.push_back({4,5});

    custome_indices.push_back({5,2});
    custome_indices.push_back({5,3});

    return custome_indices;
}



vector<Index> LatticeSiteConfiguration::configuration_6x6_2(){
    cout << "configuration_6x6_2() : line " << __LINE__ << endl;
    vector<Index> custome_indices;

    custome_indices.push_back({1,0});
    custome_indices.push_back({1,1});
    custome_indices.push_back({1,2});

    custome_indices.push_back({2,2});

    custome_indices.push_back({3,2});
    custome_indices.push_back({3,3});
    custome_indices.push_back({3,4});


    custome_indices.push_back({4,4});

    custome_indices.push_back({5,3});
    custome_indices.push_back({5,4});

    return custome_indices;
}

/**
 * 1 big cluster and 1 vertical spanning
 * no horizontal spanning
 * @return
 */
vector<Index> LatticeSiteConfiguration::configuration_6x6_3(){
    cout << "configuration_6x6_3() : line " << __LINE__ << endl;
    vector<Index> custome_indices;

    custome_indices.push_back({0,2});
    custome_indices.push_back({0,4});

    custome_indices.push_back({1,0});
    custome_indices.push_back({1,1});
    custome_indices.push_back({1,2});
    custome_indices.push_back({1,4});

    custome_indices.push_back({2,0});
    custome_indices.push_back({2,3});
    custome_indices.push_back({2,4});

    custome_indices.push_back({3,0});
    custome_indices.push_back({3,3});
    custome_indices.push_back({3,5});


    custome_indices.push_back({4,2});
    custome_indices.push_back({4,3});
    custome_indices.push_back({4,5});

    custome_indices.push_back({5,2});
    custome_indices.push_back({5,4});
    custome_indices.push_back({5,5});

    return custome_indices;
}



/**
 * 3 cluster and 1 vertical spanning
 * @return
 */
vector<Index> LatticeSiteConfiguration::configuration_6x6_4(){
    cout << "configuration_6x6_4() : line " << __LINE__ << endl;
    vector<Index> custome_indices;

    custome_indices.push_back({0,4});

    custome_indices.push_back({1,0});
    custome_indices.push_back({1,1});
    custome_indices.push_back({1,2});
    custome_indices.push_back({1,4});

    custome_indices.push_back({2,0});
    custome_indices.push_back({2,3});
    custome_indices.push_back({2,4});


    custome_indices.push_back({3,3});
    custome_indices.push_back({3,5});


    custome_indices.push_back({4,2});
    custome_indices.push_back({4,3});
    custome_indices.push_back({4,5});

    custome_indices.push_back({5,2});
    custome_indices.push_back({5,5});

    return custome_indices;
}




/**
 * Cluster Model 6x6-5.png
 *
 * 1 vertical spanning
 * Cluster travles from to to middle then again upward and then again downward
 * @return
 */
vector<Index> LatticeSiteConfiguration::configuration_6x6_5(){
    cout << "configuration_6x6_5() : line " << __LINE__ << endl;
    vector<Index> custome_indices;

    custome_indices.push_back({0,4});

    custome_indices.push_back({1,0});
    custome_indices.push_back({1,1});
    custome_indices.push_back({1,2});
    custome_indices.push_back({1,4});

    custome_indices.push_back({2,0});
    custome_indices.push_back({2,2});
    custome_indices.push_back({2,4});


    custome_indices.push_back({3,0});
    custome_indices.push_back({3,2});
    custome_indices.push_back({3,3});
    custome_indices.push_back({3,4});


    custome_indices.push_back({4,0});


    custome_indices.push_back({5,0});
    custome_indices.push_back({5,1});

    return custome_indices;
}

/***
 * It is not any particular model
 */
std::vector<Index> LatticeSiteConfiguration::configuration_6x6_0() {
    cout << "configuration_6x6_5() : line " << __LINE__ << endl;
    vector<Index> custome_indices;

    // corners
    custome_indices.push_back({0,0});
    custome_indices.push_back({0,1});
    custome_indices.push_back({0,2});

    custome_indices.push_back({1,1});


    return custome_indices;
}

/*******************************************************
 *
 *
 *******************************************************/


/**
 *
 * @param x
 * @return
 */
std::vector<BondIndex> LatticeBondConfiguration::configuration_6x6(int x) {




    return vector<BondIndex>();
}

