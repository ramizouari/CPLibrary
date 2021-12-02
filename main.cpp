//
// Created by ramizouari on 01/12/2021.
//
#include "linear_algebra.h"
#include <iostream>
#include "abstract_algebra.h"

int main()
{
    s_matrix<real,3,3> M({{0,1,1},{2,0,0},{0,0,0}});
    s_vector<real,3> v({21,0,7});
    std::cout << M.rank();
}