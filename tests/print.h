//
// Created by ramizouari on 11/04/2022.
//

#ifndef CPLIBRARY_PRINT_H
#define CPLIBRARY_PRINT_H
#include <ostream>
#include "../include/linalg/matrix.h"

template<typename T,int n>
std::ostream & operator<<(std::ostream & os,const cp::linalg::vector<T,n> & v)
{
    os << '[';
    for(auto &s:v)
        os<<s<<",";
    os  << "]";
    return os;
}

template<typename T,int n,int m>
std::ostream & operator<<(std::ostream & os,const cp::linalg::matrix<T,n,m> & v)
{
    os << '[';
    for(auto &R:v) {
        os << '[';
        for (auto &s:R)
            os << s << ",";
        os << "],";
    }
    os  << "]";
    return os;
}




#endif //CPLIBRARY_PRINT_H
