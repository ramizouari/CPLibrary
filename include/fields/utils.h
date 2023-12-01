//
// Created by ramizouari on 29/11/23.
//

#ifndef CPLIBRARY_UTILS_H
#define CPLIBRARY_UTILS_H
#include <cstdint>

namespace cp
{
    template<typename T>
    T xor_multiply(T a,T b)
    {
        T z=0;
        for(int i=0;i<sizeof(T)*8;i++)
        {
            z^=(a&1)*(b&1);
            a>>=1;
            b>>=1;
        }
        return z;
    }

}
#endif //CPLIBRARY_UTILS_H
