//
// Created by ramizouari on 26/10/23.
//

#ifndef CPLIBRARY_BITS_H
#define CPLIBRARY_BITS_H
namespace cp
{
    inline unsigned int bit_log(unsigned int n)
    {
        unsigned char a=0,b=30,r=0;
        while(a<=b)
        {
            auto c=(a+b)/2;
            if(n>>c)
                a=c+1;
            else
            {
                b=c-1;
                r=c-1;
            }
        }
        if(r && (1<<(r-1))==n)
            return r-1;
        return r;
    }

    inline unsigned int bit_floor(unsigned int n)
    {
        return 1<<bit_log(n);
    }

    inline unsigned int bit_ceil(unsigned int n)
    {
        unsigned r=1;
        while(r<n)
            r<<=1;
        return r;
    }
}

#endif //CPLIBRARY_BITS_H
