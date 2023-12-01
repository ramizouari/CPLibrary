//
// Created by ramizouari on 29/11/23.
//

#ifndef CPLIBRARY_TOWER_H
#define CPLIBRARY_TOWER_H

#include <memory>
#include <optional>
#include "algebra/abstract_algebra.h"
namespace cp
{
    template<typename R>
    struct quadratic_tower_node
    {
        std::vector<R> tower;
        std::vector<std::array<std::vector<R>,2>> modulos;
        quadratic_tower_node multiply(quadratic_tower_node &A, quadratic_tower_node &B, int a,int b, int height)
        {
            quadratic_tower_node result;
            result.tower.resize(b-a);
            if(b-a==1)
                return A.tower[a]*B.tower[a];
            result.modulos.reserve(modulos.size()-height);
            for(int i=0;i<modulos.size()-height;i++)
                result.modulos[i]=modulos[i];
            auto U=multiply(A,B,result,a,(a+b)/2,height+1);
            for(int i=a;i<(a+b)/2;i++)
                result.tower[i]+=U.tower[i];
            auto V=multiply(A,B,result,(a+b)/2,b,height+1);
            for(int i=(a+b)/2;i<b;i++)
                result.tower[i]+=V.tower[i];
            auto P= multiply(A,B,result,a,(a+b)/2,height+1);
            auto Q= multiply(A,B,result,(a+b)/2,b,height+1);
        }
        quadratic_tower_node()= default;
        quadratic_tower_node& operator+=(const quadratic_tower_node &x)
        {
            for(int i=0;i<tower.size();i++)
                tower[i]+=x.tower[i];
            return *this;
        }

        quadratic_tower_node& operator-=(const quadratic_tower_node &x)
        {
            for(int i=0;i<tower.size();i++)
                tower[i]-=x.tower[i];
            return *this;
        }

        quadratic_tower_node operator*(const quadratic_tower_node &x)
        {
            quadratic_tower_node result;
            result.modulos=modulos;
            for(int i=0;i<tower.size();i++)
                tower[i]*=x.tower[i];
            return *this;
        }
    };

    template<typename R>
    struct quadratic_tower
    {
        std::shared_ptr<quadratic_tower_node<R>> root;
        std::vector<std::shared_ptr<quadratic_tower_node<R>[2]>> modulus_towers;
    };
}

#endif //CPLIBRARY_TOWER_H
