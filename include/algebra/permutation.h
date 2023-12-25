//
// Created by ramizouari on 22/12/23.
//

#ifndef CPLIBRARY_PERMUTATION_H
#define CPLIBRARY_PERMUTATION_H
#include "abstract_algebra.h"
#include <numeric>
namespace cp
{

    struct abstract_permutation
    {
        virtual int transform(int i) const = 0;
        int operator()(int i) const
        {
            return transform(i);
        }
        virtual int size() const = 0;
        virtual ~abstract_permutation() = default;

        std::vector<std::vector<int>> cycles() const
        {
            std::vector<std::vector<int>> R;
            std::vector<bool> visited(size());
            for(int i=0;i<size();i++)
            {
                if(!visited[i])
                {
                    R.push_back({});
                    int j=i;
                    do
                    {
                        visited[j]=true;
                        R.back().push_back(j);
                        j=transform(j);
                    }while(j!=i);
                }
            }
            return R;
        }

        int number_of_cycles() const
        {
            std::vector<bool> visited(size());
            int R=0;
            for(int i=0;i<size();i++)
            {
                if(!visited[i])
                {
                    R++;
                    int j=i;
                    do
                    {
                        visited[j]=true;
                        j=transform(j);
                    }while(j!=i);
                }
            }
            return R;
        }
    };

    struct permutation : public abstract_permutation
    {
        std::vector<int> P;
        permutation(int n):P(n)
        {
            std::iota(P.begin(),P.end(),0);
        }
        permutation() = default;
        permutation(const std::vector<int> &_P):P(_P){}

        int transform(int i) const override
        {
            return i<P.size()?P[i]:i;
        }

        int size() const override
        {
            return P.size();
        }

        permutation operator*(const permutation &Q) const
        {
            permutation R(P.size());
            for(int i=0;i<P.size();i++)
                R.P[i]=P[Q(i)];
            return R;
        }

        permutation& operator*= (const permutation &Q)
        {
            return *this=*this*Q;
        }

        permutation inv() const
        {
            permutation Q(P.size());
            for(int i=0;i<P.size();i++)
                Q.P[P[i]]=i;
            return Q;
        }

        const int& operator[](int i) const
        {
            return P[i];
        }

        int& operator[](int i)
        {
            return P[i];
        }

    };
}

#endif //CPLIBRARY_PERMUTATION_H
