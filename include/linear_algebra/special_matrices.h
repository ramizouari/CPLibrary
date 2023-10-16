//
// Created by ramizouari on 08/12/22.
//

#ifndef CPLIBRARY_SPECIAL_MATRICES_H
#define CPLIBRARY_SPECIAL_MATRICES_H
#include "vector"
#include "algebra/abstract_algebra.h"

template<typename R>
struct BandMatrix
{
    std::vector<std::vector<R>> M;
    int n;
    int r;
public:
    BandMatrix(std::vector<std::vector<R>> _M):M(std::move(_M)),n(M.size()),r(M[0].size()/2){}
    R operator()(int i,int j) const
    {
        if(i>=n || i<0)
            return R{};
        else if(j-i+r <0 || j-i+r >=M[i].size())
            return R{};
        else return at(i,j);
    }
    auto &operator[](int i)
    {
        return M[i];
    }

    R& at(int i,int j)
    {
        return M[i][j-i+r];
    }

    const R& at(int i,int j) const
    {
        return M[i][j-i+r];
    }

    [[nodiscard]] std::vector<R> solve(std::vector<R> u) const
    {
        auto C=*this;
        for(int i=0;i<n;i++) for(int p=1;p<=r && i+p<n;p++)
            {
                auto w=C(i+p,i)/C(i,i);
                for(int q=i;q<=i+p+r;q++)
                    C.at(i+p, q) -= w * C(i, q);
                u[i+p]-=w*u[i];
            }

        for(int i=n-1;i>=0;i--)
        {
            auto w=C(i,i);
            u[i]/=w;
            for(int s=1;s<=r && i-s >=0;s++)
                u[i - s] -=   C(i - s,i)*u[i];
        }
        return u;
    }

    R det() const
    {
        auto C=*this;
        for(int i=0;i<n;i++) for(int p=1;p<=r && i+p<n;p++)
            {
                auto w=C(i+p,i)/C(i,i);
                for(int q=i;q<=i+p+r;q++)
                    C.at(i+p, q) -= w * C(i, q);
            }
        R w=1;
        for(int i=0;i<n;i++)
            w*=M(i,i);
        return w;
    }
};

template<typename R,bool isUpper>
struct TriangularSquareMatrix
{
    std::vector<std::vector<R>> M;
    int n;
public:
    TriangularSquareMatrix(R k,int n):n(n)
    {
        M.resize(n);
        for(int i=0;i<n;i++)
        {
            M[i].resize(n);
            M[i][i]=k;
        }
    }
    TriangularSquareMatrix(std::vector<std::vector<R>> _M):M(std::move(_M)),n(M.size()){}
    R operator()(int i,int j) const
    {
        if(i>=n || i<0)
            return R{};
        else return at(i,j);
    }
    auto &operator[](int i)
    {
        return M[i];
    }

    R& at(int i,int j)
    {
        return M[i][j];
    }

    const R& at(int i,int j) const
    {
        return M[i][j];
    }

    TriangularSquareMatrix<R,!isUpper> T() const
    {
        std::vector<std::vector<R>> A(n,std::vector<R>(n));
        for(int i=0;i<n;i++) for(int j=0;j<n;j++)
            A[i][j]=M[j][i];
        return TriangularSquareMatrix<R,!isUpper>(std::move(A));
    }

    TriangularSquareMatrix<R,!isUpper> H() const
    {
        std::vector<std::vector<R>> A(n,std::vector<R>(n));
        for(int i=0;i<n;i++) for(int j=0;j<n;j++)
                A[i][j]=conj<R,R>(M[j][i]);
        return TriangularSquareMatrix<R,!isUpper>(std::move(A));
    }

    TriangularSquareMatrix inv() const
    {
        TriangularSquareMatrix I(1,n);
        if constexpr (isUpper) for(int i=n-1;i>=0;i--)
        {
            for(int k=i;k<n;k++)
                I.M[i][k]/=M[i][i];
            for(int k=i;k<n;k++)
                for(int j=i-1;j>=0;j--)
                    I.M[j][k]-=M[j][i]*I.M[i][k];
        }
        else for(int i=0;i<n;i++)
        {
            for(int k=0;k<=i;k++)
                I.M[i][k]/=M[i][i];
            for(int k=0;k<=i;k++) for(int j=i+1;j<n;j++)
                    I.M[j][k]-=M[j][i]*I.M[i][k];
        }
        return I;
    }

    TriangularSquareMatrix pinv() const
    {
        TriangularSquareMatrix I(1,n);
        if constexpr (isUpper) for(int i=n-1;i>=0;i--)
        {
            if(M[i][i]==R{})
                continue;
            for(int k=i;k<n;k++)
                I.M[i][k]/=M[i][i];
            for(int k=i;k<n;k++)
                for(int j=i-1;j>=0;j--)
                I.M[j][k]-=M[j][i]*I.M[i][k];
        }
        else for(int i=0;i< n;i++)
        {
            if(M[i][i]==R{})
                continue;
            for(int k=0;k<=i;k++)
                I.M[i][k]/=M[i][i];

            for(int k=0;k<=i;k++) for(int j=i+1;j<n;j++)
                I.M[j][k]-=M[j][i]*I.M[i][k];
        }
        return I;
    }

    [[nodiscard]] std::vector<R> solve(std::vector<R> u) const
    {
        if constexpr (isUpper) for(int i=n-1;i>=0;i--)
        {
            u[i]/=M[i][i];
            for(int j=i-1;j>=0;j--)
                u[j]-=M[j][i]*u[i];
        }
        else for(int i=0;i<n;i++)
        {
            u[i]/=M[i][i];
            for(int j=i+1;j<n;j++)
                u[j]-=M[j][i]*u[i];
        }
        return u;
    }

    R det() const
    {
        R w=1;
        for(int i=0;i<n;i++)
            w*=M[i][i];
        return w;
    }
};

template<typename R>
using UpperTriangularSquareMatrix=TriangularSquareMatrix<R,true>;
template<typename R>
using LowerTriangularSquareMatrix=TriangularSquareMatrix<R,false>;
#endif //CPLIBRARY_SPECIAL_MATRICES_H
