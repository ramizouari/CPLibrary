//
// Created by ramizouari on 09/12/22.
//

#ifndef CPLIBRARY_SPECIAL_POLYNOMIALS_H
#define CPLIBRARY_SPECIAL_POLYNOMIALS_H

#include "vector.h"
#include "matrix.h"
#include "polynomial/polynomial.h"
namespace cp::linalg
{
    template<typename R, std::size_t extent = dynamic_extent>
    polynomial<R> faddev_lerrier_characteristic_polynomial(const matrix<R, extent> & A)
    {
        int n=A.row_dim();
        std::vector<R> S(n + 1);
        S[n] = 1;
        matrix<R> C(n,n,size_tag);
        for (int i = n - 1; i >= 0; i--)
        {
            for (int j = 0; j < n; j++)
                C[j][j] += S[i + 1];
            C = A * C;
            S[i] = -C.tr() / R(n - i);
        }
        return S;
    }


    template<typename R, std::size_t extent>
    polynomial<R> interpolation_characteristic_polynomial(matrix<R,extent> M)
    {
        int n=M.row_dim();
        std::vector<R> X(n+1), Y(n+1);
        for (int i = 0; i <= n; i++)
        {
            X[i] = i;
            Y[i] = M.det();
            for (int j = 0; j < n; j++)
                M[j][j] = M[j][j] - 1;
        }
        return newton_interpolation(X, Y);
    }

    template<typename IK, std::size_t extent>
    bool annihilable(const matrix<IK, extent> &T, const vector<IK, extent> &u,int m)
    {
        matrix<IK> A(m+1,u.dim(),size_tag);
        A[0]=u.data();
        for(int i=1;i<=m;i++)
        {
            auto v=T * A[i - 1];
            std::copy(v.begin(),v.end(),A[i].begin());
        }
        return A.nullity()!=0;
    }

    template<typename IK, std::size_t extent>
    polynomial<IK> minimal_polynomial(const matrix<IK, extent> &T, const vector<IK, extent> &u)
    {
        auto n=u.dim();
        std::vector<int> D(n+1);
        for(int i=0;i<=n;i++)
            D[i]=i;
        auto d=*std::upper_bound(D.begin(),D.end(),0,[&T,&u](const auto &x,const auto &y){return annihilable(T,u,x) < annihilable(T,u,y);});
        std::vector<int> mapper(d+1);
        std::vector<polynomial<IK>> Z(d+1);
        std::vector<vector<IK>> U(d+1);
        U[0]=u;
        for(int i=1;i<=d;i++)
            U[i] = T * U[i - 1];
        for(int i=0;i<=d;i++)
        {
            mapper[i]=i;
            Z[i].p.resize(d+1);
            Z[i][i]=1;
        }
        int r=0,t=0;
        while(r<=d && t<n)
        {
            if(std::all_of(U[mapper[r]].begin(),U[mapper[r]].end(),[](const auto &x){return is_zero(x);}))
                break;
            int s=r;
            while(s<=d && is_zero(U[mapper[s]][t]))
                s++;
            if(s==d+1)
            {
                t++;
                continue;
            }
            std::swap(mapper[s],mapper[r]);
            Z[mapper[r]]/=U[mapper[r]][t];
            for(int j=t+1;j<n;j++)
                U[mapper[r]][j]/=U[mapper[r]][t];
            U[mapper[r]][t]=1;
            for(int i=r+1;i<=d;i++)
            {
                auto w=U[mapper[i]][t];
                for(int j=t;j<n;j++)
                    U[mapper[i]][j]-=w*U[mapper[r]][j];
                U[mapper[i]][t]=0;
                Z[mapper[i]]-=w*Z[mapper[r]];
            }
            r++;
            t++;
        }
        return Z[mapper[d]]/Z[mapper[d]][Z[mapper[d]].degree()];
    }

    template<typename IK, std::size_t extent>
    polynomial<IK> minimal_polynomial(const matrix<IK, extent> &T, const std::vector<vector<IK, extent>> &U)
    {
        polynomial<IK> mu=1;
        for(auto &u:U)
            mu=lcm(mu, minimal_polynomial(T,u));
        return mu;
    }

    template<typename IK, std::size_t extent>
    polynomial<IK> minimal_polynomial(const matrix<IK, extent> &T)
    {
        std::vector<vector<IK,extent>> E;
        int n=T.row_dim();
        for(int i=0;i<n;i++)
        {
            E.emplace_back(n,size_tag);
            E.back()[i]=1;
        }
        return minimal_polynomial(T,E);
    }
}



#endif //CPLIBRARY_SPECIAL_POLYNOMIALS_H
