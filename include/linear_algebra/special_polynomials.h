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
    template<typename R>
    polynomial<R> faddev_lerrier_characteristic_polynomial(const d_matrix<R>&A)
    {
        int n=A.row_dim();
        std::vector<R> S(n + 1);
        S[n] = 1;
        d_matrix<R> C(0,m_shape{n,n});
        for (int i = n - 1; i >= 0; i--)
        {
            for (int j = 0; j < n; j++)
                C[j][j] += S[i + 1];
            C = A * C;
            S[i] = -C.tr() / R(n - i);
        }
        return S;
    }

    template<typename R,int n>
    polynomial<R> faddev_lerrier_characteristic_polynomial(const s_matrix<R,n,n>&A)
    {
        std::vector<R> S(n + 1);
        S[n] = 1;
        s_matrix<R,n,n> C;
        for (int i = n - 1; i >= 0; i--)
        {
            for (int j = 0; j < n; j++)
                C[j][j] += S[i + 1];
            C = A * C;
            S[i] = -C.tr() / R(n - i);
        }
        return S;
    }

    template<typename R>
    polynomial<R> interpolation_characteristic_polynomial(d_matrix<R> M)
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

    template<typename R,int n>
    polynomial<R> interpolation_characteristic_polynomial(s_matrix<R,n,n> M)
    {
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

    template<typename IK>
    bool annihilable(const d_matrix<IK> &T, const d_vector<IK> &u,int m)
    {
        d_matrix<IK> A(0,m_shape{m+1,(int)u.dim()});
        A[0]=(std::vector<IK>&)u;
        for(int i=1;i<=m;i++)
        {
            auto v=T * A[i - 1];
            std::copy(v.begin(),v.end(),A[i].begin());
        }
        return A.nullity()!=0;
    }

    template<typename IK>
    polynomial<IK> minimal_polynomial(const d_matrix<IK> &T, const d_vector<IK> &u)
    {
        int n=u.dim();
        std::vector<int> D(n+1);
        for(int i=0;i<=n;i++)
            D[i]=i;
        auto d=*std::upper_bound(D.begin(),D.end(),0,[&T,&u](const auto &x,const auto &y){return annihilable(T,u,x) < annihilable(T,u,y);});
        std::vector<int> mapper(d+1);
        std::vector<polynomial<IK>> Z(d+1);
        std::vector<d_vector<IK>> U(d+1);
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

    template<typename IK>
    polynomial<IK> minimal_polynomial(const d_matrix<IK> &T, const std::vector<d_vector<IK>> &U)
    {
        polynomial<IK> mu=1;
        for(auto &u:U)
            mu=lcm(mu, minimal_polynomial(T,u));
        return mu;
    }

    template<typename IK>
    polynomial<IK> minimal_polynomial(const d_matrix<IK> &T)
    {
        std::vector<d_vector<IK>> E;
        int n=T.row_dim();
        for(int i=0;i<n;i++)
        {
            E.emplace_back(v_shape{n});
            E.back()[i]=1;
        }
        return minimal_polynomial(T,E);
    }
}



#endif //CPLIBRARY_SPECIAL_POLYNOMIALS_H
