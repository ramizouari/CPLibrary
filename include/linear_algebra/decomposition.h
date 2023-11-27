//
// Created by ramizouari on 08/12/22.
//

#ifndef CPLIBRARY_DECOMPOSITION_H
#define CPLIBRARY_DECOMPOSITION_H

#include <unordered_set>
#include <set>
#include "matrix.h"
#include "special_matrices.h"

namespace linalg
{
    template<typename IK>
    std::pair<LowerTriangularSquareMatrix<IK>,UpperTriangularSquareMatrix<IK>> choleskyDecomposition(const d_matrix<IK> &P)
    {
        int n=P.row_dim();
        LowerTriangularSquareMatrix<IK> L(0,n);
        for(int i=0;i<n;i++) for(int j=0;j<=i;j++)
            {
                L[i][j]=P[i][j];
                for(int k=0;k<j;k++)
                    L[i][j]-=L[i][k]*conj<IK,IK>(L[j][k]);
                if(j<i)
                {
                    if(L[j][j]==IK(0))
                        L[i][j]=0;
                    else L[i][j]/=conj<IK,IK>(L[j][j]);
                }
                else L[i][j]=std::sqrt(std::abs(L[i][j]));
            }
        return std::make_pair(L,L.H());
    }

    template<typename IK>
    struct LDL_Decomposition_t
    {
        LowerTriangularSquareMatrix<IK> L;
        std::vector<IK> D;
        UpperTriangularSquareMatrix<IK> U;
        LDL_Decomposition_t(int n):L(1,n),D(n),U(1,n){}
    };

    template<typename IK>
    LDL_Decomposition_t<IK> LDL_Decomposition(const d_matrix<IK> &P)
    {
        int n=P.row_dim();
        LDL_Decomposition_t<IK> LDL(n);
        auto &[L,D,U]=LDL;
        for(int i=0;i<n;i++) for(int j=0;j<=i;j++)
            {
                if(j<i)
                {
                    L[i][j]=P[i][j];
                    for(int k=0;k<j;k++)
                        L[i][j]-= D[k]*L[i][k]*conj<IK,IK>(L[j][k]);
                    if(D[j]==IK(0))
                        L[i][j]=0;
                    else L[i][j]/=conj<IK,IK>(D[j]);
                }
                else
                {
                    D[i]=P[i][i];
                    for(int k=0;k<i;k++)
                        D[i]-=D[k]*L[i][k]*conj<IK,IK>(L[i][k]);
                }
            }
        U=L.H();
        return LDL;
    }


    template<typename IK>
    d_vector<IK> unit_proj(d_vector<IK> u,d_vector<IK> v)
    {
        IK w=0;
        for(int i=0;i<u.dim();i++)
            w+=conj<IK,IK>(u[i])*v[i];
        return w*u;
    }

    template<typename IK>
    std::vector<d_vector<IK>> gram_schmidt(const std::vector<d_vector<IK>>& A,real eps=1e-6)
    {
        if(A.empty())
            return A;
        int n=A[0].dim();
        std::vector<d_vector<IK>> B;
        std::set<int> J;
        for (int i = 0; i < A.size(); i++)
        {
            d_vector<IK> R(v_shape{n}), u(v_shape{n});
            for (auto& s : B)
                R += unit_proj(s,A[i]);
            u = A[i] - R;
            real N =0;
            for(auto a:u)
                N+=std::abs(a*conj<IK,IK>(a));
            if (N > eps)
            {
                u /= std::sqrt(N);
                J.insert(i);
                B.push_back(u);
            }
            else
                B.push_back(d_vector<IK>(v_shape{n}));
        }
        return B;
    }

    template<typename IK>
    std::vector<d_vector<IK>> complete_gram_schmidt(const std::vector<d_vector<IK>>& A,real eps=1e-6)
    {
        if(A.empty())
            return A;
        int n=A[0].dim();
        std::vector<d_vector<IK>> B;
        std::set<int> J;
        for (int i = 0; i < A.size(); i++)
        {
            d_vector<IK> R(v_shape{n}), u(v_shape{n});
            for (auto& s : B)
                R += unit_proj(s,A[i]);
            u = A[i] - R;
            real N =0;
            for(auto a:u)
                N+=std::abs(a*conj<IK,IK>(a));
            if (N > eps)
            {
                u /= std::sqrt(N);
                J.insert(i);
                B.push_back(u);
            }
            else B.emplace_back(v_shape{n});
        }
        int k=0;
        for(int i=0;i<B.size();i++) if(!J.contains(i))
            {
                real N=0;
                while(k<n && N<=eps)
                {
                    d_vector<IK> R(v_shape{n}), u(v_shape{n}),e(v_shape{n});
                    e[k]=1;
                    for (auto j:J)
                        R += unit_proj(B[j], e);
                    u = e - R;
                    N = 0;
                    for(auto a:u)
                        N+=std::abs(a*conj<IK,IK>(a));
                    if (N > eps)
                    {
                        u /= std::sqrt(N);
                        J.insert(i);
                        B[i]=u;
                    }
                    k++;
                }
            }
        return B;
    }

    template<typename IK>
    std::vector<d_vector<IK>> gram_schmidt(const d_matrix<IK>& A)
    {
        int n=A.row_dim(),m=A.col_dim();
        std::vector<d_vector<IK>> U(m,v_shape{n});
        for (int i = 0; i < n; i++) for(int j=0;j<m;j++)
                U[j][i]=A[j][i];
        return gram_schmidt(U);
    }


    template<typename IK>
    d_vector<IK> proj(const d_vector<IK> &u,const d_vector<IK> &v)
    {
        IK w=0,N=0;
        for(int i=0;i<u.dim();i++) {
            w += conj<IK, IK>(u[i]) * v[i];
            N+=conj<IK,IK>(u[i])*u[i];
        }
        return w/N*u;
    }

    template<typename IK>
    struct QR_Decomposition_t
    {
        d_matrix<IK> Q,R;
    };


/**
 * @issues Matrix [[0,1,1],[0,1,1],[0,1,1]] has a faulty QR Decomposition
 * */
    template<typename IK>
    QR_Decomposition_t<IK> QR_Decomposition(const d_matrix<IK>& A)
    {
        auto L=choleskyDecomposition(A.H() * A).second;
        auto R = d_matrix<IK>(L.M);
        std::vector<d_vector<IK>> G;
        auto Q = A*L.pinv().M;
        int n=A.row_dim(),m=A.col_dim();
        for (int i = 0; i < m; i++) {
            G.emplace_back(v_shape{n});
            for (int j = 0; j < n; j++)
                G.back()[j]=Q[j][i];
        }
        auto B= complete_gram_schmidt(G);
        for (int i = 0; i < n; i++)
            for(int j = 0; j < std::min(n,m); j++)
                Q[i][j] = B[j][i];
        return { Q,R };
    }

    template<typename IK>
    struct LQ_Decomposition_t
    {
        d_matrix<IK> L,Q;
    };


/**
 * @issues Matrix [[0,1,1],[0,1,1],[0,1,1]] has a faulty QR Decomposition
 * */
    template<typename IK>
    QR_Decomposition_t<IK> LQ_Decomposition(const d_matrix<IK>& A)
    {
        auto T=choleskyDecomposition(A*A.H()).first;
        auto L = d_matrix<IK>(T.M);
        std::vector<d_vector<IK>> G;
        auto Q = d_matrix<IK>(T.pinv().M)*A;
        int n=A.row_dim(),m=A.col_dim();
        for (int i = 0; i < n; i++) {
            G.emplace_back(v_shape{m});
            for (int j = 0; j < m; j++)
                G.back()[j]=Q[i][j];
        }


        auto B= gram_schmidt(G);
        for (int i = 0; i < m; i++)
            for(int j = 0; j < m; j++)
                Q[i][j] = B[i][j];
        return { L,Q };
    }

    template<typename IK>
    struct SchurDecomposition_t
    {
        d_matrix<IK> P,T;
    };

    template<typename IK>
    SchurDecomposition_t<IK> QR_Algorithm(d_matrix<IK> A,int iter)
    {
        d_matrix<IK> U(1,m_shape{(int)A.col_dim(),(int)A.row_dim()});
        while(iter--)
        {
            auto [Q,R]=QR_Decomposition(A);
            A=R*Q;
            U*=Q;

        }
        return {U,A};
    }

    template<typename IK>
    struct SVD_t
    {
        d_matrix<IK> U,D,V;
        SVD_t(int n,int m):U(1,m_shape{n,n}),D(1,m_shape{n,m}),V(1,m_shape{m,m}){}
    };

    template<typename IK>
    SVD_t<IK> SVD(d_matrix<IK> A,int iter)
    {
        int n=A.row_dim(),m=A.col_dim();
        SVD_t<IK> UDV(n,m);
        while (iter--)
        {
            auto [Q, R] = QR_Decomposition(A);
            auto [L, P] = LQ_Decomposition(R);
            A = L;
            UDV.U = UDV.U * Q;
            UDV.V = P * UDV.V;
        }
        UDV.D = A;
        return UDV;
    }

    template<typename IK>
    struct PolarDecomposition_t
    {
        d_matrix<IK> U,P;
    };

    template<typename IK>
    PolarDecomposition_t<IK> polarDecomposition(const d_matrix<IK>&A,int iter)
    {
        auto [U, D, V] = SVD(A,iter);
        return PolarDecomposition_t{ U,D * V };
    }

    template<typename IK>
    d_matrix<IK> pinv(const d_matrix<IK>&A,int iter,real eps=1e-7)
    {
        auto [U,D,V]= SVD(A,iter);
        for(int i=0;i<std::min(A.col_dim(),A.row_dim());i++)
            if(std::abs(D[i][i]) > eps)
                D[i][i]=IK(1)/D[i][i];
        return V.H()*D*U.H();
    }

    template<typename IK>
    d_matrix<IK> sqrt(const d_matrix<IK> &A,int iter)
    {
        auto [Q,D]= QR_Algorithm(A,iter);
        for(int i=0;i<std::min(D.col_dim(),D.row_dim());i++)
            D[i][i]=std::sqrt(D[i][i]);
        return Q*D*Q.H();
    }

    template<typename IK,typename Function>
    d_matrix<IK> symmetric_matrix_function(const d_matrix<IK> &A,const Function &F,int iter)
    {
        auto [Q,D]= QR_Algorithm(A,iter);
        for(int i=0;i<std::min(D.col_dim(),D.row_dim());i++)
            D[i][i]=F(D[i][i]);
        return Q*D*Q.H();
    }
}

#endif //CPLIBRARY_DECOMPOSITION_H
