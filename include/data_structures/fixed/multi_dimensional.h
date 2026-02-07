//
// Created by ramizouari on 06/01/2026.
//

#ifndef CPLIBRARY_NESTED_H
#define CPLIBRARY_NESTED_H
#include <algorithm>
#include <bitset>
#include <queue>
#include <stack>
#include <vector>

#include "algebra/bits.h"
#include "functional/zip.h"
#include "../../tensors/tensor.h"
#include "../../tensors/view.h"

namespace cp::data_structures::fixed {

    using index_vector=std::vector<std::size_t>;
    using index_span = std::span<std::size_t>;

    template<sized_random_access Q>
    Q bit_ceil_array(const Q& N)
    {
        Q res(N);
        for(auto &n:res) n=bit_ceil(n);
        return res;
    }

    template<sized_random_access Q>
    Q bit_log_array(const Q& N)
    {
        Q res(N);
        for(auto &n:res) n=bit_log(n);
        return res;
    }

    template<sized_random_access Q>
    Q build_depth(const Q& N)
    {
        Q H;
        if constexpr (std::is_same_v<Q,index_vector>)
            H.resize(N.size());
        for (int i=0;i<N.size();i++) {
            H[i]=0;
            std::size_t m=N[i];
            while(m)
            {
                m/=2;
                ++H[i];
            }
        }
        return H;
    }

    template <typename R>
    linalg::flat_tensor<R,dynamic_extent> resize_tensor(const linalg::tensor_view<R, dynamic_extent> &A,
        index_vector N, const R& value)
    {
        auto M = A.shape();
        for (int i=0;i<M.size();i++)
            M[i] = std::min(M[i],N[i]);
        linalg::flat_tensor<R,dynamic_extent> B(N, value);
        for (const auto &I : iterate_multi_index_right_up(M))
            B.at(I)=A.at(I);
        return B;
    }

    template <typename R, std::size_t extent>
    linalg::flat_tensor<R,extent> resize_tensor(const linalg::tensor_view<R, extent> &A,
        std::array<std::size_t,extent> N, const R& value)
    {
        auto M = A.shape();
        for (int i=0;i<M.size();i++)
            M[i] = std::min(M[i],N[i]);
        linalg::flat_tensor<R,extent> B(N, value);
        for (const auto &I : iterate_multi_index_right_up(M))
            B.at(I)=A.at(I);
        return B;
    }

    template<typename O, std::size_t extent = dynamic_extent>
    struct multidimensional_segment_tree
    {
        using R= O::type;
        using index_array = linalg::flat_tensor<R,extent>::index_array;
        using type = R;
        using value_type = type;
        using key_type = index_array;
        index_array N,H;
        linalg::flat_tensor<R,extent> A;
        linalg::flat_tensor<linalg::flat_tensor<R,extent>,extent> S;
        multidimensional_segment_tree(const linalg::tensor_view<R, extent> &A):N(bit_ceil_array(A.shape())),
            H(build_depth(N)),A(resize_tensor(A,N,O::neutral)), S(H,linalg::flat_tensor<R,extent>(N))
        {
            for (const auto &I: iterate_multi_index_right_up(H))
            {
                auto J=I;
                for (auto &j:J)
                    j = 1<<j;
                S.at(I) = linalg::flat_tensor<R,extent>(J);
            }
            build();
        }

        void update(index_array I,R u) {
            auto Z = H;
            for (auto &z:Z) --z;
            using quadruplet = std::tuple<index_array,index_array,std::size_t,bool>;
            std::queue<quadruplet> Q;
            Q.emplace(Z,I,0,true);
            while (!Q.empty())
            {
                auto [U,V,s,initial] = Q.front();
                Q.pop();
                if (initial) S.at(U).at(V) = A.at(V) = u;
                else
                {
                    S.at(U).at(V) = O::neutral;
                    index_vector K;
                    for (auto [e,h] : H) if (e < h-1) K.push_back(e);
                    auto X = I;
                    for (auto k:K) ++X[k];
                    for (std::size_t b=0;b<1<<K.size();b++)
                    {
                        auto Y = V;
                        std::bitset<32> B(b);
                        for (std::size_t t=0;t<K.size();t++) Y[K[t]] = 2*V[K[t]]+B[t];
                        auto & x = S.at(U).at(V);
                        x = F(S.at(X).at(), x);
                    }
                }
                for (auto t=s; t<N.size();++t) if (U[t])
                {
                    --U[t];
                    auto z = V[t]%2;
                    V[t]/=2;
                    Q.emplace(U,V,t,false);
                    ++U[t];
                    V[t]=2*V[t]+z;
                }
            }
        }

        R query(index_array L, index_array R)
        {
            for (auto &l:L) l = std::max<std::size_t>(l,0);
            for (int i=0;i<N.size();i++) L[i]=std::min(L[i],N[i]-1);
            if constexpr (extent == dynamic_extent) {
                index_array Z(H.size());
                return query(L,R,Z,N,Z);
            }
            else {
                index_array Z{};
                return query(L,R,Z,N,Z);
            }

        }
    private:
        inline static O F=O();
        void build()
        {

            for (auto I : iterate_multi_index_right_down(H)) {
                index_array admissible=N;
                std::size_t admissible_size=0;
                for (int j=0;j<I.size();j++) if (I[j] < H[j]-1)
                    admissible[admissible_size++]=j;
                if (!admissible_size) for (const auto &J: iterate_multi_index_right_up(N))
                    S.at(I).at(J) = A.at(J);
                else
                {
                    auto L = I;
                    for (auto &l:L) l = 1 << l;
                    auto P=I;
                    for (int t=0;t<admissible_size;t++)
                        ++P[admissible[t]];
                    for (const auto &J: iterate_multi_index_right_up(L))
                    {
                        auto Q=J;
                        for (int t=0;t < admissible_size;t++) Q[admissible[t]]*=2;
                        for (size_t z = 0; z< 1<<admissible_size;z++)
                        {
                            auto Z = std::bitset<32>(z);
                            for (size_t k=0;k<admissible_size;k++) if (Z[k]) Q[admissible[k]]+=1;
                            auto &x=S.at(I).at(J);
                            x  = F(x,S.at(P).at(Q));
                            for (size_t k=0;k<admissible_size;k++) if (Z[k]) Q[admissible[k]]-=1;
                        }

                    }
                }
            }
        }
        R query(index_array L,index_array R,index_array A,index_array B,index_array D, std::size_t dim=0)
        {
            if (L[dim] >= R[dim]) return O::neutral;
            while (dim < H.size() && L[dim] == A[dim] && R[dim] == B[dim])
                ++dim;
            if(dim == H.size()) {
                if constexpr (extent == dynamic_extent) {
                    index_vector K(L.size());
                    for (auto [k,l,h,depth]: zip(K,L,H,D)) k = l>>(h-1-depth);
                    return S.at(D).at(K);
                }
                else {
                    index_array K;
                    for (auto [k,l,h,depth]: zip(K,L,H,D)) k = l>>(h-1-depth);
                    return S.at(D).at(K);
                }

            }
            auto mid=(A[dim]+B[dim])/2;
            if(mid>R[dim]) {
                auto M = B;
                M[dim] = mid;
                ++D[dim];
                return query(L,R,A,M,D,dim);
            }
            if(mid<L[dim]) {
                auto M = A;
                M[dim] = mid;
                ++D[dim];
                return query(L,R,M,B,D,dim);
            }
            auto L1 = L, R2= R, A1=A,B2=B;
            L1[dim] = mid;
            R2[dim] = mid;
            ++D[dim];
            B2[dim] = mid;
            A1[dim] = mid;
            return F(query(L,R2,A,B2,D,dim),query(L1,R,A1,B,D,dim));
        }
    };


    template<typename O, std::size_t extent = dynamic_extent>
    struct multidimensional_sparse_array
    {
        using T=O::type;
        using R= O::type;
        using index_array = linalg::flat_tensor<R,extent>::index_array;
        using type = T;
        using value_type = type;
        using key_type = index_array;
        inline static O F=O();
        index_array N,H;
        linalg::flat_tensor<linalg::flat_tensor<R,extent>,extent> S;

        index_array shape_add(index_array A,std::size_t offset) {
            for (auto &a:A) a+=offset;
            return A;
        }

        multidimensional_sparse_array(const linalg::tensor_view<R, extent> &A):N(bit_ceil_array(A.shape())),
            H(bit_log_array(N)), S(shape_add(H,1),linalg::flat_tensor<R,extent>(N))
        {
            for (auto I : iterate_multi_index_right_down(H))
            {
                index_array K=N;
                for (auto [k,n,h,i]:zip(K,N,H,I))
                    k -= (1<<(h-i))-1;
                S.at(I) = linalg::flat_tensor<R,extent>(K);
            }
            for (auto I : iterate_multi_index_right_down(shape_add(H,1)))
            {
                index_array admissible=N;
                std::size_t admissible_size=0;

                for (int j=0;j<I.size();j++) if (I[j] < H[j]) admissible[admissible_size++]=j;
                if (!admissible_size) for (const auto &J: iterate_multi_index_right_up(A.shape()))
                    S.at(I).at(J) = A.at(J);
                else
                {
                    auto L = N,R=H;
                    for (int t =0; t< admissible_size;t++) {
                        auto k=admissible[t];
                        R[k]-=I[k]+1;
                        R[k] = 1<<R[k];
                        L[k] -= 2*R[k];
                        ++L[k];
                    }
                    for (const auto &J: iterate_multi_index_right_up(L))
                    {
                        auto Q=J;
                        auto P=I;
                        for (int t=0; t < admissible_size; t++)
                            ++P[admissible[t]];
                        for (std::size_t b=0;b< 1<< admissible_size;b++)
                        {
                            std::bitset<32> B(b);
                            for (std::size_t t=0;t<admissible_size;t++) if (B[t])
                                Q[admissible[t]]+=R[admissible[t]];
                            auto &x=S.at(I).at(J);
                            x = F(x,S.at(P).at(Q));
                            for (std::size_t t=0;t<admissible_size;t++) if (B[t])
                                Q[admissible[t]]-=R[admissible[t]];
                        }
                    }
                }
            }
        }

        T query(index_array L,index_array R) const
        {
            T res= O::neutral;
            for (auto [l,r]:zip(L,R)) if(l>=r)
                return O::neutral;
            auto D=R;
            for (auto [d,l]:zip(D,L)) d-=l;
            auto W=D;
            for (auto&w:W) w = bit_floor(w);
            auto B=W;
            for (auto&b:B) b = bit_log(b);
            auto Z=H;
            for (auto [z,b]:zip(Z,B)) z-=b;
            auto T = R;
            for (auto [t,w] : zip(T,W)) t-=w;
            for (std::size_t i=0;i< 1<< N.size();i++) {
                std::bitset<32> I(i);
                for (auto k=0;k<N.size();k++) T[k] = I[k]?R[k] - W[k]:L[k];
                res = F(res,S.at(Z).at(T));
            }
            return res;
        }
    };

    template<typename O, std::size_t extent = dynamic_extent>
    struct multidimensional_fenwick_tree {
        using index_array=linalg::flat_tensor<std::size_t,extent>::index_array;
        using type=O::type;
        using value_type = type;
        using key_type = index_array;
        index_array N;
        using T=O::type;
        using R=O::type;
        linalg::flat_tensor<R,extent> bit;
        inline static O F = O();

        multidimensional_fenwick_tree(index_array N):N(std::move(N)),bit(this->N,O::neutral){}
        multidimensional_fenwick_tree(const linalg::flat_tensor<R,extent> &X) : multidimensional_fenwick_tree(X.shape())
        {
            for (auto I : iterate_multi_index_right_up(N))
                update(I,X.at(I));
        }
        T sum(index_array I)
        {
            for (auto i:I) if (!~i) return O::neutral;
            T ret = O::neutral;
            std::stack<std::pair<index_array,std::size_t>> Q;
            Q.emplace(I,0);
            while (!Q.empty())
            {
                auto [J,s] = Q.top();
                ret = F(ret,bit.at(J));
                Q.pop();
                for (auto i = s; i < N.size();i++)
                {
                    auto K=J;
                    K[i] = (K[i] & K[i]+1) - 1;
                    if (~K[i]) Q.emplace(K,i); // K[i] == -1, but K[i] unsigned
                }
            }
            return ret;
        }

        T query(index_array A,index_array B)
        {
            T ret = O::neutral;
            for (auto &a:A) --a;
            for (size_t i=0; i < 1<<N.size();i++)
            {
                const auto s = N.size()-std::popcount(i);
                std::bitset<32> I(i);
                auto C=A;
                bool skip=false;
                for (auto k=0;k<N.size();k++)
                {
                    if (I[k]) C[k] = B[k];
                    else if (~A[k] == 0) {
                        skip=true;
                        break;
                    }
                }
                if (skip) continue;
                auto delta = sum(C);
                if (s%2 == 0) ret=F(ret,delta);
                else ret=F(ret,F.inv(delta));
            }
            return ret;
        }

        T sum(index_array A,index_array B)
        {
            return query(A,B);
        }

        void add(index_array I, T delta) {
            T ret = O::neutral;
            std::queue<std::pair<index_array,std::size_t>> Q;
            Q.emplace(I,0);
            while (!Q.empty())
            {
                auto [J,s] = Q.front();
                bit.at(J) = F(bit.at(J),delta);
                Q.pop();
                for (auto i = s; i < N.size();i++)
                {
                    auto K=J;
                    K[i] = K[i] | K[i]+1;
                    if (K[i] < N[i]) Q.emplace(K,i);
                }
            }
        }

        void update(index_array I, T delta) {
            add(I,F(F.inv(sum(I,I)),delta));
        }
    };
}

#endif //CPLIBRARY_NESTED_H