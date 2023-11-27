//
// Created by ramizouari on 26/10/23.
//

#ifndef CPLIBRARY_FIXED_SEGMENT_TREE_H
#define CPLIBRARY_FIXED_SEGMENT_TREE_H
#include <vector>
namespace cp::data_structures::fixed
{
    template<typename O>
    struct segment_tree
    {
        using R=typename O::type;
        using type=R;
        std::vector<std::vector<R>> S;
        std::vector<R> A;
        int n,h;
        segment_tree(const std::vector<R> &_A):A(_A)
        {
            n=bit_ceil(A.size());
            A.resize(n,O::neutral);
            int m=n;
            h=0;
            while(m)
            {
                m/=2;
                h++;
            }
            S.resize(h);
            for(int i=0;i<h;i++)
                S[i].resize(1<<i);
            build();
        }

        void update(int i,R u)
        {
            A[i]=u;
            S[h-1][i]=u;
            int m=h-2;
            i/=2;
            while(m>=0)
            {
                S[m][i]=F(S[m+1][2*i],S[m+1][2*i+1]);
                m--;
                i/=2;
            }
        }

        R query(int l,int r)
        {
            return query(std::max(l,0),std::min(r,n),0,n,0);
        }
    private:
        inline static O F=O();
        void build()
        {
            for(int i=0;i<n;i++)
                S.back()[i]=A[i];
            for(int i=h-2;i>=0;i--) for(int k=0;k<(1<<i);k++)
                    S[i][k]=F(S[i+1][2*k],S[i+1][2*k+1]);
        }
        R query(int l,int r,int a,int b,int depth)
        {
            if(l>=r)
                return O::neutral;
            if(l==a && r==b)
                return S[depth][l>>(h-1-depth)];
            int mid=(a+b)/2;
            if(mid>r)
                return query(l,r,a,mid,depth+1);
            else if(mid<l)
                return query(l,r,mid,b,depth+1);
            else
                return F(query(l,mid,a,mid,depth+1),query(mid,r,mid,b,depth+1));
        }
    };
}

#endif //CPLIBRARY_SEGMENT_TREE_H
