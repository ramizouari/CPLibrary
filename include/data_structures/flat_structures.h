//
// Created by ramizouari on 22/11/22.
//

#ifndef CPLIBRARY_FLAT_STRUCTURES_H
#define CPLIBRARY_FLAT_STRUCTURES_H

#include <array>
#include <vector>
#include <cstdint>
#include <numeric>

using integer=std::int64_t;
/*
 * This file contains the following data structures:
 * 1) static_vector<T,capacity>: a vector with a fixed capacity
 * 2) static_flat_map<K,V,M>: a map with a fixed capacity
 * 3) flat_map<K,V,M>: a map with a dynamic capacity (implemented as a vector)
 *
 * */

namespace data_structures
{
    template<typename T,std::size_t capacity>
    struct static_vector:public std::array<T,capacity>
    {
        template<typename B=T>
        struct iterator
        {
            B* ptr;
            int size;
            iterator(B*_ptr,int _size):ptr(_ptr),size(_size){}
            iterator& operator++()
            {
                ptr++;
                return *this;
            }
            iterator operator++(int)
            {
                return iterator{ptr+1,size};
            }
            iterator& operator--()
            {
                ptr--;
                return *this;
            }
            iterator operator--(int)
            {
                return iterator{ptr-1,size};
            }
            iterator& operator+=(int r)
            {
                ptr+=r;
                return *this;
            }
            iterator operator+(int r)
            {
                return iterator{ptr+r,size};
            }
            iterator& operator-=(int r)
            {
                ptr-=r;
                return *this;
            }
            iterator operator-(int r)
            {
                return iterator{ptr-r,size};
            }
            auto operator*()
            {
                return *ptr;
            }
            bool operator==(const iterator&) const=default;
        };
        using std::array<T,capacity>::array;
        using std::array<T,capacity>::operator[];
        using std::array<T,capacity>::empty;
        int _size=0;
        std::size_t size()
        {
            return _size;
        }
        void push_back(const T&X)
        {
            if(_size==capacity)
                return;
            (*this)[_size++]=X;
        }
        void pop_back()
        {
            _size--;
        }

        auto begin()
        {
            return iterator<T>(this->data(),_size);
        }
        auto end()
        {
            return iterator<T>(this->data()+_size,_size);
        }

        auto begin() const
        {
            return iterator<const T>(this->data(),_size);
        }
        auto end() const
        {
            return iterator<const T>(this->data()+_size,_size);
        }

        bool operator==(const static_vector &Y) const
        {
            return std::equal(begin(),end(),Y.begin());
        }
        std::partial_ordering operator<=>(const static_vector &Y) const
        {
            return std::lexicographical_compare_three_way(begin(),end(),Y.begin(),Y.end());
        }

    };

    template<typename V,typename Less=std::less<>,int M=3>
    struct static_flat_map:public std::array<std::pair<integer,V>,M>
    {
        using size_t=std::size_t;
        using K=integer;
        using couple=std::pair<K,V>;
    public:

        using std::array<couple,M>::array;
        bool contains(const K&a)
        {
            for(auto &[b,_]:*this)
                if(a==b)
                    return true;
            return false;
        }

        size_t count(const K&a)
        {
            size_t R=0;
            for(auto &[b,_]:*this)
                if(a==b)
                    R++;
            return R;
        }

        void erase(const K&a)
        {
            int s=0;
            while(s<M && std::array<couple,M>::operator[](s).first!=a)
                s++;
            for(int i=s;i<M;i++)
                (*this)[i]=(*this)[i+1];
        }

        auto erase(typename std::array<couple,M>::iterator it)
        {
            auto Z=this->begin();
            while(Z!=this->end() && Z->first!=it->first)
                ++Z;
            for(int i=std::distance(this->begin(),Z);i<M;i++)
                (*this)[i]=(*this)[i+1];
        }


        auto find(const K &x)
        {
            auto it=this->begin();
            while(it!=this->end())
                if(it->first==x)
                    return it;
            return it;
        }

        auto emplace(const K&k,const V&v)
        {
            auto it=find(k);
            if(it!=this->end()) {
                it->second = v;
                return std::make_pair(it,false);
            }
            else for(int i=M-1;i>0;i++)
                {
                    if(i==0|| k> std::array<couple,M>::operator[](i-1).first)
                    {
                        std::array<couple,M>::operator[](i).first=k;
                        std::array<couple,M>::operator[](i).second=v;
                        return std::make_pair(this->begin()+i,true);
                    }
                    std::array<couple,M>::operator[](i)=std::array<couple,M>::operator[](i-1);
                }
        }

        V& operator[](size_t x)
        {
            auto it=find(x);
            if(it!=(*this).end())
                return it->second;
            else {
                auto [Z,_]=emplace(x, V{});
                return Z->second;
            }
        }

    };

    template<typename V,typename Less=std::less<>,int M=3>
    struct flat_map:public std::vector<std::pair<integer,V>>
    {
        using K=integer;
        using couple=std::pair<K,V>;
    public:

        template<typename ...T>
        flat_map(T&&... X):std::vector<couple>(std::forward(X)...)
        {
            this->reserve(M);
        }
        bool contains(const K&a)
        {
            for(auto &[b,_]:*this)
                if(a==b)
                    return true;
            return false;
        }

        std::size_t count(const K&a)
        {
            std::size_t R=0;
            for(auto &[b,_]:*this)
                if(a==b)
                    R++;
            return R;
        }

        void erase(const K&a)
        {
            int s=0;
            while(s<this->size() && std::vector<couple>::operator[](s).first!=a)
                s++;
            for(int i=s;i+1<this->size();i++)
                std::vector<couple>::operator[](i)=std::vector<couple>::operator[](i+1);
            this->pop_back();
        }

        auto erase(typename std::vector<couple>::iterator it)
        {
            auto Z=this->begin();
            while(Z!=this->end() && Z->first!=it->first)
                ++Z;
            for(int i=std::distance(this->begin(),Z);i+1<this->size();i++)
                std::vector<couple>::operator[](i)=std::vector<couple>::operator[](i+1);
            this->pop_back();
            return Z;
        }


        auto find(const K &x)
        {
            auto it=this->begin();
            while(it!=this->end())
                if(it->first==x)
                    return it;
                else ++it;
            return it;
        }

        auto emplace(const K&k,const V&v)
        {
            auto it=find(k);
            if(it!=this->end()) {
                it->second = v;
                return std::make_pair(it,false);
            }
            else {
                this->template emplace_back(std::numeric_limits<K>::max(),0);
                for (int i = this->size() - 1; i >= 0; i--) {
                    if (i==0 || k > std::vector<couple>::operator[](i-1).first) {
                        std::vector<couple>::operator[](i).first = k;
                        std::vector<couple>::operator[](i).second = v;
                        return std::make_pair(this->begin() + i, true);
                    }
                    std::vector<couple>::operator[](i) = std::vector<couple>::operator[](i-1);
                }
            }
        }

        V& operator[](std::size_t x)
        {
            auto it=find(x);
            if(it!=(*this).end())
                return it->second;
            else {
                auto [Z,_]=emplace(x, V{});
                return Z->second;
            }
        }
    };
}


#endif //CPLIBRARY_FLAT_STRUCTURES_H
