//
// Created by ramizouari on 27/11/23.
//

#ifndef CPLIBRARY_VIEW_H
#define CPLIBRARY_VIEW_H
#include <array>
#include <utility>
#include <vector>
#include <numeric>

#include "utils.h"
#include "indexing.h"


namespace cp::tensors
{
    template<typename R,std::size_t Rank>
    struct tensor_subview;
    template<typename R,std::size_t Rank>
    struct tensor_view
    {
        using index_array = std::array<std::size_t, Rank>;

        virtual R& at(index_array indexes) = 0;
        virtual const R& at(index_array indexes) const = 0;
        virtual ~tensor_view()= default;
        template<typename ...Args>
        R& at(Args... args)
        {
            return at(index_array{args...});
        }
        template<typename ...Args>
        const R& at(Args... args) const
        {
            return at(index_array{args...});
        }

        template<typename ...Args>
        R& operator()(Args... args)
        {
            return at(index_array{static_cast<std::size_t>(args)...});
        }

        template<typename ...Args>
        const R& operator()(Args... args) const
        {
            return at(index_array{static_cast<std::size_t>(args)...});
        }


        const R& operator()(index_array args) const
        {
            return at(std::move(args));
        }

        R& operator()(index_array args)
        {
            return at(std::move(args));
        }

        virtual index_array shape() const = 0;
        static constexpr std::size_t rank()
        {
            return Rank;
        }
        virtual std::size_t size() const
        {
            auto s=shape();
            return std::accumulate(s.begin(),s.end(),1,std::multiplies());
        }

        template<bool is_const = false>
        struct iterator_t
        {
            using tensor_view_t = std::conditional_t<is_const,const tensor_view,tensor_view>;
            using type_t = std::conditional_t<is_const,const R,R>;


            tensor_view_t &src;
            index_array indexes;
            bool is_end=false;
            iterator_t(tensor_view_t &src,index_array indexes,bool is_end=false):src(src),indexes(indexes),is_end(is_end){}
            iterator_t& operator++()
            {
                auto shape=src.shape();
                is_end=rincrement(indexes,shape);
                return *this;
            }
            iterator_t operator++(int)
            {
                iterator_t tmp=*this;
                ++*this;
                return tmp;
            }
            bool operator==(const iterator_t& rhs) const
            {
                return is_end==rhs.is_end && indexes==rhs.indexes;
            }

            type_t& operator*()
            {
                return src.at(indexes);
            }

            const type_t& operator*() const
            {
                return src.at(indexes);
            }

            type_t* operator->()
            {
                return &src.at(indexes);
            }
        };

        using iterator = iterator_t<false>;
        using const_iterator = iterator_t<true>;

        iterator begin()
        {
            return iterator(*this,index_array{},false);
        }
        iterator end()
        {
            auto shape=this->shape();
            return iterator(*this,index_array{},true);
        }

        const_iterator begin() const
        {
            return const_iterator(*this,index_array{},false);
        }
        const_iterator end() const
        {
            auto shape=this->shape();
            return const_iterator(*this,index_array{},true);
        }

        virtual tensor_subview<R,Rank> slice(index_array start,index_array end);
        virtual tensor_subview<R,Rank> slice(index_array start,index_array end,index_array step);
    };


    template <typename R,std::size_t Rank>
    struct tensor_subview final: tensor_view<R,Rank>
    {
        using typename tensor_view<R,Rank>::index_array;

        tensor_view<R,Rank> &src;
        index_array m_start,m_end,m_step;
        tensor_subview(tensor_view<R,Rank> &src,index_array start,index_array end):src(src),m_start(start),m_end(end)
        {
            std::fill(m_step.begin(),m_step.end(),1);
        }
        tensor_subview(tensor_view<R,Rank> &src,index_array start,index_array end,
                       index_array step):src(src),m_start(start),m_end(end),m_step(step){}

        ~tensor_subview() override = default;


        R& at(index_array indexes) override
        {
            index_array new_indexes;
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i]*m_step[i] + m_start[i];
            return src.at(new_indexes);
        }

        const R& at(index_array indexes) const override
        {
            index_array new_indexes;
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i]*m_step[i] + m_start[i];
            return src.at(new_indexes);
        }

        index_array shape() const override
        {
            index_array new_shape;
            for(int i=0;i<Rank;i++)
                new_shape[i]=(m_end[i]-m_start[i]+m_step[i]-1)/m_step[i];
            return new_shape;
        }
        tensor_subview<R,Rank> slice(index_array start,index_array end) override
        {
            index_array new_start,new_end;
            for(int i=0;i<Rank;i++)
            {
                new_start[i]=this->m_start[i]+start[i]*m_step[i];
                new_end[i]=this->m_start[i]+end[i]*m_step[i];
            }
            return tensor_subview(src,new_start,new_end,m_step);
        }
        tensor_subview<R,Rank> slice(index_array start,index_array end,index_array step) override
        {
            index_array new_start,new_end,new_step;
            for(int i=0;i<Rank;i++)
            {
                new_start[i]=this->m_start[i]+start[i]*m_step[i];
                new_end[i]=this->m_start[i]+end[i]*m_step[i];
                new_step[i]=this->m_step[i]*step[i];
            }
            return tensor_subview(src,new_start,new_end,new_step);
        }
    };

    template<typename R,std::size_t Rank>
    tensor_subview<R,Rank> tensor_view<R,Rank>::slice(index_array start,index_array end)
    {
        return tensor_subview<R,Rank>(*this,start,end);
    }

    template<typename R,std::size_t Rank>
    tensor_subview<R,Rank> tensor_view<R,Rank>::slice(index_array start,index_array end,index_array step)
    {
        return tensor_subview<R,Rank>(*this,start,end,step);
    }

    template<typename R>
    struct tensor_view<R,dynamic_extent>
    {
        using index_array = std::vector<std::size_t>;

        virtual ~tensor_view()= default;

        template<typename ...Args>
        R& at(Args... args)
        {
            return at(index_array{args...});
        }
        template<typename ...Args>
        const R& at(Args... args) const
        {
            return at(index_array{args...});
        }

        template<typename ...Args>
        R& operator()(Args... args)
        {
            return at(args...);
        }

        template<typename ...Args>
        const R& operator()(Args... args) const
        {
            return at(args...);
        }
        virtual R& at(index_array indexes) = 0;
        virtual const R& at(index_array indexes) const = 0;

        [[nodiscard]] virtual index_array shape() const = 0;
        virtual std::size_t rank() const
        {
            return shape().size();
        }
        virtual std::size_t size() const
        {
            auto s=shape();
            return std::accumulate(s.begin(),s.end(),1,std::multiplies());
        }

        template<bool is_const = false>
        struct iterator_t
        {
            using tensor_view_t = std::conditional_t<is_const,const tensor_view,tensor_view>;
            using type_t = std::conditional_t<is_const,const R,R>;

            tensor_view_t &src;
            index_array indexes;
            bool is_end=false;
            iterator_t(tensor_view_t &src,index_array indexes,bool is_end=false):src(src),indexes(std::move(indexes)),is_end(is_end){}
            iterator_t& operator++()
            {
                is_end=rincrement(indexes,src.shape());
                return *this;
            }

            iterator_t operator++(int)
            {
                iterator_t tmp=*this;
                ++*this;
                return tmp;
            }
            bool operator==(const iterator_t& rhs) const
            {
                return is_end==rhs.is_end && indexes==rhs.indexes;
            }


            type_t& operator*()
            {
                return src.at(indexes);
            }

            const type_t& operator*() const
            {
                return src.at(indexes);
            }

            type_t* operator->()
            {
                return &src.at(indexes);
            }
        };
        using iterator = iterator_t<false>;
        using const_iterator = iterator_t<true>;
        iterator begin()
        {
            auto zeros=shape();
            for(auto &x:zeros)
                x=0;
            return iterator(*this,zeros,false);
        }
        iterator end()
        {
            auto zeros=shape();
            for(auto &x:zeros)
                x=0;
            return iterator(*this,zeros,true);
        }

        const_iterator begin() const
        {
            auto zeros=shape();
            for(auto &x:zeros)
                x=0;
            return const_iterator(*this,zeros,false);
        }
        const_iterator end() const
        {
            auto zeros=shape();
            for(auto &x:zeros)
                x=0;
            return const_iterator(*this,zeros,true);
        }

        virtual tensor_subview<R,dynamic_extent> slice(index_array start,index_array end);
        virtual tensor_subview<R,dynamic_extent> slice(index_array start,index_array end,index_array step);
    };

    template<typename R>
    struct tensor_subview<R,dynamic_extent>:tensor_view<R,dynamic_extent>
    {
        using typename tensor_view<R,dynamic_extent>::index_array;
        tensor_view<R,dynamic_extent> &src;
        index_array m_start,m_end,m_step;
        tensor_subview(tensor_view<R,dynamic_extent> &src,index_array start,index_array end):src(src),m_start(std::move(start)),m_end(std::move(end)),m_step(src.rank())
        {
            std::fill(m_step.begin(),m_step.end(),1);
        }
        tensor_subview(tensor_view<R,dynamic_extent> &src,index_array start,index_array end,
                       index_array step):src(src),m_start(std::move(start)),m_end(std::move(end)),m_step(std::move(step)){}

        ~tensor_subview() override = default;

        R& at(index_array indexes) override
        {
            index_array new_indexes(src.rank());
            for(int i=0;i<src.rank();i++)
                new_indexes[i]=indexes[i] + m_start[i];
            return src.at(new_indexes);
        }
        const R& at(index_array indexes) const override
        {
            index_array new_indexes(src.rank());
            for(int i=0;i<src.rank();i++)
                new_indexes[i]=indexes[i] + m_start[i];
            return src.at(new_indexes);
        }
        [[nodiscard]] index_array shape() const override
        {
            index_array new_shape(src.rank());
            for(int i=0;i<src.rank();i++)
                new_shape[i]=m_end[i]-m_start[i];
            return new_shape;
        }
        tensor_subview slice(index_array start,index_array end) override
        {
            index_array new_start(src.rank()),new_end(src.rank());
            for(int i=0;i<src.rank();i++)
            {
                new_start[i]=this->m_start[i]+start[i];
                new_end[i]=this->m_start[i]+end[i];
            }
            return tensor_subview(src,new_start,new_end);
        }

        tensor_subview slice(index_array start,index_array end,index_array step) override
        {
            auto Rank=src.rank();
            index_array new_start(Rank),new_end(Rank),new_step(Rank);
            for(int i=0;i<Rank;i++)
            {
                new_start[i]=this->m_start[i]+start[i]*m_step[i];
                new_end[i]=this->m_start[i]+end[i]*m_step[i];
                new_step[i]=this->m_step[i]*step[i];
            }
            return tensor_subview(src,new_start,new_end,new_step);
        }

    };

    template<typename R,std::size_t Rank>
    struct to_dynamic_view_t : tensor_view<R,dynamic_extent>
    {
        using index_array = std::vector<std::size_t>;
        using index_original = std::array<std::size_t,Rank>;

        tensor_view<R,Rank> &src;
        explicit to_dynamic_view_t(tensor_view<R,Rank> &src):src(src){}
        R& at(index_array indexes) override
        {
            index_original new_indexes;
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i];
            return src.at(new_indexes);
        }
        const R& at(index_array indexes) const override
        {
            index_original new_indexes;
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i];
            return src.at(new_indexes);
        }
        [[nodiscard]] index_array shape() const override
        {
            return index_array(src.shape().begin(),src.shape().end());
        }
        tensor_subview<R,dynamic_extent> slice(index_array start,index_array end) override
        {
            index_array new_start(Rank),new_end(Rank);
            for(int i=0;i<Rank;i++)
            {
                new_start[i]=start[i];
                new_end[i]=end[i];
            }
            return tensor_subview<R,dynamic_extent>(*this,new_start,new_end);
        }

        tensor_subview<R,dynamic_extent> slice(index_array start,index_array end, index_array step) override
        {
            index_array new_start(Rank),new_end(Rank),new_step(Rank);
            for(int i=0;i<Rank;i++)
            {
                new_start[i]=start[i];
                new_end[i]=end[i];
                new_step[i]=step[i];
            }
            return tensor_subview<R,dynamic_extent>(*this,new_start,new_end,new_step);
        }
    };

    template<typename R>
    struct to_dynamic_view_t<R,dynamic_extent>: tensor_view<R,dynamic_extent>{};

    template<typename R,std::size_t Rank>
    struct to_static_view_t : tensor_view<R,Rank>
    {
        using index_array = std::array<std::size_t,Rank>;
        using index_original = std::vector<std::size_t>;

        tensor_view<R,dynamic_extent> &src;
        std::array<std::size_t,Rank> _shape;
        to_static_view_t(tensor_view<R,dynamic_extent> &src,index_array _shape):src(src),_shape(_shape){}
        explicit to_static_view_t(tensor_view<R,dynamic_extent> &src):src(src)
        {
            auto s=src.shape();
            for(int i=0;i<std::min(Rank,src.rank());i++)
                _shape[i]=s[i];
        }
        R& at(index_array indexes) override
        {
            index_original new_indexes(Rank);
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i];
            return src.at(new_indexes);
        }
        const R& at(index_array indexes) const override
        {
            index_original new_indexes(Rank);
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i];
            return src.at(new_indexes);
        }
        index_array shape() const override
        {
            return _shape;
        }
    };


    template<typename R>
    tensor_subview<R,dynamic_extent> tensor_view<R,dynamic_extent>::slice(index_array start,index_array end)
    {
        return tensor_subview<R,dynamic_extent>(*this,start,end);
    }

    template<typename R>
    tensor_subview<R,dynamic_extent> tensor_view<R,dynamic_extent>::slice(index_array start,index_array end,index_array step)
    {
        return tensor_subview<R,dynamic_extent>(*this,start,end,step);
    }

    template<typename R,std::size_t Rank>
    to_dynamic_view_t<R,Rank> to_dynamic_view(tensor_view<R,Rank> &src)
    {
        return to_dynamic_view_t<R,Rank>(src);
    }

    template<typename R,std::size_t Rank>
    to_static_view_t<R,Rank> to_static_view(tensor_view<R,dynamic_extent> &src,std::array<std::size_t,Rank> _shape)
    {
        return to_static_view_t<R,Rank>(src,_shape);
    }

    template<std::size_t Rank,typename R>
    to_static_view_t<R,Rank> to_static_view(tensor_view<R,dynamic_extent> &src)
    {
        return to_static_view_t<R,Rank>(src);
    }
}

namespace cp::linalg {
    using namespace tensors;
}

#endif //CPLIBRARY_VIEW_H
