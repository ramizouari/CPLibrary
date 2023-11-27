//
// Created by ramizouari on 27/11/23.
//

#ifndef CPLIBRARY_VIEW_H
#define CPLIBRARY_VIEW_H
#include <array>
#include <utility>
#include <vector>
#include <numeric>
namespace cp::linalg
{
    inline constexpr std::size_t dynamic_extent = -1;
    template<typename R,std::size_t Rank>
    struct tensor_subview;
    template<typename R,std::size_t Rank>
    struct to_dynamic_view_t;
    template<typename R,std::size_t Rank>
    struct tensor_view
    {
        virtual R& at(std::array<std::size_t,Rank> indexes) = 0;
        virtual const R& at(std::array<std::size_t,Rank> indexes) const = 0;
        virtual ~tensor_view()= default;
        template<typename ...Args>
        R& at(Args... args)
        {
            return at(std::array<std::size_t,Rank>{args...});
        }
        template<typename ...Args>
        const R& at(Args... args) const
        {
            return at(std::array<std::size_t,Rank>{args...});
        }

        template<typename ...Args>
        R& operator()(Args... args)
        {
            return at(std::array<std::size_t,Rank>{static_cast<std::size_t>(args)...});
        }

        template<typename ...Args>
        const R& operator()(Args... args) const
        {
            return at(std::array<std::size_t,Rank>{static_cast<std::size_t>(args)...});
        }


        const R& operator()(std::array<std::size_t,Rank> args) const
        {
            return at(std::move(args));
        }

        R& operator()(std::array<std::size_t,Rank> args)
        {
            return at(std::move(args));
        }


        virtual std::array<std::size_t,Rank> shape() const = 0;
        static constexpr std::size_t rank()
        {
            return Rank;
        }
        virtual std::size_t size() const
        {
            auto s=shape();
            return std::accumulate(s.begin(),s.end(),1,std::multiplies<>());
        }

        struct iterator
        {
            tensor_view<R,Rank> &src;
            std::array<std::size_t,Rank> indexes;
            bool is_end=false;
            iterator(tensor_view<R,Rank> &src,std::array<std::size_t,Rank> indexes,bool is_end=false):src(src),indexes(indexes),is_end(is_end){}
            iterator& operator++()
            {
                for(int i=Rank-1;i>=0;i--)
                {
                    if(indexes[i]+1<src.shape()[i])
                    {
                        indexes[i]++;
                        return *this;
                    }
                    else
                        indexes[i]=0;
                }
                is_end=true;
                return *this;
            }
            iterator operator++(int)
            {
                iterator tmp=*this;
                ++(*this);
                return tmp;
            }
            bool operator==(const iterator& rhs) const
            {
                return is_end==rhs.is_end && indexes==rhs.indexes;
            }

            bool operator!=(const iterator& rhs) const
            {
                return !(*this==rhs);
            }

            R& operator*()
            {
                return src.at(indexes);
            }

            R* operator->()
            {
                return &src.at(indexes);
            }
        };
        iterator begin()
        {
            return iterator(*this,std::array<std::size_t,Rank>{},false);
        }
        iterator end()
        {
            auto shape=this->shape();
            return iterator(*this,std::array<std::size_t,Rank>{},true);
        }
        virtual tensor_subview<R,Rank> slice(std::array<std::size_t,Rank> start,std::array<std::size_t,Rank> end);
    };


    template <typename R,std::size_t Rank>
    struct tensor_subview:tensor_view<R,Rank>
    {
        tensor_view<R,Rank> &src;
        std::array<std::size_t,Rank> start,end;
        tensor_subview(tensor_view<R,Rank> &src,std::array<std::size_t,Rank> start,std::array<std::size_t,Rank> end):src(src),start(start),end(end){}

        R& at(std::array<std::size_t,Rank> indexes) override
        {
            std::array<std::size_t,Rank> new_indexes;
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i] + start[i];
            return src.at(new_indexes);
        }
        const R& at(std::array<std::size_t,Rank> indexes) const override
        {
            std::array<std::size_t,Rank> new_indexes;
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i] + start[i];
            return src.at(new_indexes);
        }
        std::array<std::size_t,Rank> shape() const override
        {
            std::array<std::size_t,Rank> new_shape;
            for(int i=0;i<Rank;i++)
                new_shape[i]=end[i]-start[i];
            return new_shape;
        }
        tensor_subview<R,Rank> slice(std::array<std::size_t,Rank> start,std::array<std::size_t,Rank> end) override
        {
            std::array<std::size_t,Rank> new_start,new_end;
            for(int i=0;i<Rank;i++)
            {
                new_start[i]=this->start[i]+start[i];
                new_end[i]=this->start[i]+end[i];
            }
            return tensor_subview<R,Rank>(src,new_start,new_end);
        }
    };

    template<typename R,std::size_t Rank>
    tensor_subview<R,Rank> tensor_view<R,Rank>::slice(std::array<std::size_t,Rank> start,std::array<std::size_t,Rank> end)
    {
        return tensor_subview<R,Rank>(*this,start,end);
    }

    template<typename R>
    struct tensor_view<R,dynamic_extent>
    {
        virtual ~tensor_view()= default;

        template<typename ...Args>
        R& at(Args... args)
        {
            return at(std::vector<std::size_t>{args...});
        }
        template<typename ...Args>
        const R& at(Args... args) const
        {
            return at(std::vector<std::size_t>{args...});
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
        virtual R& at(std::vector<std::size_t> indexes) = 0;
        virtual const R& at(std::vector<std::size_t> indexes) const = 0;

        [[nodiscard]] virtual std::vector<std::size_t> shape() const = 0;
        virtual std::size_t rank() const
        {
            return shape().size();
        }
        virtual std::size_t size() const
        {
            auto s=shape();
            return std::accumulate(s.begin(),s.end(),1,std::multiplies<>());
        }

        struct iterator
        {
            tensor_view<R,dynamic_extent> &src;
            std::vector<std::size_t> indexes;
            bool is_end=false;
            iterator(tensor_view<R,dynamic_extent> &src,std::vector<std::size_t> indexes,bool is_end=false):src(src),indexes(std::move(indexes)),is_end(is_end){}
            iterator& operator++()
            {
                for(int i=src.rank()-1;i>=0;i--)
                {
                    if(indexes[i]+1<src.shape()[i])
                    {
                        indexes[i]++;
                        return *this;
                    }
                    else
                        indexes[i]=0;
                }
                is_end=true;
                return *this;
            }
            iterator operator++(int)
            {
                iterator tmp=*this;
                ++(*this);
                return tmp;
            }
            bool operator==(const iterator& rhs) const
            {
                return is_end==rhs.is_end && indexes==rhs.indexes;
            }

            bool operator!=(const iterator& rhs) const
            {
                return !(*this==rhs);
            }

            R& operator*()
            {
                return src.at(indexes);
            }
        };
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
        virtual tensor_subview<R,dynamic_extent> slice(std::vector<std::size_t> start,std::vector<std::size_t> end);
    };

    template<typename R>
    struct vector_view : public tensor_view<R,1>
    {
        R* m_data;
        std::size_t m_size;
        vector_view(R* _data,std::size_t _size):m_data(_data),m_size(_size){}
        vector_view(std::vector<R> &v):m_data(v.data()),m_size(v.size()){}
        template<std::size_t N>
        vector_view(std::array<R,N> &v):m_data(v.data()),m_size(v.size()){}
        std::size_t size() const override
        {
            return m_size;
        }
        virtual R& at(std::size_t i)
        {
            return m_data[i];
        }
        virtual const R& at(std::size_t i) const
        {
            return m_data[i];
        }
        virtual ~vector_view()= default;
        R& at(std::array<std::size_t,1> indexes) override
        {
            return at(indexes[0]);
        }
        const R& at(std::array<std::size_t,1> indexes) const override
        {
            return at(indexes[0]);
        }
        std::array<std::size_t,1> shape() const override
        {
            return {size()};
        }
        tensor_subview<R,1> slice(std::array<std::size_t,1> start,std::array<std::size_t,1> end) override
        {
            return tensor_subview<R,1>(*this,{start[0]},{end[0]});
        }
        vector_view& operator=(const std::vector<R>& O)
        {
            for(int i=0;i<size();i++)
                at(i)=O.at(i);
            return *this;
        }
    };

    template<typename R,std::size_t Rank>
    struct to_dynamic_view_t : public tensor_view<R,dynamic_extent>
    {
        tensor_view<R,Rank> &src;
        to_dynamic_view_t(tensor_view<R,Rank> &src):src(src){}
        R& at(std::vector<std::size_t> indexes) override
        {
            std::array<std::size_t,Rank> new_indexes;
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i];
            return src.at(new_indexes);
        }
        const R& at(std::vector<std::size_t> indexes) const override
        {
            std::array<std::size_t,Rank> new_indexes;
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i];
            return src.at(new_indexes);
        }
        std::vector<std::size_t> shape() const override
        {
            return std::vector<std::size_t>(src.shape().begin(),src.shape().end());
        }
        tensor_subview<R,dynamic_extent> slice(std::vector<std::size_t> start,std::vector<std::size_t> end) override
        {
            std::array<std::size_t,Rank> new_start,new_end;
            for(int i=0;i<Rank;i++)
            {
                new_start[i]=start[i];
                new_end[i]=end[i];
            }
            return tensor_subview<R,dynamic_extent>(src,new_start,new_end);
        }
    };

    template<typename R>
    struct to_dynamic_view_t<R,dynamic_extent>: public tensor_view<R,dynamic_extent>{};

    template<typename R,std::size_t Rank>
    struct to_static_view_t : public tensor_view<R,Rank>
    {
        tensor_view<R,dynamic_extent> &src;
        std::array<std::size_t,Rank> _shape;
        to_static_view_t(tensor_view<R,dynamic_extent> &src,std::array<std::size_t,Rank> _shape):src(src),_shape(_shape){}
        to_static_view_t(tensor_view<R,dynamic_extent> &src):src(src)
        {
            auto s=src.shape();
            for(int i=0;i<std::min(Rank,src.rank());i++)
                _shape[i]=s[i];
        }
        R& at(std::array<std::size_t,Rank> indexes) override
        {
            std::vector<std::size_t> new_indexes(Rank);
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i];
            return src.at(new_indexes);
        }
        const R& at(std::array<std::size_t,Rank> indexes) const override
        {
            std::vector<std::size_t> new_indexes(Rank);
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i];
            return src.at(new_indexes);
        }
        std::array<std::size_t,Rank> shape() const override
        {
            return _shape;
        }
    };

    template<typename R>
    struct tensor_subview<R,dynamic_extent>:tensor_view<R,dynamic_extent>
    {
        tensor_view<R,dynamic_extent> &src;
        std::vector<std::size_t> start,end;
        tensor_subview(tensor_view<R,dynamic_extent> &src,std::vector<std::size_t> start,std::vector<std::size_t> end):src(src),start(start),end(end){}

        R& at(std::vector<std::size_t> indexes) override
        {
            std::vector<std::size_t> new_indexes(src.rank());
            for(int i=0;i<src.rank();i++)
                new_indexes[i]=indexes[i] + start[i];
            return src.at(new_indexes);
        }
        const R& at(std::vector<std::size_t> indexes) const override
        {
            std::vector<std::size_t> new_indexes(src.rank());
            for(int i=0;i<src.rank();i++)
                new_indexes[i]=indexes[i] + start[i];
            return src.at(new_indexes);
        }
        std::vector<std::size_t> shape() const override
        {
            std::vector<std::size_t> new_shape(src.rank());
            for(int i=0;i<src.rank();i++)
                new_shape[i]=end[i]-start[i];
            return new_shape;
        }
        tensor_subview<R,dynamic_extent> slice(std::vector<std::size_t> start,std::vector<std::size_t> end) override
        {
            std::vector<std::size_t> new_start(src.rank()),new_end(src.rank());
            for(int i=0;i<src.rank();i++)
            {
                new_start[i]=this->start[i]+start[i];
                new_end[i]=this->start[i]+end[i];
            }
            return tensor_subview<R,dynamic_extent>(src,new_start,new_end);
        }
    };

    template<typename R>
    tensor_subview<R,dynamic_extent> tensor_view<R,dynamic_extent>::slice(std::vector<std::size_t> start,std::vector<std::size_t> end)
    {
        return tensor_subview<R,dynamic_extent>(*this,start,end);
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

#endif //CPLIBRARY_VIEW_H
