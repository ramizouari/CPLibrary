//
// Created by ramizouari on 27/10/23.
//

#ifndef CPLIBRARY_GENERAL_H
#define CPLIBRARY_GENERAL_H
#include <span>
#include <variant>
#include <vector>

template<typename T>
struct view_or_value
{
    std::variant<std::vector<T>,std::span<const T>> value;
    view_or_value(std::vector<T> _value):value(std::move(_value)){}
    view_or_value(std::span<const T> _value):value(std::move(_value)){}
    std::span<const T> get() const
    {
        if(value.index()==0) return std::span<const T>(std::get<0>(value).data(),std::get<0>(value).size());
        return std::get<1>(value);
    }
    operator std::span<const T>() const
    {
        return get();
    }

    auto begin() const
    {
        return get().begin();
    }
    auto end() const
    {
        return get().end();
    }

};



namespace graph
{
    template<typename T>
    struct AbstractGraph
    {
        [[nodiscard]] virtual int size() const = 0;
        virtual view_or_value<T> adjacentNodes(const T&u, bool direction) const =0;
        virtual view_or_value<T> adjacentNodes(const T&u) const =0;
        virtual view_or_value<T> nodes() const =0;
    };

    template<typename T,typename Weight>
    struct AbstractWeightedGraph
    {
        [[nodiscard]] virtual int size() const = 0;
        using AdjacentType=std::pair<T,Weight>;
        virtual view_or_value<AdjacentType> adjacentNodes(const T&u, bool direction) const =0;
        virtual view_or_value<AdjacentType> adjacentNodes(const T&u) const =0;
        virtual view_or_value<T> nodes() const =0;
    };
}

#endif //CPLIBRARY_GENERAL_H
