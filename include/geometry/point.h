//
// Created by ramizouari on 15/04/24.
//

#ifndef CPLIBRARY_POINT_H
#define CPLIBRARY_POINT_H
#include <complex>
#include "algebra/abstract_algebra.h"
#include "linear_algebra/vector.h"

namespace cp::geometry
{
    template<std::floating_point Float>
    using vector = linalg::s_vector<Float,2>;

    template<std::floating_point Float>
    struct point
    {
        vector<Float> coordinates;
        point() = default;
        point(Float x, Float y) : coordinates({x,y}){}
        explicit point(vector<Float> && coordinates) : coordinates(std::move(coordinates)){}
        explicit point(vector<Float> &coordinates) : coordinates(coordinates){}
        explicit point(std::complex<Float> z) : coordinates({z.real(),z.imag()}){}
        Float& x()
        {
            return coordinates[0];
        }

        const Float& x() const
        {
            return coordinates[0];
        }

        Float& y()
        {
            return coordinates[1];
        }

        const Float& y() const
        {
            return coordinates[1];
        }

        point& operator+=(const vector<Float> &z)
        {
            coordinates+=z;
            return *this;
        }

        point& operator-=(const vector<Float> &z)
        {
            coordinates-=z;
            return *this;
        }

        point operator+(const vector<Float> & z) const
        {
            return point(coordinates+z);
        }

        point operator-(const vector<Float> & z) const
        {
            return point(coordinates-z);
        }

        point& operator*=(Float scale)
        {
            coordinates*=scale;
            return *this;
        }

        point& operator/=(Float scale)
        {
            coordinates/=scale;
            return *this;
        }

        point operator*(Float scale) const
        {
            return point(scale * coordinates);
        }

        point operator/(Float scale) const
        {
            return point(coordinates / scale);
        }

        explicit operator std::complex<Float>() const
        {
            return std::complex<Float>(x(),y());
        };

        vector<Float> operator-(const point & o) const
        {
            return coordinates - o.coordinates;
        }
    };

    template<std::floating_point Float>
    Float norm(const vector<Float> & u)
    {
        return std::hypot(u[0],u[1]);
    }

    template<std::floating_point Float>
    Float distance(const point<Float> & p, const point<Float> & q)
    {
        return norm(q-p);
    }

    template<std::floating_point Float>
    point<Float> operator*(Float scale, point<Float> & p)
    {
        return  p * scale;
    }
}

#endif //CPLIBRARY_POINT_H
