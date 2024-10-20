//
// Created by ramizouari on 15/04/24.
//

#ifndef CPLIBRARY_SHAPES_H
#define CPLIBRARY_SHAPES_H
#include "point.h"
#include <optional>
#include <variant>

namespace cp::geometry
{
    template<std::floating_point Float>
    struct circle
    {
        point<Float> center;
        Float radius;
        circle(const point<Float> &center, Float radius) : center(center),radius(radius){}
    };


    template<std::floating_point Float>
    vector<Float> normal(const vector<Float> &p)
    {
        return vector<Float>({-p[1],p[0]});
    }

    template<std::floating_point Float>
    struct line
    {
        point<Float> p;
        vector<Float> u;

        vector<Float> normal() const
        {
            return geometry::normal(u);
        }
    };

    template<std::floating_point Float>
    struct line_segment
    {
        point<Float> p,q;
        line_segment(const point<Float> &p, const point<Float> & q) : p(p),q(q){}
    };

    template<std::floating_point Float>
    line<Float> bisector(const line_segment<Float> &P)
    {
        auto dP= P.q - P.p;
        return line(P.p + dP/2, normal(dP));
    }

    template<std::floating_point Float>
    Float det(const vector<Float> & u, const vector<Float> & v)
    {
        return u[0]*v[1] - u[1]*v[0];
    }

    template<std::floating_point Float>
    bool collinear(const vector<Float> &p1, const vector<Float> &p2, const vector<Float> & p3, Float eps= cp::epsilon)
    {
        return det(p3 - p1,p2-p1) < eps;
    }

    template<std::floating_point Float>
    using ll_intersection  = std::variant<point<Float>,line<Float>,std::nullopt_t>;

    template<std::floating_point Float>
    ll_intersection<Float> intersection(const line<Float> &P, const line<Float> & Q, Float eps = cp::epsilon)
    {
        if(std::abs(det(P.u,Q.u)) < eps)
        {
            if(det(P.p - Q.p,P.u) < eps)
                return P;
            else
                return std::nullopt;
        }

    }

    template<std::floating_point Float>
    struct simple_polygon
    {
        std::vector<point<Float>> vertices;
        Float signed_area() const
        {
            Float A=0;
            int n=vertices.size();
            for(int i=0;i<n;i++)
                A+=det(vertices[i],vertices[(i+1)%n]);
            A/=2;
        }

        Float area() const
        {
            return std::abs(signed_area());
        }

    };
}

#endif //CPLIBRARY_SHAPES_H
