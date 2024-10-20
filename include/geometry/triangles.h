//
// Created by ramizouari on 15/04/24.
//

#ifndef CPLIBRARY_TRIANGLES_H
#define CPLIBRARY_TRIANGLES_H
#include "point.h"
#include "shapes.h"
namespace cp::geometry
{
    template<std::floating_point Float>
    struct triangle
    {
        std::array<point<Float>,3> points;
        triangle(const point<Float> &A,const point<Float> &B, const point<Float> & C) : points({A,B,C})
        {
        }

        circle<Float> incircle() const
        {
            std::array<Float,3> alpha;
            point<Float> center;
            Float perimeter=0,radius;
            for(int i=0;i<3;i++)
            {
                alpha[i] = distance(points[(i + 1) % 3], points[(i + 2) % 3]);
                perimeter+=alpha[i];
                center+=alpha[i]*points[i];
            }
            center /= perimeter;
            Float semi_perimeter = perimeter/2;
            radius =1;
            for(int i=0;i<3;i++)
                radius *= (semi_perimeter-points[i]);
            radius /=semi_perimeter;
            radius = std::sqrt(radius);
            return circle(center,radius);
        }

        circle<Float> circumcircle() const
        {

        }

    };
}

#endif //CPLIBRARY_TRIANGLES_H
