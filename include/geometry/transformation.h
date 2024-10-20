//
// Created by ramizouari on 15/04/24.
//

#ifndef CPLIBRARY_TRANSFORMATION_H
#define CPLIBRARY_TRANSFORMATION_H
#include "point.h"
#include "linear_algebra/matrix.h"
#include "shapes.h"

namespace cp::geometry
{

    enum direction : bool
    {
        direct, opposite, indirect = opposite
    };

    template<std::floating_point Float>
    class transformation
    {
    public:
        virtual point<Float> transform(const point<Float> & p) const = 0;
        point<Float> operator()(const point<Float> & p) const
        {
            return transform(p);
        }

        [[nodiscard]] virtual direction orientation() const
        {
            point<Float> e1(1,0),e2(0,1);
            return static_cast<direction>(det(transform(e1), transform(e2)) < 0);
        }
    };

    template<std::floating_point Float>
    struct invertible_transformation : public transformation<Float>
    {
        virtual point<Float> inverse_transform(const point<Float> &p ) const = 0;
    };

    template<std::floating_point Float>
    class affine_transformation : transformation<Float>
    {
        cp::linalg::s_matrix<Float,2,2> M;
        point<Float> f0;
    public:
        affine_transformation(const point<Float> &f0,const point<Float> &f1, const point<Float> &f2):f0(f0)
        {
            for(int i=0;i<2;i++)
            {
                M[i][0] = f1 - f0;
                M[i][1] = f2 - f0;
            }
        }

        point<Float> transform(const point<Float> & p) const override
        {
            return f0 + M * p.coordinates;
        }

        [[nodiscard]] direction orientation() const override
        {
            return static_cast<direction>(M.det() < 0);
        }
    };

    template<std::floating_point Float>
    class linear_transformation : transformation<Float>
    {
        cp::linalg::s_matrix<Float,2,2> M;
    public:
        linear_transformation(const point<Float> &f1, const point<Float> &f2)
        {
            for(int i=0;i<2;i++)
            {
                M[i][0] = f1.coordinates[i];
                M[i][1] = f2.coordinates[i];
            }
        }

        point<Float> transform(const point<Float> & p) const override
        {
            return point(M * p.coordinates);
        }

        [[nodiscard]] direction orientation() const override
        {
            return static_cast<direction>(M.det() < 0);
        }

        operator affine_transformation<Float> () const
        {
            point<Float> f1,f2;
            for(int i=0;i<2;i++)
            {
                f1.coordinates[i]=M[i][0];
                f2.coordinates[i]=M[i][1];
            }
            return affine_transformation<Float>(point<Float>{},f1,f2);
        }
    };

    template<std::floating_point Float>
    struct rotation : invertible_transformation<Float>
    {
        vector<Float> center;
        Float angle;
        explicit rotation(Float angle, const vector<Float> & center={}) : center(center),angle(angle){}

        point<Float> transform(const point<Float> & p) const override
        {
            auto u=p-center;
            std::complex<Float> z(u[0],u[1]);
            return center + point<Float>(std::polar(1,angle) * z);
        }

        point<Float> inverse_transform(const point<Float> & p) const override
        {
            auto u=p-center;
            std::complex<Float> z(u[0],u[1]);
            return center + point<Float>(std::polar(1,-angle) * z);
        }

        rotation inverse() const
        {
            return rotation(-angle,center);
        }

        rotation operator-() const
        {
            return inverse();
        }

        [[nodiscard]] direction orientation() const override
        {
            return direct;
        }
    };

    template<std::floating_point Float>
    struct homothety : invertible_transformation<Float>
    {
        vector<Float> center;
        Float scale;
        explicit homothety(Float scale, const vector<Float> & center={}) : center(center),scale(scale){}

        point<Float> transform(const point<Float> & p) const override
        {
            return center + scale*(p-center);
        }

        point<Float> inverse_transform(const point<Float> & p) const override
        {
            return center + 1/scale*(p-center);
        }

        homothety inverse() const
        {
            return rotation(1/scale,center);
        }


        [[nodiscard]] direction orientation() const override
        {
            return direct;
        }
    };

    template<std::floating_point Float>
    struct similitude : transformation<Float>
    {
        std::complex<Float> a,b;
        direction o;
        similitude(Float scale,Float angle, vector<Float> t,
                   direction o = direct) : a(std::polar(scale,angle)), b(t[0],t[1]),o(o){}
        point<Float> transform(const point<Float> & p) const override
        {
            std::complex<Float> z(p.x(),p.y());
            if(o==opposite) z=std::conj(z);
            return point<Float>(a*z+b);
        }

        [[nodiscard]] direction orientation() const
        {
            return o;
        }
    };

    template<std::floating_point Float>
    struct direct_similitude : invertible_transformation<Float>
    {
        point<Float> center;
        Float scale,angle;
        direct_similitude(Float scale, Float angle,point<Float> center= {}) : scale(scale),angle(angle),center(center)
        {
            if(scale <= 0)
                throw std::runtime_error("Invalid similitude: scale<=0");
        }
        point<Float> transform(const point<Float> & p) const override
        {
            auto u=p-center;
            std::complex<Float> z(u[0],u[1]);
            return center + point<Float>(std::polar(scale,angle) * z);
        }

        point<Float> inverse_transform(const point<Float> & p) const override
        {
            auto u=p-center;
            std::complex<Float> z(u[0],u[1]);
            return center + point<Float>(std::polar(1/scale,-angle) * z);
        }

        direct_similitude inverse() const
        {
            return direct_similitude(1/scale, -angle, center);
        }

        [[nodiscard]] direction orientation() const override
        {
            return direct;
        }
    };
}

#endif //CPLIBRARY_TRANSFORMATION_H
