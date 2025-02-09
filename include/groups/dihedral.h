//
// Created by ramizouari on 21/12/24.
//

#ifndef DIHEDRAL_H
#define DIHEDRAL_H
#include <algebra/types.h>

namespace cp::groups {
    template<integer order>
    struct dihedral {
        integer index{};
        dihedral() =default;
        dihedral(integer index) : index((index%order+order)%order) {}

        bool direction() const {
            return index & 1;
        }

        bool operator==(const dihedral& rhs) const = default;

        dihedral operator*(const dihedral& rhs) const {
            integer res;
            if (direction())
                res = 2*((index / 2 + rhs.index /2 )%order);
            else
                res =2*((order + index / 2 - rhs.index /2 )%order);
            return dihedral(2*res+(rhs.direction()^direction()));
        }

        dihedral& operator*=(const dihedral& rhs)  {
            return *this = *this * rhs;
        }

        dihedral inv() const {
            return dihedral(2*order - 2*(index/2));
        }

        dihedral& operator/(const dihedral& rhs)  {
            return *this * rhs.inv();
        }
        dihedral& operator/=(const dihedral& rhs)  {
            return *this *= rhs.inv();
        }

    };
}

#endif //DIHEDRAL_H
