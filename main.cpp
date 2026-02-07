#include <chrono>
#include <iostream>
#include <random>

#include "algebra/binary_operation.h"
#include "data_structures/fixed/multi_dimensional.h"
#include "data_structures/fixed/recursive.h"
#include "data_structures/fixed/segment_tree.h"

#include "data_structures/statistic_tree/key_statistics.h"
#include "nt/dirichelet.h"
#include "include/tensors/tensor.h"
#include "nt/modular_arithmetic.h"
#include "rings/identities.h"
#include "linalg/matrix.h"

// using IK = double;
using Ring = cp::d_cyclic;
enum arithmetic_expression_rule
{
    E,
    M,
    T,
    N,
    COUNT,
    Expression=E,
    Multiplier=M,
    Term=T,
    Number=N,
};

using IK = cp::d_cyclic;

using triplet = std::array<std::size_t,3>;

struct triplet_right_op : cp::binary_operation<triplet>
{
    triplet reduce(const triplet &a, const triplet &b) const override
    {
        if (b==neutral)
            return a;
        return b;
    }

    inline static triplet neutral = {0,0,0};


};


constexpr auto W=64, Q=1000;

struct triplet_left_op : cp::binary_operation<triplet>
{
    triplet reduce(const triplet &a, const triplet &b) const override
    {
        if (a==neutral)
            return b;
        return a;
    }


    inline static triplet neutral = {W,W,W};


};

int main()
{
    cp::linalg::flat_tensor<triplet,3> Z({W,W,W});
    for (std::size_t i=0;i<W;i++) for (std::size_t j=0;j<W;j++) for (std::size_t k=0;k<W;k++)
        Z.at({i,j,k})={i,j,k};
    cp::data_structures::fixed::multidimensional_segment_tree<triplet_left_op,3> ST(Z);
    std::size_t a,b,c,d,e,f;
    while (std::cin >> a >> b >> c >> d >> e >> f)
    {
        auto [i,j,k] = ST.query({a,b,c},{d,e,f});
        std::cout << i << ' ' << j << ' ' << k << '\n';
    }

}
