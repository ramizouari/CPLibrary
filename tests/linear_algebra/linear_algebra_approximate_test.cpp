//
// Created by ramizouari on 11/04/2022.
//
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <complex>
#include <topology/topology.h>
#include <functional/zip.h>
using test_types=boost::mpl::list<double,long double,float,std::complex<long double>,std::complex<double>,std::complex<float>> ;
using real_types=boost::mpl::list<double,long double,float>;
#include "../print.h"
#include <random>

constexpr double err=1e-4;

constexpr int N=100;

constexpr double failure_tolerance=.05;

double rel_error(double a, double b)
{
    return std::max(std::abs(a-b)/std::abs(a),std::abs(a-b)/std::abs(b));
}

template<typename T,int n,int m>
bool approx(const s_matrix<T,n,m>&A,const s_matrix<T,n,m>&B,double err)
{
    for(auto&& [R1,R2]:zip(A,B)) for(auto &&[a,b]:zip(R1,R2))
            if(rel_error(a,b)>err)
                return false;
    return true;
}

template<typename T,int n,int m>
double l2_distance(const s_matrix<T,n,m>&A,const s_matrix<T,n,m>&B)
{
    double res=0;
    for(auto&& [R1,R2]:zip(A,B)) for(auto &&[a,b]:zip(R1,R2))
            res+=std::pow(a-b,2);
    return std::sqrt(res);
}

template<int n,int m>
s_matrix<real,n,m> random_matrix()
{
    std::random_device dev;
    std::mt19937_64 rng(dev());
    using T = real;
    constexpr double std = .175;
    std::normal_distribution<double> d(0, std);
    s_matrix<T, n, m> U;
    for (auto &R: U)
        for (auto &c: R)
            c = d(rng);
    return U;
}

template<int rank>
bool test_solver() {
    std::random_device dev;
    std::mt19937_64 rng(dev());
    constexpr int dimension = 100;
    using T = real;
    L2_inner_product<T, s_vector<T, 100>> innerProduct;
    constexpr double std = .175;
    constexpr double singularity = 1;
    std::normal_distribution<double> d(0, std), h(0, singularity);
    s_matrix<T, dimension, rank> U;
    s_matrix<T, rank, dimension> V;
    s_vector<T, dimension> v;
    for (auto &R: U)
        for (auto &c: R)
            c = d(rng);
    for (auto &R: V)
        for (auto &c: R)
            c = d(rng);
    auto A = U * V;
    for (auto &c: v)
        c = d(rng);
    return  innerProduct.metric(A*A.solve(v), v)< err;
}

BOOST_AUTO_TEST_SUITE(test_vector_real)
    BOOST_AUTO_TEST_SUITE(test_vector_constructor)
        BOOST_AUTO_TEST_CASE_TEMPLATE(test_vector_create,T,test_types)
        {
            s_vector<T, 5> v1({1, 2, 3, 4, 5});
            BOOST_CHECK_EQUAL(v1[0], T(1));
            BOOST_CHECK_EQUAL(v1[1], T(2));
            BOOST_CHECK_EQUAL(v1[2], T(3));
            BOOST_CHECK_EQUAL(v1[3], T(4));
            BOOST_CHECK_EQUAL(v1[4], T(5));
        }
    BOOST_AUTO_TEST_SUITE_END()
    BOOST_AUTO_TEST_SUITE(test_vector_operator)
        BOOST_AUTO_TEST_CASE_TEMPLATE(test_add_vector,T,test_types) {
            s_vector<T,5> v1({1,2,3,4,5});
            s_vector<T,5> v2({3,7,-1,-2,2});
            s_vector<T,5> v3({4,9,2,2,7});
            BOOST_CHECK_EQUAL(v1+v2,v3);
        }

        BOOST_AUTO_TEST_CASE_TEMPLATE(test_diff_vector,T,test_types) {
            s_vector<T,5> v1({1,2,3,4,5});
            s_vector<T,5> v2({3,7,-1,-2,2});
            s_vector<T,5> v3({-2,-5,4,6,3});
            BOOST_CHECK_EQUAL(v1-v2,v3);
        }

        BOOST_AUTO_TEST_CASE_TEMPLATE(test_scalar_vector,T,test_types)
        {
            s_vector<T,5> v1({1,2,3,4,5});
            s_vector<T,5> v2({3,6,9,12,15});
            BOOST_CHECK_EQUAL(T(3)*v1,v2);
        }

        BOOST_AUTO_TEST_CASE_TEMPLATE(test_vector_negate,T,test_types)
        {
            s_vector<T,5> v1({1,2,3,4,5});
            s_vector<T,5> v2({-1,-2,-3,-4,-5});
            BOOST_CHECK_EQUAL(-v1,v2);
        }

        BOOST_AUTO_TEST_CASE_TEMPLATE(test_vector_self_add,T,test_types)
        {
            s_vector<T,5> v1({1,2,3,4,5});
            s_vector<T,5> v2({-1,-2,-3,-4,-5});
            s_vector<T,5>&v3=v1+=v2;
            s_vector<T,5> v4;
            BOOST_CHECK_EQUAL(std::addressof(v1),std::addressof(v3));
            BOOST_CHECK_EQUAL(v1,v4);
        }
        BOOST_AUTO_TEST_CASE_TEMPLATE(test_vector_self_diff,T,test_types)
        {
            s_vector<T,5> v1({1,2,3,4,5});
            s_vector<T,5> v2=v1;
            s_vector<T,5>&v3=v1-=v2;
            s_vector<T,5> v4;
            BOOST_CHECK_EQUAL(std::addressof(v1),std::addressof(v3));
            BOOST_CHECK_EQUAL(v1,v4);
        }
    BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test_matrix_real)

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_matrix_create,T,test_types){
        s_matrix<T,2,2> m1({{1,2},{3,4}});
        s_matrix<T,2,2> m2({{1,2},{3,4}});

        BOOST_CHECK_EQUAL(m1[0][0],T(1));
        BOOST_CHECK_EQUAL(m1[0][1],T(2));
        BOOST_CHECK_EQUAL(m1[1][0],T(3));
        BOOST_CHECK_EQUAL(m1[1][1],T(4));
    }
    BOOST_AUTO_TEST_CASE_TEMPLATE(test_matrix_add,T,test_types) {
        s_matrix<T,2,2> M1({{1,2},{3,4}});
        s_matrix<T,2,2> M2({{3,7},{-1,-2}});
        s_matrix<T,2,2> M3({{4,9},{2,2}});
        BOOST_CHECK_EQUAL(M1+M2, M3);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_matrix_diff,T,test_types) {
        s_matrix<T,2,2> M1({{1,2},{3,4}});
        s_matrix<T,2,2> M2({{3,7},{-1,-2}});
        s_matrix<T,2,2> M3({{-2,-5},{4,6}});
        BOOST_CHECK_EQUAL(M1-M2, M3);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_matrix_scalar,T,test_types)
    {
        s_matrix<T,2,2> M1({{1,2},{3,4}});        s_matrix<int,2,2> m1({{1,2},{3,4}});

        s_matrix<T,2,2> M2({{3,6},{9,12}});
        BOOST_CHECK_EQUAL(T(3)*M1,M2);
    }
    BOOST_AUTO_TEST_CASE_TEMPLATE(test_matrix_negate,T,test_types)
    {
        s_matrix<T,2,2> M1({{1,2},{3,4}});
        s_matrix<T,2,2> M2({{-1,-2},{-3,-4}});
        BOOST_CHECK_EQUAL(-M1,M2);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_matrix_mult,T,test_types)
    {
        s_matrix<T,2,2> M1({{1,2},{3,4}});
        s_matrix<T,2,2> M2({{1,2},{3,4}});
        s_matrix<T,2,2> M3({{7,10},{15,22}});
        BOOST_CHECK_EQUAL(M1*M2,M3);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_matrix_self_mult,T,test_types)
    {
        s_matrix<T,2,2> M1({{1,2},{3,4}});
        s_matrix<T,2,2> M2({{1,2},{3,4}});
        s_matrix<T,2,2> &M3=M1*=M2;
        BOOST_CHECK_EQUAL(std::addressof(M3),std::addressof(M1));
        BOOST_CHECK_EQUAL(M1,M2*M2);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_matrix_pow,T,real_types)
    {

        s_matrix<T,4,4> A({{1,2,3,4},{3,1,2,4},{1,4,3,1},{5,3,1,2}});
        s_matrix<T,4,4> B({{T(271373769), T(257707107), T(229627054), T(291157093)},
                           {T(274537505), T(260711810), T(232288922), T(294525263)},
                           {T(236529470), T(224618247), T(200110798), T(253716959)},
                           {T(299967605), T(284862168), T(253759121), T(321726562)}});
        bool approx_equal=true;
        double err=1e-4;

        BOOST_CHECK(approx(pow(A,9),B,err));
    }

    BOOST_AUTO_TEST_SUITE(test_matrix_linalg)
        BOOST_AUTO_TEST_CASE_TEMPLATE(test_matrix_vector,T,test_types)
        {
            s_matrix<T,3,4> A({{1,2,3,4},{3,1,2,4},{1,4,3,1}});
            s_vector<T,4> u({2,1,3,1});
            s_vector<T,3> v({17,17,16});
            BOOST_CHECK_EQUAL(A*u,v);
        }
        BOOST_AUTO_TEST_CASE_TEMPLATE(test_matrix_det,T,real_types)
        {
            s_matrix<T,4,4> A({{1,2,3,4},{3,1,2,4},{1,4,3,1},{5,3,1,2}});
            BOOST_CHECK_CLOSE(A.det(),35,err);
        }

        BOOST_AUTO_TEST_CASE_TEMPLATE(test_matrix_trace,T,real_types)
        {
            s_matrix<T,4,4> A({{1,2,3,4},{3,1,2,4},{1,4,3,1},{5,3,1,2}});
            BOOST_CHECK_CLOSE(A.tr(),7,err);
        }

        BOOST_AUTO_TEST_CASE_TEMPLATE(test_matrix_inv,T,real_types)
        {
            s_matrix<T,4,4> A({{1,2,3,4},{3,1,2,4},{1,4,3,1},{5,3,1,2}});
            BOOST_CHECK(std::abs(l2_distance(A * A.inv(), decltype(A)((T)1)))<err);
        }

        BOOST_AUTO_TEST_CASE_TEMPLATE(test_matrix_inv_rand,T,real_types)
        {
            std::random_device dev;
            std::mt19937_64 rng(dev());
            constexpr int dimension=100;
            constexpr double std=.175;
            std::normal_distribution<double> d(0,std);
            for(int i=0;i<N;i++) {
                s_matrix<T, dimension, dimension> A;
                for (auto &R: A)
                    for (auto &c: R)
                        c = d(rng);
                if(A.det()!=0)
                    BOOST_CHECK(std::abs(l2_distance(A * A.inv(), decltype(A)((T)1)))<err);
            }
        }
        BOOST_AUTO_TEST_SUITE(test_solve)
            BOOST_AUTO_TEST_CASE(test_solve)
            {
                std::random_device dev;
                std::mt19937_64 rng(dev());
                constexpr int dimension=100;
                using T=real;
                L2_inner_product<T,s_vector<T,100>> innerProduct;
                constexpr double std=.175;
                std::normal_distribution<double> d(0,std);
                for(int i=0;i<N;i++) {
                    s_matrix<T, dimension, dimension> A;
                    s_vector<T,dimension> v;
                    for (auto &R: A)
                        for (auto &c: R)
                            c = d(rng);
                    for(auto &c:v)
                        c=d(rng);
                    BOOST_CHECK(innerProduct.metric(A*A.solve(v), v)<err);
                }
            }
            BOOST_AUTO_TEST_CASE(test_solve_singular)
            {
                std::random_device dev;
                std::mt19937_64 rng(dev());
                constexpr int dimension=100;
                using T=real;
                L2_inner_product<T,s_vector<T,100>> innerProduct;
                constexpr double std=.175;
                std::normal_distribution<double> d(0,std);
                for(int i=0;i<N;i++) {
                    s_matrix<T, dimension, dimension> A;
                    s_vector<T,dimension> v;
                    for (auto &R: A)
                        for (auto &c: R)
                            c = d(rng);
                    for(auto &c:v)
                        c=d(rng);
                    BOOST_CHECK(innerProduct.metric(A*A.solve(v), v)<err);
                }
            }
            BOOST_AUTO_TEST_CASE(test_solve_rand)
            {
                std::random_device dev;
                std::mt19937_64 rng(dev());
                constexpr int dimension=100;
                using T=real;
                L2_inner_product<T,s_vector<T,100>> innerProduct;
                constexpr double std=.175;
                std::normal_distribution<double> d(0,std);
                for(int i=0;i<N;i++) {
                    s_matrix<T, dimension, dimension> A;
                    s_vector<T,dimension> v;
                    for (auto &R: A)
                        for (auto &c: R)
                            c = d(rng);
                    for(auto &c:v)
                        c=d(rng);
                    BOOST_CHECK(innerProduct.metric(A*A.solve(v), v)<err);
                }
            }
            BOOST_AUTO_TEST_CASE(test_solve_rand_rank50)
            {
                int success=0;
                for(int i=0;i<N;i++) {
                    bool test= test_solver<50>();
                    success+=test;
                    BOOST_WARN(test);
                }
                BOOST_CHECK(static_cast<double>(success)/N>1-failure_tolerance);
            }

            BOOST_AUTO_TEST_CASE(test_solve_rand_rank100)
            {
                int success=0;
                for(int i=0;i<N;i++) {
                    bool test= test_solver<100>();
                    success+=test;
                    BOOST_WARN(test);
                }
                BOOST_CHECK(static_cast<double>(success)/N>1-failure_tolerance);
            }

            BOOST_AUTO_TEST_CASE(test_solve_rand_rank90)
            {
                int success=0;
                for(int i=0;i<N;i++) {
                    bool test= test_solver<90>();
                    success+=test;
                    BOOST_WARN(test);
                }
                BOOST_CHECK(static_cast<double>(success)/N>1-failure_tolerance);
            }
            BOOST_AUTO_TEST_CASE(test_solve_rand_rank0)
            {
                int success=0;
                for(int i=0;i<N;i++) {
                    bool test= test_solver<0>();
                    success+=test;
                    BOOST_WARN(test);
                }
                BOOST_CHECK(static_cast<double>(success)/N>1-failure_tolerance);
            }
            BOOST_AUTO_TEST_CASE(test_solve_rand_rank99)
            {
                int success=0;
                for(int i=0;i<N;i++) {
                    bool test= test_solver<99>();
                    success+=test;
                    BOOST_WARN(test);
                }
                BOOST_CHECK(static_cast<double>(success)/N>1-failure_tolerance);
            }
            BOOST_AUTO_TEST_CASE(test_rank_rand_rank99)
            {
                int success=0;
                for(int i=0;i<N;i++) {
                    int rank=(random_matrix<100,99>()* random_matrix<99,100>()).rank();
                    success+=rank==99;
                    BOOST_WARN_EQUAL(rank,99);
                }
                BOOST_CHECK(static_cast<double>(success)/N>1-failure_tolerance);
            }
        BOOST_AUTO_TEST_SUITE_END()
    BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()


