#include <iostream>
#include <chrono>
#include <numeric>
#include <complex>
#include <array>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <map>
#include "linear_algebra/tensor.h"
#include "signals/fft.h"
#include "signals/multi_fft.h"
#include "signals/hadamard.h"

unsigned bit_ceil(unsigned x)
{
    unsigned r=1;
    while(r<x)
        r<<=1;
    return r;
}

enum class FFTNormalization
{
    None,
    Sqrt,
    Normalized
};
using real=double;
using integer=long long;
using complex=std::complex<real>;
real constexpr pi=3.14159265359;

using Real=real;
using Complex=std::complex<Real>;


struct ComplexArrayView : public cp::linalg::tensor_view<Complex,1>
{
    ComplexArrayView(Complex* src,int n):src(src),n(n){}
    ComplexArrayView(std::vector<Complex>& src):src(src.data()),n(src.size()){}
    ComplexArrayView(std::vector<Complex>&& src):src(src.data()),n(src.size()){}
    template<size_t m>
    ComplexArrayView(std::array<Complex,m>& src):src(src.data()),n(src.size()){}
    template<size_t m>
    ComplexArrayView(std::array<Complex,m>&& src):src(src.data()),n(src.size()){}
    Complex& at(int i)
    {
        return src[i];
    }
    const Complex& at(int i) const
    {
        return src[i];
    }
    Complex& at(std::array<std::size_t,1> indexes) override
    {
        return src[indexes[0]];
    }

    const Complex& at(std::array<std::size_t,1> indexes) const override
    {
        return src[indexes[0]];
    }

    Complex* src;
    std::size_t n;
    std::array<std::size_t,1> shape() const override
    {
        return {n};
    }

    ComplexArrayView slice(int i,int j) const
    {
        return ComplexArrayView(src+i,j-i);
    }
};

template<typename R>
struct VectorSubView : public cp::linalg::tensor_view<R,1>
{
    cp::linalg::tensor_view<R,1>* src;
    int i;
    std::size_t n;
    VectorSubView(cp::linalg::tensor_view<R,1>& src,int i,int n):i(i),n(n)
    {
        auto S=dynamic_cast<VectorSubView<R>*>(&src);
        if(S)
        {
            this->src=S->src;
            this->i+=S->i;
        }
        else
            this->src=&src;
    }
    R& at(int j)
    {
        return src->at(i+j);
    }
    const R& at(int j) const
    {
        return src->at(i+j);
    }
    R& at(std::array<std::size_t,1> indexes) override
    {
        return src->at(i+indexes[0]);
    }

    const R& at(std::array<std::size_t,1> indexes) const override
    {
        return src->at(i+indexes[0]);
    }
    std::array<std::size_t,1> shape() const override
    {
        return {n};
    }
};


void inplace_fft2(cp::linalg::tensor_view<Complex,1> & a, bool inverse, FFTNormalization normalized = FFTNormalization::Sqrt)
{
    int n = a.size();

    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            swap(a(i), a(j));
    }

    for (int len = 2; len <= n; len <<= 1) {
        Real ang = 2 * pi / len * (inverse ? -1 : 1);
        Complex wlen = std::polar<Real>(1.,ang);
        for (int i = 0; i < n; i += len) {
            Complex w(1);
            for (int j = 0; j < len / 2; j++) {
                Complex u = a(i+j), v = a(i+j+len/2) * w;
                a(i+j) = u + v;
                a(i+j+len/2) = u - v;
                w *= wlen;
            }
        }
    }
    Real r;
    switch (normalized)
    {
        case FFTNormalization::None:
            r=1;
            break;
        case FFTNormalization::Sqrt:
            r=std::sqrt(n);
            break;
        case FFTNormalization::Normalized:
            r=n;
            break;
    }
    if(normalized!=FFTNormalization::None) for (Complex & x : a)
            x /= r;
}

void inplace_fft2(cp::linalg::tensor_view<Complex,1> && a, bool inverse, FFTNormalization normalized = FFTNormalization::Sqrt)
{
    inplace_fft2(a,inverse,normalized);
}

integer smallest_prime_divisor(integer n, const std::vector<integer> &primes)
{
    for(auto p:primes)
    {
        if(p*p>n)
            break;
        if(n%p==0)
            return p;
    }
    return n;
}

integer smallest_prime_divisor(integer n)
{
    if(n>10)
        throw std::runtime_error("n too large for this problem");
    static std::vector<integer> primes={2,3,5,7};
    return smallest_prime_divisor(n,primes);
}

int fft_permutation_index(int k,int n,int p)
{
    return (k%p)*(n/p)+(k/p);
}

void permute_elements(cp::linalg::tensor_view<Complex,1> &a,integer p)
{
    std::vector<bool> visited(a.size());
    for(int i=0;i<a.size();i++) if(!visited[i])
    {
        int j=i;
        visited[j]=true;
        do
        {
            j=fft_permutation_index(j,a.size(),p);
            if(j!=i)
                std::swap(a(i),a(j));
            visited[j]=true;
        }while(j!=i);
    }
}

void inplace_fft(cp::linalg::tensor_view<Complex,1> && a, bool inverse, FFTNormalization normalized = FFTNormalization::Sqrt);


void inplace_fft(cp::linalg::tensor_view<Complex,1> & a, bool inverse, FFTNormalization normalized = FFTNormalization::Sqrt)
{
    if(a.size()==1)
        return;
    int n = a.size();
    auto p=smallest_prime_divisor(n);
    permute_elements(a,p);
    for(int i=0;i<p;i++)
        inplace_fft(VectorSubView(a,i*n/p,n/p),inverse,FFTNormalization::None);
    Complex w=std::polar<Real>(1.,2*pi/n*(inverse?-1:1));
    Complex h=std::polar(1.,2*pi/p * (inverse?-1:1));
    Complex t=1;
    for(int i=0;i<n/p;i++)
    {
        std::vector<Complex> U(p);
        Complex u=1;
        Complex v=1;
        Complex z=1;
        for(int j=0;j<p;j++)
        {
            U[j] = a(i + j * n / p) * z;
            z*=t;
        }
        z=1;
        for(int j=0;j<p;j++)
        {
            a(i + j * n / p)=0;
            for (int k = 0; k < p; k++)
            {
                a(i + j * n / p) += U[k] * z;
                z *= u;
            }
            u*=w;
        }
    }
    Real r;
    switch (normalized)
    {
        case FFTNormalization::None:
            r=1;
            break;
        case FFTNormalization::Sqrt:
            r=std::sqrt(n);
            break;
        case FFTNormalization::Normalized:
            r=n;
            break;
    }
    if(normalized!=FFTNormalization::None) for (Complex & x : a)
            x /= r;
}

void inplace_fft(cp::linalg::tensor_view<Complex,1> && a, bool inverse, FFTNormalization normalized)
{
    inplace_fft(a,inverse,normalized);
}

template<typename R>
R mod_operator(R x, R y)
{
    return (x%y+y)%y;
}

using Integer=integer;

//Bluestein's algorithm
std::vector<Complex> general_fft(const cp::linalg::tensor_view<Complex,1> &a, bool inverse, FFTNormalization normalized = FFTNormalization::Sqrt)
{
    unsigned int n = a.size();
    if(bit_ceil(n)==n)
    {
        std::vector<Complex> b(n);
        for(int i=0;i<n;i++)
            b[i]=a(i);
        inplace_fft2(ComplexArrayView(b),inverse);
        return b;
    }
    //m>=2*n
    auto m=bit_ceil(2*n-1);
    std::vector<Complex> w(m);
    Real ang = pi  / n * (inverse ? 1 : -1);
    w[1]=std::polar<Real>(1.,ang);
    w[0]=1;
    for(int i=2;i<m;i++)
        w[i]=w[i-1]*w[1];
    std::vector<Complex> A(m),B(m),W(m);
    for(size_t i=0;i<n;i++)
    {
        auto r=mod_operator<Integer>(i*i,2*n);
        W[i]=w[r];
        A[i]=a(i)*w[2*n-r];
        B[i]=w[r];
    }
    for(size_t i=1;i<n;i++)
        B[m-i]=B[i];
    inplace_fft2(ComplexArrayView(A),false,FFTNormalization::None);
    inplace_fft2(ComplexArrayView(B),false,FFTNormalization::None);
    std::vector<Complex> C(m);
    for(size_t i=0;i<m;i++)
        C[i]=A[i]*B[i];
    inplace_fft2(ComplexArrayView(C),true,FFTNormalization::Normalized);
    for(size_t i=0;i<n;i++)
        C[i]*=std::conj(W[i]);
    C.resize(n);
    Real r;
    switch (normalized)
    {
        case FFTNormalization::None:
            r=1;
            break;
        case FFTNormalization::Sqrt:
            r=std::sqrt(n);
            break;
        case FFTNormalization::Normalized:
            r=n;
            break;
    }
    if(normalized!=FFTNormalization::None) for (Complex & x : C)
            x /= r;
    return C;
}

namespace algebra
{
    using IC=Complex;
}

template<size_t Rank, typename R>
struct FlatView : public cp::linalg::tensor_view<Complex,1>
{
    cp::linalg::flat_tensor<R,Rank>& src;
    FlatView(cp::linalg::flat_tensor<R,Rank>& src):src(src){}
    R& at(std::array<std::size_t,1> indexes)
    {
        return at(indexes[0]);
    }
    R& at(size_t r)
    {
        return src.at(r);
    }

    const R& at(size_t r) const
    {
        return src.at(r);
    }

    const R& at(std::array<std::size_t,1> indexes) const
    {
        return at(indexes[0]);
    }

    std::array<std::size_t,1> shape() const
    {
        return {(std::size_t)src.size()};
    }
};

template<size_t Rank,typename R>
struct TensorProjectionView : public cp::linalg::tensor_view<Complex,1>
{
    size_t k;
    cp::linalg::flat_tensor<R,Rank>& src;
    std::array<std::size_t,Rank> fixed;
    TensorProjectionView(cp::linalg::flat_tensor<R,Rank>& src,size_t k,
                         std::array<std::size_t,Rank> fixed):src(src),k(k),fixed(fixed){}
    R& at(std::array<std::size_t,1> indexes)
    {
        return at(indexes[0]);
    }


    R& at(size_t r)
    {
        auto H=fixed;
        H[k]=r;
        return src.at(H);
    }

    const R& at(size_t r) const
    {
        auto H=fixed;
        H[k]=r;
        return src.at(H);
    }

    const R& at(std::array<std::size_t,1> indexes) const
    {
        return at(indexes[0]);
    }

    std::array<std::size_t,1> shape() const
    {
        return {(std::size_t)src.m_shape[k]};
    }
};

template<typename R,size_t Rank>
R amplitude(const std::array<R,Rank> &X)
{
    return std::accumulate(X.begin(),X.end(),0,std::bit_xor<R>());
}

template<typename R1,typename R2,size_t Rank>
R2 phase(const std::array<R1,Rank> &X, R2 theta)
{
    R2 r{};
    for(int i=Rank-1;i>=0;i--)
        r=X[i]-r;
    return r*theta;
}

template<typename R>
R amplitude(const std::vector<R> &X)
{
    return std::accumulate(X.begin(),X.end(),0,std::bit_xor<R>());
}

template<typename R1,typename R2>
R2 phase(const std::vector<R1> &X, R2 theta)
{
    R2 r{};
    for(int i=X.size()-1;i>=0;i--)
        r=X[i]-r;
    return r*theta;
}

template<typename R,size_t Rank>
std::array<R,Rank> & increment(std::array<R,Rank> &X, const std::array<R,Rank> &shape)
{
    for(int i=0;i<Rank;i++)
    {
        X[i]++;
        if(X[i]<shape[i])
            break;
        X[i]=0;
    }
    return X;
}

template<typename R>
std::vector<R> & increment(std::vector<R> &X, const std::vector<R> &shape)
{
    for(int i=0;i<shape.size();i++)
    {
        X[i]++;
        if(X[i]<shape[i])
            break;
        X[i]=0;
    }
    return X;
}

int main()
{
    int T;
    std::cin >> T;
    if(T) while(T--)
    {
        std::array<std::size_t,5> shape,idx{};
        for(auto &s:shape)
            std::cin >> s;
        Real alpha;
        std::cin >> alpha;
        auto t1=std::chrono::high_resolution_clock::now();
        cp::linalg::flat_tensor<Complex, 5> A(shape);
        do
        {
            A(idx)=std::polar<Real>(amplitude(idx),phase(idx,alpha));
            increment(idx,shape);
        }while(std::any_of(idx.begin(),idx.end(),[&](int i){return i!=0;}));

        std::shared_ptr<cp::signals::bluestein_fft<Complex>> F=std::make_shared<cp::signals::bluestein_fft<Complex>>();
        cp::signals::multi_fft<Complex> M(F);
        M.transform(A,false);
        auto sum=std::accumulate(A.begin(),A.end(),0.,[](Real x,Complex y){return x+std::abs(y.real());});
        auto K=std::accumulate(shape.begin(),shape.end(),1,std::multiplies<>());
        auto t2=std::chrono::high_resolution_clock::now();
        std::cout << std::setprecision(6) << sum/K << std::endl;
        std::cout << std::endl;
        std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << std::endl;
    }
    else
    {
        std::size_t n;
        std::cin >> n;
        cp::linalg::flat_tensor<integer,1> A({n});
        for(auto &a:A)
            a=1;
        cp::signals::faster_hadamard<integer > F;
        auto t1=std::chrono::high_resolution_clock::now();
        F.transform(A,false,cp::signals::FFTNormalization::None);
        auto t2=std::chrono::high_resolution_clock::now();
        for(auto &a:A)
            std::cout << a << '\n';
        std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << std::endl;
    }
}