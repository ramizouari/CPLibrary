//
// Created by ramizouari on 30/11/2021.
//
#ifndef __LINEAR__ALGEBRA__
#define __LINEAR__ALGEBRA__
#include <vector>
#include <array>
#include "algebra/abstract_algebra.h"
#include "polynomial/polynomial.h"
#include "vector.h"


namespace cp::linalg
{
    struct m_shape
    {
        int n,m;
    };

/**
 * @brief Matrix:
* @details This is the union of R^(n*m) for all n and m
* @Requirements
* <strong>R</strong> is a commutative ring.
* @formal it is the set of matrices over the commutative ring <strong>R</strong> <br>
 * In fact, It is the union of L(R^a,R^b) over all a and b where L(E,F) is the set of matrices acting on E with values over F
*/
    template<ring R>
    struct matrix : tensor_view<R,2>
    {
        std::vector<vector<R>> M;

        using base_field=R;
        using base_ring=R;

        const auto& data() const
        {
            return M;
        }
        auto &data()
        {
            return M;
        }

        matrix() = default;


        matrix(std::size_t rows,std::size_t cols,size_tag_t tag) : M(rows,vector<R>(cols,tag)) {}


        template<std::convertible_to<R> O>
        matrix(const O &k):matrix(1,1,size_tag)
        {
            M[0][0] = k;
        }
        matrix(std::vector<vector<R>> &&_M):M(std::move(_M)){}
        matrix(const std::vector<std::vector<R>> &_M) :M(_M.begin(),_M.end())
        {
        }

        static matrix eye(std::size_t n) {
            matrix I(n,n,size_tag);
            for (int i=0;i<n;i++)
                I[i][i]=1;
            return I;
        }

        std::array<std::size_t,2> shape() const override {
            return {rows(),cols()};
        }

        std::size_t size() const override
        {
            return rows()*cols();
        }

        R& at(std::array<std::size_t,2> I) override
        {
            return M[I[0]][I[1]];
        }

        const R& at(std::array<std::size_t,2> I) const override
        {
            return M[I[0]][I[1]];
        }

        std::size_t rows() const
        {
            return M.size();
        }

        std::size_t cols() const
        {
            return M.empty()?0:M[0].size();
        }

        vector<R>& operator[](int k)
        {
            return M[k];
        }

        R tr() const
        {
            int m=cols(),n=rows();
            R r{};
            for(int i=0;i<std::min(n,m);i++) r+=M[i][i];
            return r;
        }

        matrix T() const
        {
            int m=cols(),n=rows();
            matrix P(m,n);
            for(int i=0;i<n;i++) for(int j=0;j<m;j++)
                P.M[j][i]=M[i][j];
            return P;
        }

        matrix H() const
        {
            int m = cols(), n = rows();
            matrix P(m,n);
            for (int i = 0; i < n; i++) for (int j = 0; j < m; j++)
                    P.M[j][i] = conj<R,R>(M[i][j]);
            return P;
        }

        const auto& operator[](int k) const
        {
            return M[k];
        }

        matrix &operator+=(const matrix &O)
        {
            if (O.rows() == 1 && O.cols() == 1) for (int i=0;i<std::min(rows(),cols());i++)
                M[i][i]+=O.M[0][0];
            else for(int i=0;i<std::min(rows(),O.rows());i++) for(int j=0;j<std::min(cols(),O.cols());j++)
                M[i][j]+=O.M[i][j];
            return *this;
        }

        matrix &operator-=(const matrix &O)
        {
            if (O.rows() == 1 && O.cols() == 1) for (int i=0;i<std::min(rows(),cols());i++)
                M[i][i]-=O.M[0][0];
            else for(int i=0;i<std::min(rows(),O.rows());i++) for(int j=0;j<std::min(cols(),O.cols());j++)
                M[i][j]-=O.M[i][j];
            return *this;
        }

        matrix &operator-=(const R &k)
        {
            for(int i=0;i<std::min(rows(),cols());i++)
                M[i][i]-=k;
            return *this;
        }

        matrix operator-() const
        {
            auto N=*this;
            for(auto &row:N.M) for(auto &s:row) s=-s;
            return N;
        }

        auto & operator*=(R k)
        {
            for(auto &row:M) for(auto &u:row)
                    u*=k;
            return *this;
        }


        matrix operator*(const matrix &B) const
        {
            if (std::min({rows(),B.rows(),cols(),B.cols()}) == 0) return matrix(rows(),B.cols(),size_tag);
            if (cols() == B.rows()) {
                auto n=rows(),p=cols(),m=B.cols();
                matrix C(n,m,size_tag);
                for (int i=0;i<n;i++) for (int k=0;k<p;k++) for (int j=0;j<m;j++)
                    C.M[i][j] += M[i][k] * B.M[k][j];
                return C;
            }
            if (rows() == 1 && cols() == 1)
                return M[0][0] * B;
            if (B.rows() == 1 && B.cols() == 1)
                return *this * B.M[0][0];
            throw std::runtime_error("matrix::operator*(): matrix does not have same size");
        }

        matrix& operator*=(const matrix &B)
        {
            return *this = *this * B;
        }

        vector<R> operator*(const vector<R> &u) const
        {
            int n=rows(),m=cols();
            vector<R> v(n,size_tag);
            for(int j=0;j<m;j++) for(int i=0;i<n;i++)
                v[i]+=M[i][j]*u[j];
            return v;
        }

        matrix &operator/=(R k)
        {
            for(auto &row:M) for(auto &u:row)
                    u/=k;
            return *this;
        }

        matrix operator/(R k) const
        {
            auto N=*this;
            return N/=k;
        }

        auto& operator/=(const matrix& O)
        {
            return *this *= O.inv();
        }

        auto operator/(const matrix &O) const
        {
            return (*this) * O.inv();
        }

        bool operator==(const matrix &O) const {
            if (O.rows() == 0 && O.cols() == 0) {
                if (rows() != cols())
                    return false;
                for (int i=0;i< rows();i++) if (!is_zero(M[i][i])) return false;
                return true;
            }
            if (O.rows() == 1 && O.cols() == 1) {
                if (O.rows() != cols())
                    return false;
                for (int i=0;i<rows();i++)
                    if (!is_zero(M[i][i]-O[0][0])) return false;
                return true;
            }
            if (rows()!= O.rows() || cols()!= O.cols())
                return false;
            for (int i=0;i<rows();i++) for (int j=0;j<cols();j++)
                if (!is_zero(M[i][j]-O.M[i][j])) return false;
            return true;
        }

        auto begin()
        {
            return M.begin();
        }

        auto begin() const
        {
            return M.cbegin();
        }

        auto end()
        {
            return M.end();
        }

        auto end() const
        {
            return M.cend();
        }

        struct MatSolve
        {
            matrix E,Y;
            std::optional<matrix> X; // Results of the row echelon form
            size_t rank; // Rank of E
            bool dir; // Number of swaps modulo 2 in row_echelon is 0. Denotes that E and the original matrix have the same direction
            std::vector<std::pair<int,int>> mapper; // Mapping denoting the pivot cells

            MatSolve(matrix E_,matrix Y_,size_t rank,bool dir,std::vector<std::pair<int,int>> mapper):
                E(std::move(E_)),Y(std::move(Y_)),rank(rank),mapper(std::move(mapper)),dir(dir),X(std::make_optional<matrix>(E.cols(),Y.cols(),size_tag))
            {
            }

            // Solve the system, in general
            void solve_general()
            {
                int n=E.rows(),l=Y.cols();
                for (int i=rank;i<n;i++) for (int j=0;j<l;j++) if (!is_zero(Y[i][j])) {
                    X= std::nullopt; // No solution
                    return;
                }
                for (int k=mapper.size()-1;k>=0;k--)
                {
                    auto [r,i] = mapper[k];
                    for (int j=0;j<r;j++)
                    {
                        auto w = E[j][i]/E[r][i];
                        E[j][i]=0;
                        Y[j]-=w*Y[r];
                    }
                }
                for (auto [r,i]: mapper)
                    (*X)[i]=Y[r]/E[r][i];
            }

            // Solve the system, assuming that E is square and det E != 0
            void solve_invertible()
            {
                int n=E.rows();
                X=Y;
                for (int i=n-1;i>=0;i--)
                {
                    if (is_zero(E[i][i])) throw std::invalid_argument("Matrix is not invertible");
                    for (int j=0;j<i;j++) {
                        auto w = E[j][i]/E[i][i];
                        E[j][i]=0;
                        Y[j]-=w*Y[i];
                    }
                }
                for (int i=0;i<n;i++)
                    (*X)[i]=Y[i]/E[i][i]; // Intentional
            }
        };

        // Rule for selecting the pivot element: Gauss-Jordan Rule
        inline static std::function<int(const matrix&,int,int)> pivot_rule= [](const matrix& X,int rnk,int col) {
            int r=rnk;
            int n=X.rows();
            while (r<n && is_zero(X[r][col])) r++;
            return r;
        };

        MatSolve row_echelon(matrix Y) const
        {
            size_t rnk=0;
            auto E=*this;
            auto n=rows(),m=cols();
            bool dir=true;
            std::vector<std::pair<int,int>> mapper;
            for (int i=0;i<m && rnk < n;i++)
            {
                int r=pivot_rule(E,rnk,i);
                if (r==n) continue;
                mapper.emplace_back(rnk,i);
                if (r!=rnk)
                {
                    std::swap(E[rnk],E[r]);
                    std::swap(Y[rnk],Y[r]);
                    dir=!dir;
                }
                for (int j=rnk+1;j<n;j++)
                {
                    auto w = E[j][i]/E[rnk][i];
                    E[j][i]=0;
                    for (int k=i+1;k<m;k++) E[j][k]-=w*E[rnk][k];
                    Y[j]-=w*Y[rnk];
                }
                rnk++;
            }
            return {E,Y,rnk,dir,mapper};
        }

        // Solve AX = Y, where A and Y are known
        std::optional<matrix> solve(const matrix& O, bool invertible=false) const {
            invertible = invertible && rows() == cols();
            auto dec=row_echelon(O);
            if (invertible)
                dec.solve_invertible();
            else dec.solve_general();
            return dec.X;
        }

        // Solve Ax = Y, where A and b are known
        // 1. If A is known to be invertible, set
        std::optional<vector<R>> solve(const vector<R> &V , bool invertible =false) const {
            matrix Y(V.size(),1,size_tag);
            for (int i=0;i<V.size();i++) Y[i][0]=V[i];
            auto X=solve(Y,invertible);
            if (!X.has_value()) return std::nullopt; // If no solutions, return the empty vector
            vector<R> U(cols(),size_tag);
            for (int i=0;i<cols();i++) U[i]=(*X)[i][0];
            return U;
        }

        // Basis of vectors (e_1,..,e_r) such that Ae_k=0 for all k
        matrix null_basis() const {
            int n=rows(),m=cols();
            matrix Z(*this);
            for (int i=0;i<m;i++) {
                Z.M.emplace_back(m,size_tag);
                Z[rows()+i][i]=1;
            }
            auto C = Z.T().row_echelon(matrix(m,0,size_tag)).E;
            matrix B;
            for (int i=0;i<m;i++) if (all_of(C[i].begin(),C[i].begin()+n,[](auto x) {return is_zero(x);}))
            {
                vector<R> u(m,size_tag);
                for (int j=0;j<m;j++) u[j] = C[i][n+j];
                B.M.push_back(u);
            }
            return B;
        }

        // Basis induced by the matrix
        matrix image_basis() const {
            auto [E,_1,_2,rnk,_3] = row_echelon(*this,false);
            E.resize(rnk);
            return E;
        }

        // Inverse of a square matrix.
        matrix inv() const
        {
            return *solve(matrix::eye(rows()),true);
        }

        // Determinant of a square matrix
        R det() const
        {
            auto dec = row_echelon(matrix(rows(),0,size_tag));
            R w=dec.dir?1:-1;
            for (int i=0;i<rows();i++) w*=dec.E[i][i];
            return w;
        }

        // Rank of a matrix:
        // 1. The number of independent columns/rows in that matrix
        // 2. The dimension of vector space induced by the matrix
        // 3. Size of the basis generated by the matrix
        size_t rank() const {
            return row_echelon(matrix(rows(),0,size_tag)).rank;
        }

        // Nullity of the matrix:
        // 1. Size of the basis of elements that annihilate the matrix: Ax = 0
        // 2. S
        size_t nullity() const {
            return cols() - rank();
        }

    };

    template<typename Mat,typename R>
    concept ToMatrix=std::convertible_to<Mat,matrix<R>>;

    template<typename Mat,typename R>
    concept ToMatrixProper=std::convertible_to<Mat,matrix<R>> && !std::same_as<Mat,matrix<R>>;


    template<ring R ,ToMatrix<R> O>
    matrix<R> operator+(const matrix<R> &A,const O &B)
    {
        auto C=A;
        return C+=B;
    }

    template<ring R ,ToMatrixProper<R> O>
    matrix<R> operator+(const O & A,const matrix<R> & B)
    {
        matrix<R> C=A;
        return C+=B;
    }

    template<ring R ,ToMatrix<R> O>
    matrix<R> operator-(const matrix<R> &A,const O &B)
    {
        auto C=A;
        return C-=B;
    }

    template<ring R ,ToMatrixProper<R> O>
    matrix<R> operator-(const O & A,const matrix<R> & B)
    {
        matrix<R> C=A;
        return C-=B;
    }

    template<ring R ,std::convertible_to<R> O>
    matrix<R> operator*(const O &k,const matrix<R> &A)
    {
        auto C=A;
        return C*=k;
    }

    template<ring R ,std::convertible_to<R> O>
    matrix<R> operator*(const matrix<R> &A,const O &k)
    {
        auto C=A;
        return C*=k;
    }


    template<ring R ,std::convertible_to<R> O>
    matrix<R> operator/(const matrix<R> &A,const O &k)
    {
        auto C=A;
        return C/=k;
    }

/**
 * @details Static Matrix
 * It is an element of the vector space L(R^n,R^m)
 * @tparam R the base commutative ring
 * @tparam n the number of rows
 * @tparam m the number of columns
 * @Requirements R is a commutative ring
 * @Requirements n>=0
 * @Requirements m>=0
 * @Formal
 * <ul>
 * <li> Multiplication between 2 matrices is defined if the shapes are compatible
 * <li> Multiplication between a matrix and a vector is defined if the shapes are compatible
 * <li> It is an associative algebra if n=m
 * </ul>
 * */

    template<typename R,int n,int m>
    class s_matrix
    {
        std::array<std::array<R,m>,n> M;
    public:
        using base_field=R;
        using base_ring=R;
        bool operator==(const s_matrix&) const = default;
        s_matrix(R k=0)
        {
            for(int i=0;i<n;i++) for(int j=0;j<m;j++)
                    M[i][j]=i==j?k:0;
        }
        s_matrix(std::array<std::array<R,m>,n> _M):M(std::move(_M)){}

        inline static constexpr int row_dim()
        {
            return n;
        }

        s_matrix(const std::vector<std::array<R,m>> &_M)
        {
            int counter=0;
            for(int i=0;i<n;i++) for(int j=0;j<m;j++)
                    M[i][j]=_M[i][j];
        }

        inline static constexpr int col_dim()
        {
            return m;
        };

        auto& operator[](int k)
        {
            return M[k];
        }

        const auto& operator[](int k) const
        {
            return M[k];
        }

        R tr() const
        {
            R r=0;
            for(int i=0;i<std::min(n,m);i++)
                r+=M[i][i];
            return r;
        }

        s_matrix<R,m,n> T() const
        {
            s_matrix<R,m,n> P;
            for(int i=0;i<n;i++) for(int j=0;j<m;j++)
                    P.M[j][i]=M[i][j];
            return P;
        }

        s_matrix<R, m, n> H() const
        {
            s_matrix<R, m, n> P;
            for (int i = 0; i < n; i++) for (int j = 0; j < m; j++)
                    P.M[j][i] = conj(M[i][j]);
            return P;
        }

        auto &operator+=(const s_matrix &O)
        {
            for(int i=0;i<n;i++) for(int j=0;j<m;j++)
                    M[i][j]+=O.M[i][j];
            return *this;
        }

        auto &operator-=(const s_matrix &O)
        {
            for(int i=0;i<n;i++) for(int j=0;j<m;j++)
                    M[i][j]-=O.M[i][j];
            return *this;
        }

        auto operator+(const s_matrix &O) const
        {
            auto N=*this;
            return N+=O;
        }

        auto operator-(const s_matrix &O) const
        {
            auto N=*this;
            return N-=O;
        }

        auto operator-() const
        {
            auto N=*this;
            for(auto &row:N.M) for(auto &s:row)
                    s=-s;
            return N;
        }

        template<int p>
        s_matrix<R,n,p> operator*(const s_matrix<R,m,p> &O) const
        {
            s_matrix<R,n,p> N(0);
            for(int i=0;i<n;i++) for(int k=0;k<m;k++) for(int j=0;j<p;j++)
                        N[i][j]+=M[i][k]*O[k][j];
            return N;
        }

        auto &operator*=(const s_matrix &O)
        {
            static_assert(n==m);
            auto N=(*this)*O;
            M.swap(N.M);
            return *this;
        }

        auto & operator*=(R k) {
            for (auto &row: M)
                for (auto &u: row)
                    u *= k;
            return *this;
        }
        s_vector<R,n> operator*(const s_vector<R,m> &u) const
        {
            s_vector<R,n> v;
            for(int j=0;j<m;j++) for(int i=0;i<n;i++)
                    v[i]+=M[i][j]*u[j];
            return v;
        }

        auto &operator/=(R k)
        {
            for(auto &row:M) for(auto &u:row)
                    u/=k;
            return *this;
        }

        auto operator/(R k) const
        {
            auto N=*this;
            return N/=k;
        }

        auto& operator/=(const s_matrix &O)
        {
            return *this*=O.inv();
        }

        auto operator/(const s_matrix &O) const
        {
            return (*this) * O.inv();
        }

        auto begin()
        {
            return M.begin();
        }

        auto begin() const
        {
            return M.cbegin();
        }

        auto end()
        {
            return M.end();
        }

        auto end() const
        {
            return M.cend();
        }

        auto row_echelon_form() const
        {
            auto P=*this;
            int s=0;
            for(int i=0;i<n;i++)
            {
                int p=s;
                while(p<n && is_zero(P.M[p][i]))
                    p++;
                if(p==n)
                    continue;
                P.M[p].swap(P.M[s]);
                R w=P.M[s][i];
                for(int j=s+1;j<n;j++)
                {
                    R r=P.M[j][i]/w;
                    for (int k = i; k < m; k++)
                        P.M[j][k]-=r*P.M[i][k];
                }
                s++;
            }
            return P;
        }

        int rank() const
        {
            auto E=row_echelon_form();
            int r=0;
            for(int i=0,j=0;i<n&&j<m;j++)
                if(! is_zero(E.M[i][j]))
                {
                    r++;
                    i++;
                }
            return r;
        }

        int nullity() const
        {
            return row_dim()-rank();
        }

        R det() const
        {
            static_assert(n==m);
            auto P=*this;
            bool invert=false;
            for(int i=0;i<n;i++)
            {
                int p=i;
                while(p<n && is_zero(P.M[p][i]))
                    p++;
                if(p==n)
                    continue;
                if(p!=i)
                {
                    std::swap(P.M[p], P.M[i]);
                    invert=!invert;
                }
                R w=P.M[i][i];
                for(int j=i+1;j<n;j++)
                {
                    R r=P.M[j][i]/w;
                    for (int k = i; k < m; k++)
                        P.M[j][k]-=r*P.M[i][k];
                }
            }
            R d=1;
            for(int i=0;i<n;i++)
                d*=P.M[i][i];
            return invert?-d:d;
        }

        s_matrix inv() const
        {
            static_assert(n==m);
            s_matrix P=*this,Q(1);
            for(int i=0;i<n;i++)
            {
                int p=i;
                while(p<n && is_zero(P.M[p][i]))
                    p++;
                if(p==n)
                    continue;
                std::swap(P.M[p], P.M[i]);
                std::swap(Q.M[p],Q.M[i]);
                R w=P.M[i][i];
                for(int j=i+1;j<n;j++)
                {
                    R r=P.M[j][i]/w;
                    for (int k = 0; k < m; k++)
                    {
                        P.M[j][k] -= r*P.M[i][k];
                        Q.M[j][k] -= r*Q.M[i][k];
                    }
                }
            }
            for(int i=n-1;i>=0;i--)
            {
                R w=P.M[i][i];
                for(int j=0;j<n;j++)
                    Q.M[i][j]/=w;
                for(int k=i-1;k>=0;k--)
                {
                    R r=P.M[k][i];
                    for (int j = 0; j < n; j++)
                        Q.M[k][j] -= r*Q.M[i][j];
                }
            }
            return Q;
        }

        s_vector<R,m> solve(s_vector<R,n> A) const
        {
            static_assert(n==m);
            s_matrix P=*this;
            for(int i=0;i<n;i++)
            {
                int p=i;
                while(p<n && P.M[p][i]==R(0))
                    p++;
                if(p==n)
                    continue;
                std::swap(P.M[p], P.M[i]);
                std::swap(A[p],A[i]);
                R w=P.M[i][i];
                for(int j=i+1;j<n;j++)
                {
                    if(is_zero(w))
                        continue;
                    R r=P.M[j][i]/w;
                    for (int k = 0; k < m; k++)
                        P.M[j][k] -= r*P.M[i][k];
                    A[j]-=r*A[i];
                }
            }
            for(int i=n-1;i>=0;i--)
            {
                R w=P.M[i][i];
                if(w==R(0))
                    continue;
                A[i]/=w;
                for(int k=i-1;k>=0;k--)
                {
                    R r=P.M[k][i];
                    A[k] -= r*A[i];
                }
            }
            return A;
        }

    };

    template<typename R,int n,int m>
    s_matrix<R,n,m> operator*(const R&a,const s_matrix<R,n,m> &M)
    {
        auto N=M;
        return N*=a;
    }
}
#endif // __LINEAR__ALGEBRA__