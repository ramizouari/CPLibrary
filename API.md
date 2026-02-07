# CPLibrary API Documentation

This document provides a comprehensive guide to the CPLibrary API, organized by module. For detailed mathematical background, see the documentation in the `doc/` directory.

## Table of Contents

- [Core Concepts](#core-concepts)
- [Data Structures](#data-structures)
- [Algebra](#algebra)
- [Linear Algebra](#linear-algebra)
- [Number Theory](#number-theory)
- [Polynomials](#polynomials)
- [FFT & Signal Processing](#fft--signal-processing)
- [Graph Algorithms](#graph-algorithms)
- [Combinatorics](#combinatorics)
- [Machine Learning](#machine-learning)
- [Functional Programming](#functional-programming)
- [Topology & Optimization](#topology--optimization)
- [String Algorithms](#string-algorithms)

## Core Concepts

### Namespace

All library functions and types are in the `cp` namespace:

```cpp
namespace cp {
    // All library code here
}
```

### Type Aliases

Common type aliases used throughout the library:

```cpp
using integer = long long;
using real = long double;
```

### Include Pattern

Include only what you need:

```cpp
#include "algebra/binary_operation.h"
#include "data_structures/segment_tree.h"
#include "linear_algebra/matrix.h"
// etc.
```

## Data Structures

### Segment Tree

**Header**: `data_structures/segment_tree.h`, `data_structures/range_queries.h`

Segment tree supporting range queries over associative binary operations.

#### Fixed-Size Segment Tree

```cpp
#include "data_structures/fixed/segment_tree.h"

// Create segment tree with size N and binary operation
template<typename T, typename BinaryOp>
class segment_tree {
public:
    segment_tree(size_t n, BinaryOp op);
    segment_tree(const std::vector<T>& data, BinaryOp op);

    // Query range [l, r)
    T query(size_t l, size_t r) const;

    // Update element at index i
    void update(size_t i, const T& value);

    // Get element at index i
    T operator[](size_t i) const;
};
```

#### Dynamic Segment Tree

```cpp
#include "data_structures/dynamic/segment_tree.h"

// Similar interface to fixed segment tree
// but supports dynamic resizing
```

#### Usage Example

```cpp
#include "data_structures/range_queries.h"
#include <vector>
#include <algorithm>

// Range sum query
std::vector<int> data = {1, 3, 5, 7, 9};
auto seg_tree = cp::make_segment_tree(data, std::plus<>());

int sum = seg_tree.query(1, 4);  // Sum of elements [1, 4)

// Range minimum query
auto min_tree = cp::make_segment_tree(data,
    [](int a, int b) { return std::min(a, b); });

int min_val = min_tree.query(0, 5);
```

### Fenwick Tree (Binary Indexed Tree)

**Header**: `data_structures/fenwick_tree.h`

Fenwick tree for range queries over group operations (requires inverse).

```cpp
#include "data_structures/fixed/fenwick_tree.h"

template<typename T, typename Group>
class fenwick_tree {
public:
    fenwick_tree(size_t n, Group group_op);

    // Query prefix sum [0, i)
    T prefix_query(size_t i) const;

    // Query range [l, r)
    T query(size_t l, size_t r) const;

    // Add value to element at index i
    void add(size_t i, const T& delta);

    // Update element at index i
    void update(size_t i, const T& value);
};
```

#### Usage Example

```cpp
#include "data_structures/fenwick_tree.h"

// Fenwick tree for sum queries
std::vector<int> data = {1, 2, 3, 4, 5};
cp::fenwick_tree<int, std::plus<>> ft(data.size(), std::plus<>());

for (size_t i = 0; i < data.size(); i++) {
    ft.add(i, data[i]);
}

int sum = ft.query(1, 4);  // Sum of elements [1, 4)
ft.add(2, 10);  // Add 10 to element at index 2
```

### Sparse Table

**Header**: `data_structures/fixed/sparse_array.h`

Sparse table for range queries over associative idempotent operations.

```cpp
template<typename T, typename BinaryOp>
class sparse_table {
public:
    sparse_table(const std::vector<T>& data, BinaryOp op);

    // Query range [l, r]
    T query(size_t l, size_t r) const;
};
```

#### Usage Example

```cpp
#include "data_structures/fixed/sparse_array.h"

// Range minimum query (static)
std::vector<int> data = {5, 2, 8, 1, 9, 3};
cp::sparse_table<int, decltype([](int a, int b){ return std::min(a, b); })>
    st(data, [](int a, int b){ return std::min(a, b); });

int min_val = st.query(1, 4);  // Minimum in range [1, 4]
```

### Order Statistic Tree

**Header**: `data_structures/statistic_tree.h`

Self-balancing tree with order statistics.

```cpp
template<typename Key, typename Value, typename Compare = std::less<Key>>
class order_statistic_tree {
public:
    // Insert key-value pair
    void insert(const Key& key, const Value& value);

    // Find k-th smallest element (0-indexed)
    std::pair<Key, Value> select(size_t k) const;

    // Count elements less than key
    size_t rank(const Key& key) const;

    // Erase key
    void erase(const Key& key);

    // Find value by key
    Value* find(const Key& key);

    size_t size() const;
};
```

### B-Tree

**Header**: `data_structures/b_tree.h`

B-tree implementation with guaranteed O(log n) performance.

```cpp
template<typename Key, typename Value, size_t M = 64>
class b_tree {
public:
    void insert(const Key& key, const Value& value);
    Value* find(const Key& key);
    void erase(const Key& key);
    size_t size() const;
};
```

## Algebra

### Binary Operations

**Header**: `algebra/binary_operation.h`

Generic binary operations for algebraic structures.

```cpp
namespace cp {
    // Common binary operations
    template<typename T>
    struct plus_op {
        T operator()(const T& a, const T& b) const { return a + b; }
        T identity() const { return T(0); }
    };

    template<typename T>
    struct multiplies_op {
        T operator()(const T& a, const T& b) const { return a * b; }
        T identity() const { return T(1); }
    };

    template<typename T>
    struct min_op {
        T operator()(const T& a, const T& b) const {
            return std::min(a, b);
        }
    };

    template<typename T>
    struct max_op {
        T operator()(const T& a, const T& b) const {
            return std::max(a, b);
        }
    };
}
```

### Fast Exponentiation

**Header**: `algebra/abstract_algebra.h`

Fast exponentiation over monoids.

```cpp
template<typename T, typename BinaryOp>
T pow(T base, integer exponent, BinaryOp op);

// For multiplicative monoids
template<typename T>
T pow(T base, integer exponent);
```

#### Usage Example

```cpp
#include "algebra/abstract_algebra.h"

// Fast exponentiation
int result = cp::pow(2, 10);  // 2^10 = 1024

// Matrix exponentiation
cp::d_matrix<int> M = {{1, 1}, {1, 0}};
auto M_n = cp::pow(M, 100);  // Fibonacci matrix power
```

### Extended Euclidean Algorithm

```cpp
// Extended GCD over integral domains
template<typename T>
struct egcd_result {
    T gcd;
    T x, y;  // Bézout coefficients: gcd = x*a + y*b
};

template<typename T>
egcd_result<T> extended_gcd(const T& a, const T& b);
```

### Ring Extensions

**Header**: `rings/quadratic/dynamic.h`, `rings/quadratic/static.h`

Quadratic extensions of rings (e.g., ℤ[√2], ℚ(√-1)).

```cpp
template<typename Ring, typename D>
class quadratic_extension {
public:
    quadratic_extension(const Ring& a, const Ring& b); // a + b√D

    // Arithmetic operations
    quadratic_extension operator+(const quadratic_extension& other) const;
    quadratic_extension operator*(const quadratic_extension& other) const;

    // Conjugate: a + b√D → a - b√D
    quadratic_extension conjugate() const;

    // Inverse (if base ring is a field)
    quadratic_extension inverse() const;
};
```

## Linear Algebra

### Vectors

**Header**: `linear_algebra/vector.h`

#### Dynamic Vector

```cpp
template<typename T>
class d_vector {
public:
    // Constructors
    d_vector();
    d_vector(size_t n);
    d_vector(size_t n, const T& value);
    d_vector(std::initializer_list<T> init);
    d_vector(const std::vector<T>& vec);

    // Access
    T& operator[](size_t i);
    const T& operator[](size_t i) const;

    size_t size() const;

    // Arithmetic operations
    d_vector operator+(const d_vector& other) const;
    d_vector operator-(const d_vector& other) const;
    T operator*(const d_vector& other) const;  // Dot product
    d_vector operator*(const T& scalar) const;

    // Norm
    T norm() const;
    T squared_norm() const;
};
```

#### Static Vector

```cpp
template<typename T, size_t N>
class s_vector {
public:
    // Similar interface to d_vector
    // but with compile-time size
    static constexpr size_t dimension = N;
};
```

#### Usage Example

```cpp
#include "linear_algebra/vector.h"

// Dynamic vector
cp::d_vector<double> v1 = {1.0, 2.0, 3.0};
cp::d_vector<double> v2 = {4.0, 5.0, 6.0};

auto v3 = v1 + v2;           // Vector addition
double dot = v1 * v2;        // Dot product
double norm = v1.norm();     // Euclidean norm

// Static vector (compile-time size)
cp::s_vector<int, 3> sv = {1, 2, 3};
```

### Matrices

**Header**: `linear_algebra/matrix.h`

#### Dynamic Matrix

```cpp
template<typename T>
class d_matrix {
public:
    // Constructors
    d_matrix();
    d_matrix(size_t rows, size_t cols);
    d_matrix(size_t rows, size_t cols, const T& value);
    d_matrix(std::initializer_list<std::initializer_list<T>> init);
    d_matrix(const std::vector<std::vector<T>>& mat);

    // Access
    T& operator()(size_t i, size_t j);
    const T& operator()(size_t i, size_t j) const;

    size_t rows() const;
    size_t cols() const;

    // Arithmetic operations
    d_matrix operator+(const d_matrix& other) const;
    d_matrix operator-(const d_matrix& other) const;
    d_matrix operator*(const d_matrix& other) const;  // Matrix multiplication
    d_vector<T> operator*(const d_vector<T>& vec) const;  // Matrix-vector
    d_matrix operator*(const T& scalar) const;

    // Transpose
    d_matrix transpose() const;

    // Identity matrix
    static d_matrix identity(size_t n);
};
```

#### Matrix Operations

```cpp
// Determinant
template<typename T>
T determinant(const d_matrix<T>& mat);

// Matrix inverse
template<typename T>
d_matrix<T> inverse(const d_matrix<T>& mat);

// Solve linear system Ax = b
template<typename T>
d_vector<T> solve(const d_matrix<T>& A, const d_vector<T>& b);

// Characteristic polynomial
template<typename T>
polynomial<T> characteristic_polynomial(const d_matrix<T>& mat);

// Rank
template<typename T>
size_t rank(const d_matrix<T>& mat);
```

#### Usage Example

```cpp
#include "linear_algebra/matrix.h"

// Create matrix
cp::d_matrix<double> A = {
    {1, 2, 3},
    {4, 5, 6},
    {7, 8, 9}
};

// Matrix operations
auto At = A.transpose();
auto B = A + A;
auto C = A * A;

// Linear system solving
cp::d_vector<double> b = {1, 2, 3};
auto x = cp::solve(A, b);  // Solve Ax = b

// Determinant
double det = cp::determinant(A);

// Identity matrix
auto I = cp::d_matrix<double>::identity(3);
```

## Number Theory

### Modular Arithmetic

**Header**: `nt/modular_arithmetic.h`

#### Static Modulus

```cpp
template<integer Mod>
class modular {
public:
    // Constructors
    modular();
    modular(integer value);

    // Get value
    integer value() const;

    // Arithmetic operations
    modular operator+(const modular& other) const;
    modular operator-(const modular& other) const;
    modular operator*(const modular& other) const;
    modular operator/(const modular& other) const;

    // Inverse
    modular inverse() const;

    // Power
    modular pow(integer exp) const;
};

// Common moduli
using mod1e9_7 = modular<1000000007>;
using mod998244353 = modular<998244353>;
```

#### Dynamic Modulus

```cpp
class dynamic_modular {
public:
    dynamic_modular(integer value, integer modulus);

    // Similar interface to static modular
    static void set_modulus(integer mod);
};
```

#### Usage Example

```cpp
#include "nt/modular_arithmetic.h"

// Static modulus
using mint = cp::modular<1000000007>;

mint a = 123456789;
mint b = 987654321;

mint c = a + b;
mint d = a * b;
mint e = a.pow(12345);
mint f = a.inverse();  // Modular inverse
mint g = a / b;        // Modular division

// Factorial modulo p
mint factorial(int n) {
    mint result = 1;
    for (int i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}
```

### Number Theory Functions

**Header**: `nt/number_theory.h`

```cpp
// Greatest Common Divisor
template<typename T>
T gcd(T a, T b);

// Least Common Multiple
template<typename T>
T lcm(T a, T b);

// Extended Euclidean Algorithm
template<typename T>
egcd_result<T> extended_gcd(T a, T b);

// Modular inverse
integer mod_inverse(integer a, integer mod);

// Chinese Remainder Theorem
integer crt(const std::vector<integer>& remainders,
            const std::vector<integer>& moduli);
```

### Primality Testing

**Header**: `nt/primality.h`

```cpp
// Miller-Rabin primality test
bool is_prime(integer n, int iterations = 20);

// Pollard's Rho factorization
std::vector<integer> factorize(integer n);

// Prime factorization with multiplicities
std::map<integer, integer> prime_factorization(integer n);
```

### Sieve of Eratosthenes

**Header**: `nt/sieve.h`

```cpp
class prime_sieve {
public:
    prime_sieve(integer n);

    // Check if prime
    bool is_prime(integer p) const;

    // Get all primes up to n
    const std::vector<integer>& primes() const;

    // Smallest prime factor
    integer spf(integer n) const;

    // Prime factorization
    std::vector<integer> factorize(integer n) const;
};
```

#### Usage Example

```cpp
#include "nt/sieve.h"

// Generate primes up to 1,000,000
cp::prime_sieve sieve(1000000);

bool is_17_prime = sieve.is_prime(17);
const auto& primes = sieve.primes();

// Factor 12345
auto factors = sieve.factorize(12345);
```

### Multiplicative Functions

**Header**: `nt/dirichelet.h`

```cpp
// Euler's totient function
integer euler_phi(integer n);

// Number of divisors
integer divisor_count(integer n);

// Sum of divisors
integer divisor_sum(integer n);

// Möbius function
integer mobius(integer n);

// Dirichlet convolution
template<typename F, typename G>
auto convolve(F f, G g);
```

## Polynomials

### Polynomial Class

**Header**: `polynomial/polynomial.h`

```cpp
template<typename T>
class polynomial {
public:
    // Constructors
    polynomial();
    polynomial(const std::vector<T>& coefficients);
    polynomial(std::initializer_list<T> coefficients);

    // Degree
    size_t degree() const;

    // Coefficient access
    T& operator[](size_t i);
    const T& operator[](size_t i) const;

    // Evaluation
    T eval(const T& x) const;

    // Arithmetic operations
    polynomial operator+(const polynomial& other) const;
    polynomial operator-(const polynomial& other) const;
    polynomial operator*(const polynomial& other) const;
    polynomial operator*(const T& scalar) const;

    // Division (over fields)
    polynomial operator/(const polynomial& other) const;
    polynomial operator%(const polynomial& other) const;

    // Derivative
    polynomial derivative() const;

    // Integration
    polynomial integral() const;
};
```

#### Usage Example

```cpp
#include "polynomial/polynomial.h"

// Define polynomial 1 + 2x + 3x²
cp::polynomial<double> p = {1, 2, 3};

// Evaluate at x = 5
double result = p.eval(5.0);  // 1 + 2*5 + 3*25 = 86

// Polynomial arithmetic
cp::polynomial<double> q = {4, 5};
auto sum = p + q;
auto product = p * q;

// Derivative
auto dp = p.derivative();  // 2 + 6x
```

### Fast Polynomial Multiplication

**Header**: `polynomial/fft.h`

```cpp
// Multiply polynomials using FFT
template<typename T>
std::vector<T> multiply_fft(const std::vector<T>& a,
                             const std::vector<T>& b);

// For modular arithmetic (NTT)
template<integer Mod>
std::vector<modular<Mod>> multiply_ntt(
    const std::vector<modular<Mod>>& a,
    const std::vector<modular<Mod>>& b);
```

### Formal Power Series

**Header**: `polynomial/formal_series.h`

```cpp
template<typename T>
class formal_series {
public:
    // Exponential
    formal_series exp() const;

    // Logarithm
    formal_series log() const;

    // Power
    formal_series pow(integer k) const;

    // Inverse
    formal_series inverse() const;

    // Composition
    formal_series compose(const formal_series& other) const;
};
```

## FFT & Signal Processing

### Fast Fourier Transform

**Header**: `signals/fft.h`

```cpp
// 1D FFT
template<typename T>
std::vector<std::complex<T>> fft(const std::vector<T>& input);

template<typename T>
std::vector<T> ifft(const std::vector<std::complex<T>>& input);

// 2D FFT
template<typename T>
std::vector<std::vector<std::complex<T>>> fft2(
    const std::vector<std::vector<T>>& input);

// Multidimensional FFT
template<typename T>
class multi_fft {
public:
    // FFT on tensors
    tensor<std::complex<T>> forward(const tensor<T>& input);
    tensor<T> inverse(const tensor<std::complex<T>>& input);
};
```

### Number Theoretic Transform

**Header**: `signals/ntt.h`

```cpp
// NTT for modular arithmetic
template<integer Mod>
std::vector<modular<Mod>> ntt(const std::vector<modular<Mod>>& input);

template<integer Mod>
std::vector<modular<Mod>> intt(const std::vector<modular<Mod>>& input);
```

### Fast Hadamard Transform

```cpp
// Fast Hadamard Transform
template<typename T>
std::vector<T> fht(const std::vector<T>& input);
```

## Graph Algorithms

### Tree Operations

**Header**: `graph/tree/tree.h`

```cpp
template<typename T>
class tree {
public:
    // Add edge
    void add_edge(int u, int v, const T& weight = T(1));

    // Root the tree
    void root(int r);

    // Parent of node
    int parent(int u) const;

    // Children of node
    const std::vector<int>& children(int u) const;

    // Depth of node
    int depth(int u) const;

    // Subtree size
    int subtree_size(int u) const;
};
```

### Heavy-Light Decomposition

**Header**: `graph/tree/range_queries.h`

```cpp
template<typename T, typename BinaryOp>
class heavy_light_decomposition {
public:
    heavy_light_decomposition(const tree<T>& t, BinaryOp op);

    // Path query from u to v
    T path_query(int u, int v) const;

    // Path update from u to v
    void path_update(int u, int v, const T& value);

    // Subtree query
    T subtree_query(int u) const;
};
```

#### Usage Example

```cpp
#include "graph/tree/tree.h"
#include "graph/tree/range_queries.h"

// Build tree
cp::tree<int> t(10);
t.add_edge(0, 1, 5);
t.add_edge(0, 2, 3);
t.add_edge(1, 3, 2);
// ... more edges

t.root(0);

// Heavy-Light Decomposition
cp::heavy_light_decomposition<int, std::plus<>> hld(t, std::plus<>());

int path_sum = hld.path_query(3, 5);
hld.path_update(3, 5, 10);
```

### Union-Find (Disjoint Set Union)

**Header**: `graph/union_find.h`

```cpp
class union_find {
public:
    union_find(size_t n);

    // Find root with path compression
    int find(int x);

    // Union by rank
    void unite(int x, int y);

    // Check if same set
    bool same(int x, int y);

    // Number of components
    int count_components() const;
};
```

### 2-SAT Solver

**Header**: `graph/2sat.h`

```cpp
class two_sat {
public:
    two_sat(int n);  // n variables

    // Add clause (x ∨ y)
    void add_clause(int x, bool x_value, int y, bool y_value);

    // Solve
    bool solve();

    // Get assignment
    bool value(int x) const;
};
```

## Combinatorics

### Binomial Coefficients

**Header**: `combinatorics/binomial.h`

```cpp
template<typename T>
class binomial_table {
public:
    binomial_table(size_t max_n);

    // Get C(n, k)
    T operator()(size_t n, size_t k) const;
};

// For modular arithmetic
template<integer Mod>
binomial_table<modular<Mod>> make_binomial_mod(size_t max_n);
```

### Factorial

**Header**: `combinatorics/factorial.h`

```cpp
template<typename T>
class factorial_table {
public:
    factorial_table(size_t max_n);

    // Get n!
    T factorial(size_t n) const;

    // Get 1/n!
    T inverse_factorial(size_t n) const;
};
```

### Stirling Numbers

**Header**: `combinatorics/stirling.h`

```cpp
// Stirling numbers of the first kind
template<typename T>
std::vector<std::vector<T>> stirling_first(size_t max_n);

// Stirling numbers of the second kind
template<typename T>
std::vector<std::vector<T>> stirling_second(size_t max_n);
```

### Partitions

**Header**: `combinatorics/partitions.h`

```cpp
// Number of partitions of n
template<typename T>
T partition_count(size_t n);

// Partition function table
template<typename T>
std::vector<T> partition_table(size_t max_n);
```

## Machine Learning

### Linear Regression

**Header**: `ml/ml.h`

```cpp
class linear_regression {
public:
    // Fit model
    linear_regression& fit(const d_matrix<real>& X,
                           const d_vector<real>& y);

    // Predict
    d_vector<real> predict(const d_matrix<real>& X) const;

    // Score (R² score)
    real score(const d_matrix<real>& X,
               const d_vector<real>& y) const;

    // Coefficients
    const d_vector<real>& coefficients() const;
    real intercept() const;
};
```

### Logistic Regression

```cpp
class logistic_regression {
public:
    logistic_regression& fit(const d_matrix<real>& X,
                              const d_vector<real>& y,
                              size_t max_iter = 1000,
                              real learning_rate = 0.01);

    d_vector<real> predict_proba(const d_matrix<real>& X) const;
    d_vector<int> predict(const d_matrix<real>& X) const;

    real score(const d_matrix<real>& X, const d_vector<real>& y) const;
};
```

### Multiclass Logistic Regression

```cpp
class multilogistic_regression {
public:
    multilogistic_regression& fit(const d_matrix<real>& X,
                                   const d_vector<real>& y);

    d_vector<real> predict(const d_matrix<real>& X) const;
    real score(const d_matrix<real>& X, const d_vector<real>& y) const;
    real error(const d_matrix<real>& X, const d_vector<real>& y) const;
};
```

### K-Nearest Neighbors

```cpp
template<typename Metric>
class k_nearest_neighbour_classifier {
public:
    size_t k = 5;  // Number of neighbors

    void fit(const d_matrix<real>& X, const d_vector<real>& y);
    d_vector<real> predict(const d_matrix<real>& X) const;
    real score(const d_matrix<real>& X, const d_vector<real>& y) const;
};

template<typename Metric>
class k_nearest_neighbour_regressor {
    // Similar interface for regression
};
```

#### Usage Example

```cpp
#include "ml/ml.h"
#include "topology/topology.h"

// Load data
cp::d_matrix<real> X_train = /* ... */;
cp::d_vector<real> y_train = /* ... */;

// Linear regression
cp::linear_regression lr;
lr.fit(X_train, y_train);
auto y_pred = lr.predict(X_train);
real score = lr.score(X_train, y_train);

// KNN with L2 norm
cp::k_nearest_neighbour_classifier<cp::L2_norm<cp::d_vector<real>>> knn;
knn.k = 5;
knn.fit(X_train, y_train);
auto predictions = knn.predict(X_train);
```

## Functional Programming

### Functional Operations

**Header**: `functional/functional.h`

```cpp
// Map operation
template<typename Func, typename Container>
auto map(Func f, const Container& c);

// Filter operation
template<typename Pred, typename Container>
auto filter(Pred p, const Container& c);

// Reduce operation
template<typename BinaryOp, typename Container>
auto reduce(BinaryOp op, const Container& c);

// Foreach
template<typename Func, typename Container>
void foreach(Func f, Container& c);
```

### Zip

**Header**: `functional/zip.h`

```cpp
// Zip multiple containers
template<typename... Containers>
auto zip(Containers&... containers);

// Usage in range-based for
for (auto [a, b, c] : cp::zip(vec1, vec2, vec3)) {
    // Process tuples
}
```

### Pointwise Operations

```cpp
// Pointwise unary operation
template<typename Func, typename Vector>
auto pointwise(Func f, const Vector& v);

// Pointwise binary operation
template<typename BinaryOp, typename Vector>
auto pointwise(BinaryOp op, const Vector& v1, const Vector& v2);
```

## Topology & Optimization

### Metrics and Norms

**Header**: `topology/topology.h`

```cpp
// L1 norm (Manhattan distance)
template<typename Vector>
struct L1_norm {
    typename Vector::value_type operator()(const Vector& v) const;
};

// L2 norm (Euclidean distance)
template<typename Vector>
struct L2_norm {
    typename Vector::value_type operator()(const Vector& v) const;
};

// L-infinity norm
template<typename Vector>
struct Linf_norm {
    typename Vector::value_type operator()(const Vector& v) const;
};

// Inner product
template<typename T, typename Vector>
struct L2_inner_product {
    T operator()(const Vector& v1, const Vector& v2) const;
};
```

### Optimization

**Header**: `topology/optimisation.h`

```cpp
// Newton-Raphson method
template<typename Func, typename Derivative>
real newton_raphson(Func f, Derivative df, real x0,
                    real tol = 1e-10, size_t max_iter = 100);

// Gradient descent
template<typename Func, typename Gradient>
d_vector<real> gradient_descent(Func f, Gradient grad,
                                 const d_vector<real>& x0,
                                 real learning_rate = 0.01,
                                 size_t max_iter = 1000);
```

### Simplex Method

**Header**: `topology/simplex.h`

```cpp
class simplex_solver {
public:
    // Maximize c^T x subject to Ax <= b, x >= 0
    simplex_solver(const d_matrix<real>& A,
                   const d_vector<real>& b,
                   const d_vector<real>& c);

    // Solve
    bool solve();

    // Get solution
    const d_vector<real>& solution() const;
    real objective_value() const;
};
```

## String Algorithms

### Rabin-Karp String Matching

**Header**: `string/string.h`

```cpp
class rabin_karp {
public:
    rabin_karp(const std::string& pattern);

    // Find all occurrences
    std::vector<size_t> find_all(const std::string& text) const;

    // Check if pattern exists
    bool exists(const std::string& text) const;
};
```

## Advanced Features

### Parser & Formal Verification

**Header**: `parser/StatefulParser.h`

```cpp
// Grammar definition
class grammar {
public:
    void add_production(const std::string& non_terminal,
                       const std::vector<std::string>& production);

    void set_start_symbol(const std::string& symbol);
};

// LR Parser
class lr_parser {
public:
    lr_parser(const grammar& g);

    // Parse input
    bool parse(const std::vector<std::string>& tokens);
};
```

### Generic Type Traits

Many templates work with custom types if they satisfy certain requirements:

```cpp
// Ring requirements
template<typename T>
concept Ring = requires(T a, T b) {
    { a + b } -> std::convertible_to<T>;
    { a - b } -> std::convertible_to<T>;
    { a * b } -> std::convertible_to<T>;
    { T(0) };  // Additive identity
    { T(1) };  // Multiplicative identity
};

// Field requirements (Ring + division)
template<typename T>
concept Field = Ring<T> && requires(T a, T b) {
    { a / b } -> std::convertible_to<T>;
};
```

## Performance Considerations

### Compile-Time vs Runtime

- Use **static vectors/matrices** (`s_vector`, `s_matrix`) when size is known at compile time
- Use **dynamic vectors/matrices** (`d_vector`, `d_matrix`) for runtime-sized data

### Modular Arithmetic

- Use **static modulus** (`modular<Mod>`) when modulus is known at compile time (faster)
- Use **dynamic modulus** when modulus varies at runtime

### Memory Layout

- Matrices are stored in row-major order for cache efficiency
- Segment trees use implicit array representation for memory efficiency

## Error Handling

The library uses assertions and exceptions for error handling:

```cpp
// Assertions for debug builds
assert(index < size);

// Exceptions for runtime errors
if (matrix_not_invertible) {
    throw std::runtime_error("Matrix is not invertible");
}
```

## Best Practices

1. **Include only what you need** - The library is header-only; include specific headers
2. **Use appropriate data structures** - Choose fixed vs dynamic based on your needs
3. **Leverage compile-time optimization** - Use templates and constexpr when possible
4. **Check mathematical requirements** - Ensure your types satisfy operation requirements
5. **Test thoroughly** - Use the provided test suite as examples

## See Also

- **[README.md](README.md)** - Project overview and quick start
- **[EXAMPLES.md](EXAMPLES.md)** - Detailed usage examples
- **[ARCHITECTURE.md](ARCHITECTURE.md)** - Design decisions and architecture
- **[doc/](doc/)** - Module-specific documentation with mathematical background

---

*For questions or clarifications, please refer to the source code or open an issue.*
