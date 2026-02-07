# Linear Algebra Module Documentation

## Table of Contents
1. [Overview](#overview)
2. [Core Concepts](#core-concepts)
3. [Vector Class](#vector-class)
4. [Matrix Class](#matrix-class)
5. [Matrix Operations](#matrix-operations)
6. [Linear System Solving](#linear-system-solving)
7. [Matrix Decompositions](#matrix-decompositions)
8. [Special Matrices](#special-matrices)
9. [Characteristic and Minimal Polynomials](#characteristic-and-minimal-polynomials)
10. [API Reference](#api-reference)
11. [Usage Examples](#usage-examples)
12. [Complexity Analysis](#complexity-analysis)
13. [Supported Types](#supported-types)

---

## Overview

The Linear Algebra module provides a comprehensive implementation of vectors, matrices, and related operations for competitive programming. It leverages C++23 features including concepts, template metaprogramming, and compile-time optimizations.

### Key Features

- **Generic Design**: Works with any ring or field type (integers, floats, complex numbers, modular arithmetic)
- **Static and Dynamic Sizing**: Compile-time sized arrays for performance, runtime sizing for flexibility
- **Complete Matrix Operations**: Addition, multiplication, transpose, inverse, determinant, rank, nullity
- **Linear System Solvers**: General and specialized solvers for AX = B
- **Matrix Decompositions**: Cholesky, LDL, QR, LQ, SVD, Schur, Polar
- **Orthogonalization**: Gram-Schmidt algorithm with completion
- **Polynomial Computation**: Characteristic and minimal polynomials
- **Special Matrix Types**: Band matrices, triangular matrices
- **Zero-Overhead Abstractions**: Template-based design with no runtime penalty

### Module Location

```
include/linear_algebra/
├── vector.h                  # Vector class and operations
├── matrix.h                  # Matrix class and operations
├── decomposition.h           # Matrix decompositions (Cholesky, QR, SVD, etc.)
├── special_matrices.h        # Band and triangular matrices
├── special_polynomials.h     # Characteristic and minimal polynomials
├── utils.h                   # Utility definitions
├── tensor.h                  # Tensor base class
└── view.h                    # Tensor view abstraction
```

---

## Core Concepts

### Ring and Field Concepts

The module uses C++23 concepts to constrain template parameters:

```cpp
template<typename R>
concept ring = requires(R a, R b) {
    { a + b } -> std::convertible_to<R>;
    { a - b } -> std::convertible_to<R>;
    { a * b } -> std::convertible_to<R>;
    { -a } -> std::convertible_to<R>;
};
```

**Requirements:**
- A **ring** is an algebraic structure with addition, subtraction, multiplication, and additive inverse
- A **field** additionally has division (e.g., real numbers, complex numbers, modular arithmetic with prime modulus)

### Static vs Dynamic Extent

The module supports both compile-time and runtime sizing:

```cpp
constexpr std::size_t dynamic_extent = -1;

// Static extent (compile-time size)
vector<int, 5> v1;              // Size fixed at 5
matrix<double, 3, 3> m1;        // 3×3 matrix

// Dynamic extent (runtime size)
vector<int> v2(10, size_tag);   // Size determined at runtime
matrix<double> m2(3, 4, size_tag);  // 3×4 matrix
```

**Advantages of Static Extent:**
- Compile-time optimizations
- Stack allocation (faster)
- Better cache locality
- Type safety for size mismatches

**Advantages of Dynamic Extent:**
- Flexible sizing
- Memory efficiency for variable-sized data
- Easier interfacing with runtime input

---

## Vector Class

### Structure

```cpp
template<ring R, std::size_t ext = dynamic_extent>
struct vector : tensor_view<R, 1> {
    using container = std::conditional_t<
        ext == dynamic_extent,
        std::vector<R>,
        std::array<R, ext>
    >;

    container u{};

    // Constructors
    vector();
    vector(std::size_t n, size_tag_t);  // Dynamic extent only
    vector(std::initializer_list<R> init);

    // Element access
    R& operator[](std::size_t i);
    const R& operator[](std::size_t i) const;

    // Dimensions
    std::size_t dim() const;
    std::size_t size() const;

    // Operations
    vector& operator+=(const vector& other);
    vector& operator-=(const vector& other);
    vector& operator*=(const R& k);
    vector& operator/=(const R& k);
    vector operator-() const;

    // Iterators
    auto begin();
    auto end();
};
```

### Type Aliases

```cpp
// Dynamic extent vectors
template<typename R>
using d_vector = vector<R, dynamic_extent>;

// Static extent vectors
template<typename R, std::size_t n>
using s_vector = vector<R, n>;
```

### Supported Operations

| Operation | Syntax | Description | Complexity |
|-----------|--------|-------------|------------|
| Addition | `v1 + v2` | Element-wise addition | O(n) |
| Subtraction | `v1 - v2` | Element-wise subtraction | O(n) |
| Scalar Multiplication | `k * v` or `v * k` | Multiply each element by scalar | O(n) |
| Scalar Division | `v / k` | Divide each element by scalar | O(n) |
| Negation | `-v` | Negate all elements | O(n) |
| Inner Product | `inner_product(v1, v2)` | Dot product | O(n) |
| Norm | `norm(v)` | Euclidean norm | O(n) |

---

## Matrix Class

### Structure

```cpp
template<ring R, std::size_t ext1 = dynamic_extent, std::size_t ext2 = ext1>
struct matrix : tensor_view<R, 2> {
    using vector = vector<R, ext2>;
    using container = std::conditional_t<
        ext1 == dynamic_extent,
        std::vector<vector>,
        std::array<vector, ext1>
    >;

    container M{};

    // Constructors
    matrix();
    matrix(std::size_t rows, std::size_t cols, size_tag_t);
    matrix(std::initializer_list<std::initializer_list<R>> init);

    // Factory methods
    static matrix eye(std::size_t n);  // Identity matrix

    // Dimensions
    std::size_t rows() const;
    std::size_t cols() const;
    std::size_t size() const;

    // Element access
    vector& operator[](int i);
    const vector& operator[](int i) const;
    R& at(std::array<std::size_t, 2> I);

    // Basic operations
    matrix& operator+=(const matrix& O);
    matrix& operator-=(const matrix& O);
    matrix& operator*=(const matrix& B);
    matrix& operator*=(R k);
    matrix& operator/=(R k);
    matrix operator-() const;

    // Matrix operations
    matrix operator*(const matrix& B) const;
    vector operator*(const vector& u) const;

    // Linear algebra operations
    R tr() const;                       // Trace
    matrix T() const;                   // Transpose
    matrix H() const;                   // Hermitian conjugate
    matrix inv() const;                 // Inverse
    R det() const;                      // Determinant
    size_t rank() const;                // Rank
    size_t nullity() const;             // Nullity

    // Solving linear systems
    std::optional<matrix> solve(const matrix& O, bool invertible = false) const;
    std::optional<vector> solve(const vector& V, bool invertible = false) const;

    // Basis computation
    matrix null_basis() const;          // Null space basis
    matrix image_basis() const;         // Column space basis
};
```

### Type Aliases

```cpp
// Dynamic extent matrices
template<typename R>
using d_matrix = matrix<R, dynamic_extent, dynamic_extent>;

// Static extent matrices
template<typename R, std::size_t n, std::size_t m>
using s_matrix = matrix<R, n, m>;

// Square matrices
template<typename R, std::size_t n>
using square_matrix = matrix<R, n, n>;
```

---

## Matrix Operations

### Basic Operations

#### Matrix Addition and Subtraction

```cpp
// Element-wise addition
matrix<R> C = A + B;

// In-place addition
A += B;

// Subtraction
matrix<R> D = A - B;
A -= B;

// Special case: scalar as 1×1 matrix
A += k;  // Adds k to diagonal elements
```

**Implementation:**
```cpp
template<ring R, std::size_t ext1, std::size_t ext2>
matrix<R, ext1, ext2>& matrix<R, ext1, ext2>::operator+=(const matrix& O) {
    if (O.rows() == 1 && O.cols() == 1) {
        for (int i = 0; i < std::min(rows(), cols()); i++)
            M[i][i] += O.M[0][0];
    } else {
        for (int i = 0; i < std::min(rows(), O.rows()); i++)
            for (int j = 0; j < std::min(cols(), O.cols()); j++)
                M[i][j] += O.M[i][j];
    }
    return *this;
}
```

**Complexity:**
- Time: O(n × m) for n×m matrices
- Space: O(1) for in-place, O(n × m) for new matrix

#### Matrix Multiplication

Standard matrix multiplication using the definition:
```
C[i][j] = Σ(k=0 to p-1) A[i][k] × B[k][j]
```

```cpp
template<ring R, std::size_t ext1, std::size_t ext2>
matrix<R> matrix<R, ext1, ext2>::operator*(const matrix& B) const {
    auto n = rows(), p = cols(), m = B.cols();
    matrix C(n, m, size_tag);
    for (int i = 0; i < n; i++)
        for (int k = 0; k < p; k++)
            for (int j = 0; j < m; j++)
                C.M[i][j] += M[i][k] * B.M[k][j];
    return C;
}
```

**Complexity:**
- Time: O(n × m × p) for (n×p) × (p×m) matrices
- Space: O(n × m) for result matrix

#### Matrix-Vector Multiplication

```cpp
vector<R> matrix<R>::operator*(const vector& u) const {
    int n = rows(), m = cols();
    vector v(n, size_tag);
    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++)
            v[i] += M[i][j] * u[j];
    return v;
}
```

**Complexity:**
- Time: O(n × m) for n×m matrix
- Space: O(n) for result vector

#### Scalar Multiplication and Division

```cpp
// Scalar multiplication
matrix<R> B = k * A;
matrix<R> C = A * k;
A *= k;

// Scalar division (requires field)
matrix<R> D = A / k;
A /= k;
```

**Complexity:**
- Time: O(n × m)
- Space: O(1) for in-place, O(n × m) for new matrix

### Advanced Operations

#### Transpose

Swaps rows and columns: AT[i][j] = A[j][i]

```cpp
matrix<R> matrix<R>::T() const {
    int m = cols(), n = rows();
    matrix P(m, n, size_tag);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            P.M[j][i] = M[i][j];
    return P;
}
```

**Complexity:**
- Time: O(n × m)
- Space: O(n × m)

#### Hermitian Conjugate

For complex matrices, H = (AT)* (transpose and conjugate each element)

```cpp
matrix<R> matrix<R>::H() const {
    int m = cols(), n = rows();
    matrix P(m, n, size_tag);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            P.M[j][i] = conj<R, R>(M[i][j]);
    return P;
}
```

**Complexity:**
- Time: O(n × m)
- Space: O(n × m)

#### Trace

Sum of diagonal elements: tr(A) = Σ A[i][i]

```cpp
R matrix<R>::tr() const {
    int m = cols(), n = rows();
    R r{};
    for (int i = 0; i < std::min(n, m); i++)
        r += M[i][i];
    return r;
}
```

**Complexity:**
- Time: O(min(n, m))
- Space: O(1)

#### Determinant

Computed using row echelon form with Gaussian elimination.

**Algorithm:**
1. Perform row echelon form (REF) reduction
2. Track row swaps (changes sign)
3. Multiply diagonal elements of REF

```cpp
R matrix<R>::det() const {
    auto dec = row_echelon(matrix<R>(rows(), 0, size_tag));
    R w = dec.dir ? 1 : -1;
    for (int i = 0; i < rows(); i++)
        w *= dec.E[i][i];
    return w;
}
```

**Complexity:**
- Time: O(n³) for n×n matrix
- Space: O(n²)

**Properties:**
- det(AB) = det(A) × det(B)
- det(AT) = det(A)
- det(kA) = kn × det(A) for n×n matrix
- det(A) ≠ 0 ⟺ A is invertible

#### Matrix Inverse

Solves AX = I to find A⁻¹.

```cpp
matrix<R> matrix<R>::inv() const {
    return *solve(matrix::eye(rows()), true);
}
```

**Complexity:**
- Time: O(n³) for n×n matrix
- Space: O(n²)

**Requirements:**
- Matrix must be square (n×n)
- det(A) ≠ 0 (invertible)

**Properties:**
- A × A⁻¹ = I
- (AB)⁻¹ = B⁻¹A⁻¹
- (AT)⁻¹ = (A⁻¹)T

#### Rank and Nullity

**Rank:** Number of linearly independent rows/columns (dimension of column space)

```cpp
size_t matrix<R>::rank() const {
    return row_echelon(matrix<R>(rows(), 0, size_tag)).rank;
}
```

**Nullity:** Dimension of null space (solutions to Ax = 0)

```cpp
size_t matrix<R>::nullity() const {
    return cols() - rank();
}
```

**Rank-Nullity Theorem:**
For an n×m matrix A:
```
rank(A) + nullity(A) = m
```

**Complexity:**
- Time: O(n² × m) for n×m matrix
- Space: O(n × m)

---

## Linear System Solving

The module provides robust solvers for linear systems of the form **AX = B**.

### General Solver

Solves AX = B where:
- A is n×m
- B is n×p
- X is m×p (solution)

```cpp
std::optional<matrix<R>> solve(const matrix<R>& B, bool invertible = false) const;
std::optional<vector<R>> solve(const vector<R>& b, bool invertible = false) const;
```

**Returns:**
- `std::optional` containing solution if exists
- `std::nullopt` if no solution exists

### Algorithm: Gaussian Elimination with Row Echelon Form

#### Step 1: Row Echelon Form (REF)

Transform augmented matrix [A | B] to REF using elementary row operations:

```cpp
template<std::size_t ext3>
MatSolve<ext3> row_echelon(matrix<R, ext2, ext3> Y) const {
    size_t rnk = 0;
    auto E = *this;
    auto n = rows(), m = cols();
    bool dir = true;
    std::vector<std::pair<int, int>> mapper;

    for (int i = 0; i < m && rnk < n; i++) {
        // Find pivot
        int r = pivot_rule(E, rnk, i);
        if (r == n) continue;  // No pivot in this column

        mapper.emplace_back(rnk, i);

        // Swap rows if needed
        if (r != rnk) {
            std::swap(E[rnk], E[r]);
            std::swap(Y[rnk], Y[r]);
            dir = !dir;
        }

        // Eliminate below pivot
        for (int j = rnk + 1; j < n; j++) {
            auto w = E[j][i] / E[rnk][i];
            E[j][i] = 0;
            for (int k = i + 1; k < m; k++)
                E[j][k] -= w * E[rnk][k];
            Y[j] -= w * Y[rnk];
        }
        rnk++;
    }

    return {E, Y, rnk, dir, mapper};
}
```

**Pivot Rule (Gauss-Jordan):**
```cpp
inline static std::function<int(const matrix&, int, int)> pivot_rule =
    [](const matrix& X, int rnk, int col) {
        int r = rnk;
        int n = X.rows();
        while (r < n && is_zero(X[r][col])) r++;
        return r;
    };
```

#### Step 2: Back Substitution (General Case)

```cpp
void solve_general() {
    int n = E.rows(), l = Y.cols();

    // Check for inconsistency
    for (int i = rank; i < n; i++)
        for (int j = 0; j < l; j++)
            if (!is_zero(Y[i][j])) {
                X = std::nullopt;  // No solution
                return;
            }

    // Back substitution
    for (int k = mapper.size() - 1; k >= 0; k--) {
        auto [r, i] = mapper[k];
        for (int j = 0; j < r; j++) {
            auto w = E[j][i] / E[r][i];
            E[j][i] = 0;
            Y[j] -= w * Y[r];
        }
    }

    // Extract solution
    for (auto [r, i] : mapper)
        (*X)[i] = Y[r] / E[r][i];
}
```

#### Step 3: Specialized Solver (Invertible Case)

For square, invertible matrices (det(A) ≠ 0):

```cpp
void solve_invertible() {
    int n = E.rows();
    X = Y;

    // Back substitution
    for (int i = n - 1; i >= 0; i--) {
        if (is_zero(E[i][i]))
            throw std::invalid_argument("Matrix is not invertible");

        for (int j = 0; j < i; j++) {
            auto w = E[j][i] / E[i][i];
            E[j][i] = 0;
            Y[j] -= w * Y[i];
        }
    }

    // Normalize
    for (int i = 0; i < n; i++)
        (*X)[i] = Y[i] / E[i][i];
}
```

### Complexity Analysis

| Operation | Time Complexity | Space Complexity |
|-----------|----------------|------------------|
| REF Formation | O(n² × m) | O(n × m) |
| Back Substitution (General) | O(r × m × p) where r = rank | O(m × p) |
| Back Substitution (Invertible) | O(n² × p) | O(n × p) |
| **Total (General)** | **O(n² × m + r × m × p)** | **O((n + m) × (m + p))** |
| **Total (Invertible)** | **O(n² × (m + p))** | **O(n × (m + p))** |

### Null Space Basis

Find basis vectors e₁, ..., eₖ such that A × eᵢ = 0:

```cpp
matrix<R> null_basis() const {
    int n = rows(), m = cols();
    matrix Z(*this);

    // Augment with identity
    for (int i = 0; i < m; i++) {
        Z.M.emplace_back(m, size_tag);
        Z[rows() + i][i] = 1;
    }

    // Transpose and perform REF
    auto C = Z.T().row_echelon(matrix<R>(m, 0, size_tag)).E;

    // Extract null space vectors
    matrix<R> B;
    for (int i = 0; i < m; i++) {
        if (all_of(C[i].begin(), C[i].begin() + n,
                   [](auto x) { return is_zero(x); })) {
            vector u(m, size_tag);
            for (int j = 0; j < m; j++)
                u[j] = C[i][n + j];
            B.M.push_back(u);
        }
    }
    return B;
}
```

**Complexity:**
- Time: O(m² × n)
- Space: O(m × (n + m))

### Image (Column Space) Basis

Find basis for the column space (range) of A:

```cpp
matrix image_basis() const {
    auto dec = row_echelon(*this, false);
    dec.E.resize(dec.rank);
    return dec.E;
}
```

**Complexity:**
- Time: O(n² × m)
- Space: O(n × m)

### Example: Solving Linear Systems

```cpp
// Example 1: Square invertible system
matrix<double> A({{1, 2}, {3, 4}});
vector<double> b({5, 6});
auto x = A.solve(b, true);  // Solve Ax = b
// x = [-4, 4.5]

// Example 2: Overdetermined system (no unique solution)
matrix<double> A({{1, 2}, {3, 4}, {5, 6}});  // 3×2
vector<double> b({5, 11, 17});
auto x = A.solve(b);  // May not have solution

// Example 3: Matrix equation AX = B
matrix<double> A({{1, 2}, {3, 4}});
matrix<double> B({{5, 6}, {7, 8}});
auto X = A.solve(B, true);  // Solve AX = B
```

---

## Matrix Decompositions

### Cholesky Decomposition

For positive-definite Hermitian matrix A: **A = L × L\***

where L is lower triangular and L* is its conjugate transpose.

```cpp
template<typename IK>
std::pair<LowerTriangularSquareMatrix<IK>, UpperTriangularSquareMatrix<IK>>
choleskyDecomposition(const d_matrix<IK>& P);
```

**Algorithm:**

```cpp
for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++) {
        L[i][j] = P[i][j];
        for (int k = 0; k < j; k++)
            L[i][j] -= L[i][k] * conj(L[j][k]);

        if (j < i) {
            if (L[j][j] == 0)
                L[i][j] = 0;
            else
                L[i][j] /= conj(L[j][j]);
        } else {
            L[i][j] = sqrt(abs(L[i][j]));
        }
    }
```

**Complexity:**
- Time: O(n³)
- Space: O(n²)

**Requirements:**
- Matrix must be Hermitian: A = A*
- Matrix must be positive-definite: xᵀAx > 0 for all x ≠ 0

**Applications:**
- Efficient solver for positive-definite systems
- Matrix square root computation
- Monte Carlo simulations
- Optimization algorithms

### LDL Decomposition

For Hermitian matrix A: **A = L × D × L\***

where:
- L is lower triangular with 1s on diagonal
- D is diagonal
- L* is conjugate transpose of L

```cpp
template<typename IK>
struct LDL_Decomposition_t {
    LowerTriangularSquareMatrix<IK> L;
    std::vector<IK> D;
    UpperTriangularSquareMatrix<IK> U;
};

template<typename IK>
LDL_Decomposition_t<IK> LDL_Decomposition(const d_matrix<IK>& P);
```

**Algorithm:**

```cpp
for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++) {
        if (j < i) {
            L[i][j] = P[i][j];
            for (int k = 0; k < j; k++)
                L[i][j] -= D[k] * L[i][k] * conj(L[j][k]);

            if (D[j] == 0)
                L[i][j] = 0;
            else
                L[i][j] /= conj(D[j]);
        } else {
            D[i] = P[i][i];
            for (int k = 0; k < i; k++)
                D[i] -= D[k] * L[i][k] * conj(L[i][k]);
        }
    }
U = L.H();
```

**Complexity:**
- Time: O(n³)
- Space: O(n²)

**Advantages over Cholesky:**
- No square roots required
- Works with semi-definite matrices
- Numerically more stable

### Gram-Schmidt Orthogonalization

Converts a set of linearly independent vectors into an orthonormal basis.

```cpp
template<typename IK>
std::vector<d_vector<IK>> gram_schmidt(
    const std::vector<d_vector<IK>>& A,
    real eps = 1e-6
);
```

**Classical Algorithm:**

```
for i = 0 to m-1:
    u_i = A_i
    for j = 0 to i-1:
        u_i = u_i - proj(B_j, A_i)  // Remove component in direction of B_j

    if ||u_i|| > eps:
        B_i = u_i / ||u_i||  // Normalize
    else:
        B_i = 0  // Dependent vector
```

**Projection:**
```cpp
template<typename IK>
d_vector<IK> unit_proj(d_vector<IK> u, d_vector<IK> v) {
    IK w = 0;
    for (int i = 0; i < u.dim(); i++)
        w += conj(u[i]) * v[i];
    return w * u;
}
```

**Complete Gram-Schmidt:**

Adds standard basis vectors to complete the orthonormal basis to full dimension.

```cpp
template<typename IK>
std::vector<d_vector<IK>> complete_gram_schmidt(
    const std::vector<d_vector<IK>>& A,
    real eps = 1e-6
);
```

**Complexity:**
- Time: O(m × n²) for m vectors in ℝⁿ
- Space: O(m × n)

**Applications:**
- QR decomposition
- Orthonormal basis construction
- Projection onto subspaces
- Numerical stability in computations

### QR Decomposition

For any matrix A: **A = Q × R**

where:
- Q is orthogonal (QᵀQ = I)
- R is upper triangular

```cpp
template<typename IK>
struct QR_Decomposition_t {
    d_matrix<IK> Q, R;
};

template<typename IK>
QR_Decomposition_t<IK> QR_Decomposition(const d_matrix<IK>& A);
```

**Algorithm:**

1. Compute R = Cholesky(AᵀA)
2. Compute Q = A × R⁻¹
3. Orthogonalize Q using complete Gram-Schmidt

```cpp
auto L = choleskyDecomposition(A.H() * A).second;
auto R = d_matrix<IK>(L.M);
auto Q = A * L.pinv().M;

// Extract columns of Q
for (int i = 0; i < m; i++) {
    G.emplace_back(v_shape{n});
    for (int j = 0; j < n; j++)
        G.back()[j] = Q[j][i];
}

// Orthogonalize
auto B = complete_gram_schmidt(G);
```

**Complexity:**
- Time: O(n² × m) for n×m matrix
- Space: O(n × m)

**Applications:**
- Least squares problems
- Eigenvalue computation (QR algorithm)
- Linear regression
- Numerical stability

### LQ Decomposition

For any matrix A: **A = L × Q**

where:
- L is lower triangular
- Q is orthogonal

Similar to QR but uses rows instead of columns.

```cpp
template<typename IK>
struct LQ_Decomposition_t {
    d_matrix<IK> L, Q;
};
```

### Singular Value Decomposition (SVD)

For any n×m matrix A: **A = U × Σ × Vᵀ**

where:
- U is n×n orthogonal matrix (left singular vectors)
- Σ is n×m diagonal matrix (singular values)
- V is m×m orthogonal matrix (right singular vectors)

```cpp
template<typename IK>
struct SVD_t {
    d_matrix<IK> U, D, V;
};

template<typename IK>
SVD_t<IK> SVD(d_matrix<IK> A, int iter);
```

**Algorithm (Iterative QR):**

```cpp
while (iter--) {
    auto [Q, R] = QR_Decomposition(A);
    auto [L, P] = LQ_Decomposition(R);
    A = L;
    U = U * Q;
    V = P * V;
}
D = A;
```

**Complexity:**
- Time: O(iter × n² × m) for n×m matrix
- Space: O(n × m)

**Applications:**
- Principal Component Analysis (PCA)
- Data compression
- Matrix approximation
- Pseudoinverse computation
- Image processing

### Pseudoinverse (Moore-Penrose Inverse)

For any matrix A (not necessarily square or invertible): **A⁺**

Computed using SVD:
```
A = UΣVᵀ
A⁺ = VΣ⁺Uᵀ
```

where Σ⁺ has reciprocals of non-zero singular values.

```cpp
template<typename IK>
d_matrix<IK> pinv(const d_matrix<IK>& A, int iter, real eps = 1e-7) {
    auto [U, D, V] = SVD(A, iter);
    for (int i = 0; i < std::min(A.col_dim(), A.row_dim()); i++)
        if (std::abs(D[i][i]) > eps)
            D[i][i] = IK(1) / D[i][i];
    return V.H() * D * U.H();
}
```

**Complexity:**
- Time: O(iter × n² × m)
- Space: O(n × m)

**Properties:**
- AA⁺A = A
- A⁺AA⁺ = A⁺
- (AA⁺)ᵀ = AA⁺
- (A⁺A)ᵀ = A⁺A

### Schur Decomposition

For square matrix A: **A = QTQᵀ**

where:
- Q is orthogonal
- T is upper triangular (Schur form)

```cpp
template<typename IK>
struct SchurDecomposition_t {
    d_matrix<IK> P, T;
};

template<typename IK>
SchurDecomposition_t<IK> QR_Algorithm(d_matrix<IK> A, int iter);
```

**Algorithm (QR Iteration):**

```cpp
d_matrix<IK> U = matrix::eye(A.rows());
while (iter--) {
    auto [Q, R] = QR_Decomposition(A);
    A = R * Q;
    U *= Q;
}
return {U, A};
```

**Complexity:**
- Time: O(iter × n³)
- Space: O(n²)

**Applications:**
- Eigenvalue computation (diagonal of T)
- Matrix functions
- Stability analysis

### Polar Decomposition

For any matrix A: **A = U × P**

where:
- U is unitary
- P is positive semi-definite Hermitian

```cpp
template<typename IK>
struct PolarDecomposition_t {
    d_matrix<IK> U, P;
};

template<typename IK>
PolarDecomposition_t<IK> polarDecomposition(const d_matrix<IK>& A, int iter);
```

**Algorithm:**

```cpp
auto [U, D, V] = SVD(A, iter);
return {U, D * V};
```

**Complexity:**
- Time: O(iter × n³)
- Space: O(n²)

### Matrix Square Root

For positive semi-definite matrix A: **√A**

```cpp
template<typename IK>
d_matrix<IK> sqrt(const d_matrix<IK>& A, int iter);
```

**Algorithm:**

```cpp
auto [Q, D] = QR_Algorithm(A, iter);
for (int i = 0; i < std::min(D.col_dim(), D.row_dim()); i++)
    D[i][i] = std::sqrt(D[i][i]);
return Q * D * Q.H();
```

### Matrix Functions

General framework for applying functions to symmetric matrices:

```cpp
template<typename IK, typename Function>
d_matrix<IK> symmetric_matrix_function(
    const d_matrix<IK>& A,
    const Function& F,
    int iter
);
```

**Algorithm:**

```cpp
auto [Q, D] = QR_Algorithm(A, iter);
for (int i = 0; i < std::min(D.col_dim(), D.row_dim()); i++)
    D[i][i] = F(D[i][i]);
return Q * D * Q.H();
```

**Examples:**
- exp(A) - Matrix exponential
- log(A) - Matrix logarithm
- sin(A), cos(A) - Trigonometric functions
- A^α - Fractional powers

---

## Special Matrices

### Band Matrix

Matrix with non-zero elements only near the diagonal.

```cpp
template<typename R>
struct BandMatrix {
    std::vector<std::vector<R>> M;
    int n;  // Size
    int r;  // Bandwidth

    R operator()(int i, int j) const;
    R& at(int i, int j);

    std::vector<R> solve(std::vector<R> u) const;
    R det() const;
};
```

**Storage:**
Only stores elements within bandwidth r:
- Elements M[i][j] where |i - j| ≤ r

**Efficient Operations:**
```cpp
// Solve Ax = b for band matrix
std::vector<R> BandMatrix::solve(std::vector<R> u) const {
    auto C = *this;

    // Forward elimination (only within band)
    for (int i = 0; i < n; i++)
        for (int p = 1; p <= r && i + p < n; p++) {
            auto w = C(i + p, i) / C(i, i);
            for (int q = i; q <= i + p + r; q++)
                C.at(i + p, q) -= w * C(i, q);
            u[i + p] -= w * u[i];
        }

    // Back substitution
    for (int i = n - 1; i >= 0; i--) {
        auto w = C(i, i);
        u[i] /= w;
        for (int s = 1; s <= r && i - s >= 0; s++)
            u[i - s] -= C(i - s, i) * u[i];
    }
    return u;
}
```

**Complexity:**
- Solve: O(n × r²) instead of O(n³)
- Determinant: O(n × r²)
- Space: O(n × r) instead of O(n²)

**Applications:**
- Finite difference methods
- Spline interpolation
- Tridiagonal systems

### Triangular Matrices

#### Lower Triangular Matrix

Non-zero elements only on or below diagonal.

```cpp
template<typename R>
using LowerTriangularSquareMatrix = TriangularSquareMatrix<R, false>;
```

#### Upper Triangular Matrix

Non-zero elements only on or above diagonal.

```cpp
template<typename R>
using UpperTriangularSquareMatrix = TriangularSquareMatrix<R, true>;
```

#### Operations

```cpp
template<typename R, bool isUpper>
struct TriangularSquareMatrix {
    std::vector<std::vector<R>> M;
    int n;

    // Basic operations
    R operator()(int i, int j) const;
    R& at(int i, int j);

    // Matrix operations
    TriangularSquareMatrix<R, !isUpper> T() const;     // Transpose
    TriangularSquareMatrix<R, !isUpper> H() const;     // Hermitian
    TriangularSquareMatrix inv() const;                // Inverse
    TriangularSquareMatrix pinv() const;               // Pseudoinverse

    // Linear system
    std::vector<R> solve(std::vector<R> u) const;
    R det() const;
};
```

**Efficient Solver:**

For upper triangular:
```cpp
for (int i = n - 1; i >= 0; i--) {
    if (!is_zero(M[i][i])) {
        u[i] /= M[i][i];
        for (int j = i - 1; j >= 0; j--)
            u[j] -= M[j][i] * u[i];
    } else {
        u[i] = 0;
    }
}
```

For lower triangular:
```cpp
for (int i = 0; i < n; i++) {
    if (!is_zero(M[i][i])) {
        u[i] /= M[i][i];
        for (int j = i + 1; j < n; j++)
            u[j] -= M[j][i] * u[i];
    }
}
```

**Complexity:**
- Solve: O(n²) instead of O(n³)
- Inverse: O(n²)
- Determinant: O(n) (product of diagonal)
- Space: O(n²) but can be optimized to O(n(n+1)/2)

**Properties:**
- Product of triangular matrices is triangular
- Inverse of triangular matrix is triangular
- det(T) = ∏ T[i][i]

---

## Characteristic and Minimal Polynomials

### Characteristic Polynomial

For n×n matrix A, the characteristic polynomial is:
```
p(λ) = det(A - λI)
```

Roots are the eigenvalues of A.

#### Faddeev-Leverrier Algorithm

Computes characteristic polynomial using traces.

```cpp
template<typename R>
polynomial<R> faddev_lerrier_characteristic_polynomial(const d_matrix<R>& A);
```

**Algorithm:**

```cpp
std::vector<R> S(n + 1);
S[n] = 1;
d_matrix<R> C(0, m_shape{n, n});

for (int i = n - 1; i >= 0; i--) {
    // Add S[i+1] to diagonal
    for (int j = 0; j < n; j++)
        C[j][j] += S[i + 1];

    // Multiply by A
    C = A * C;

    // Compute next coefficient
    S[i] = -C.tr() / R(n - i);
}

return S;
```

**Complexity:**
- Time: O(n⁴)
- Space: O(n²)

**Formula:**
```
S_k = -1/k × tr(A × C_{k-1})
where C_k = A^k + S_{k+1}A^{k-1} + ... + S_n I
```

#### Interpolation Method

Uses Newton interpolation on det(A - λI).

```cpp
template<typename R>
polynomial<R> interpolation_characteristic_polynomial(d_matrix<R> M);
```

**Algorithm:**

```cpp
std::vector<R> X(n + 1), Y(n + 1);
for (int i = 0; i <= n; i++) {
    X[i] = i;
    Y[i] = M.det();
    for (int j = 0; j < n; j++)
        M[j][j] = M[j][j] - 1;  // M := M - I
}
return newton_interpolation(X, Y);
```

**Complexity:**
- Time: O(n⁴) (n+1 determinants)
- Space: O(n²)

**Advantage:**
- Simpler implementation
- Works well with polynomial interpolation

### Minimal Polynomial

The minimal polynomial μ(λ) is the monic polynomial of smallest degree such that:
```
μ(A) = 0
```

It divides the characteristic polynomial.

#### Algorithm

```cpp
template<typename IK>
polynomial<IK> minimal_polynomial(const d_matrix<IK>& T, const d_vector<IK>& u);
```

**Steps:**

1. Find smallest degree d such that vectors {u, Tu, T²u, ..., Tᵈu} are linearly dependent
2. Find coefficients c₀, ..., cᵈ such that:
   ```
   c₀u + c₁Tu + ... + cᵈTᵈu = 0
   ```
3. Minimal polynomial is: μ(λ) = c₀ + c₁λ + ... + cᵈλᵈ

```cpp
// Find degree
auto d = *std::upper_bound(D.begin(), D.end(), 0,
    [&T, &u](const auto& x, const auto& y) {
        return annihilable(T, u, x) < annihilable(T, u, y);
    });

// Build vectors U[i] = T^i × u
for (int i = 1; i <= d; i++)
    U[i] = T * U[i - 1];

// Gaussian elimination to find coefficients
// [Implementation details...]

return Z[mapper[d]] / Z[mapper[d]][Z[mapper[d]].degree()];
```

**For full matrix:**

```cpp
template<typename IK>
polynomial<IK> minimal_polynomial(const d_matrix<IK>& T) {
    std::vector<d_vector<IK>> E;
    int n = T.row_dim();

    // Use standard basis vectors
    for (int i = 0; i < n; i++) {
        E.emplace_back(v_shape{n});
        E.back()[i] = 1;
    }

    return minimal_polynomial(T, E);
}
```

**Complexity:**
- Time: O(n⁴)
- Space: O(n²)

**Properties:**
- μ(λ) divides p(λ) (characteristic polynomial)
- Roots of μ(λ) are eigenvalues of A
- deg(μ) ≤ n
- For diagonalizable matrices: μ(λ) = ∏(λ - λᵢ) where λᵢ are distinct eigenvalues

---

## API Reference

### Vector API

```cpp
namespace cp::linalg {

// Vector class
template<ring R, std::size_t ext = dynamic_extent>
struct vector {
    // Constructors
    vector();
    vector(std::size_t n, size_tag_t);
    vector(std::initializer_list<R> init);

    // Dimensions
    std::size_t dim() const;
    std::size_t size() const;

    // Element access
    R& operator[](std::size_t i);
    const R& operator[](std::size_t i) const;

    // Arithmetic operations
    vector& operator+=(const vector& other);
    vector& operator-=(const vector& other);
    vector& operator*=(const R& k);
    vector& operator/=(const R& k);
    vector operator-() const;

    // Comparison
    bool operator==(const vector& other) const;

    // Iterators
    auto begin();
    auto end();
    auto begin() const;
    auto end() const;
};

// Free functions
template<ring R, std::size_t ext>
vector<R, ext> operator+(const vector<R, ext>& a, const vector<R, ext>& b);

template<ring R, std::size_t ext>
vector<R, ext> operator-(const vector<R, ext>& a, const vector<R, ext>& b);

template<ring R, std::size_t ext>
vector<R, ext> operator*(const R& k, const vector<R, ext>& v);

template<ring R, std::size_t ext>
vector<R, ext> operator*(const vector<R, ext>& v, const R& k);

template<ring R, std::size_t ext>
vector<R, ext> operator/(const vector<R, ext>& v, const R& k);

// Type aliases
template<typename R>
using d_vector = vector<R, dynamic_extent>;

template<typename R, std::size_t n>
using s_vector = vector<R, n>;

} // namespace cp::linalg
```

### Matrix API

```cpp
namespace cp::linalg {

// Matrix class
template<ring R, std::size_t ext1 = dynamic_extent, std::size_t ext2 = ext1>
struct matrix {
    // Constructors
    matrix();
    matrix(std::size_t rows, std::size_t cols, size_tag_t);
    matrix(std::initializer_list<std::initializer_list<R>> init);

    // Factory methods
    static matrix eye(std::size_t n);  // Identity matrix

    // Dimensions
    std::size_t rows() const;
    std::size_t cols() const;
    std::size_t size() const;
    std::array<std::size_t, 2> shape() const;

    // Element access
    vector<R, ext2>& operator[](int i);
    const vector<R, ext2>& operator[](int i) const;
    R& at(std::array<std::size_t, 2> I);
    const R& at(std::array<std::size_t, 2> I) const;

    // Basic operations
    matrix& operator+=(const matrix& O);
    matrix& operator-=(const matrix& O);
    matrix& operator*=(const matrix& B);
    matrix& operator*=(R k);
    matrix& operator/=(R k);
    matrix operator-() const;

    // Matrix operations
    matrix operator*(const matrix& B) const;
    vector<R, ext1> operator*(const vector<R, ext2>& u) const;
    matrix operator/(const matrix& O) const;

    // Linear algebra operations
    R tr() const;                    // Trace
    matrix<R, ext2, ext1> T() const; // Transpose
    matrix<R, ext2, ext1> H() const; // Hermitian conjugate
    matrix inv() const;              // Inverse
    R det() const;                   // Determinant
    size_t rank() const;             // Rank
    size_t nullity() const;          // Nullity

    // Solving linear systems
    template<std::size_t ext3>
    std::optional<matrix<R, ext2, ext3>> solve(
        const matrix<R, ext1, ext3>& O,
        bool invertible = false
    ) const;

    std::optional<vector<R, ext2>> solve(
        const vector<R, ext1>& V,
        bool invertible = false
    ) const;

    // Basis computation
    matrix<R> null_basis() const;
    matrix image_basis() const;

    // Comparison
    bool operator==(const matrix& O) const;

    // Iterators
    auto begin();
    auto end();
    auto begin() const;
    auto end() const;

    // Pivot rule (customizable)
    inline static std::function<int(const matrix&, int, int)> pivot_rule;
};

// Free functions
template<ring R, std::size_t ext1, std::size_t ext2>
matrix<R, ext1, ext2> operator+(const matrix<R, ext1, ext2>& A, const matrix<R, ext1, ext2>& B);

template<ring R, std::size_t ext1, std::size_t ext2>
matrix<R, ext1, ext2> operator-(const matrix<R, ext1, ext2>& A, const matrix<R, ext1, ext2>& B);

template<ring R, std::size_t ext1, std::size_t ext2>
matrix<R, ext1, ext2> operator*(const R& k, const matrix<R, ext1, ext2>& A);

template<ring R, std::size_t ext1, std::size_t ext2>
matrix<R, ext1, ext2> operator/(const matrix<R, ext1, ext2>& A, const R& k);

// Type aliases
template<typename R>
using d_matrix = matrix<R, dynamic_extent, dynamic_extent>;

template<typename R, std::size_t n, std::size_t m>
using s_matrix = matrix<R, n, m>;

} // namespace cp::linalg
```

### Decomposition API

```cpp
namespace cp::linalg {

// Cholesky decomposition: A = L × L*
template<typename IK>
std::pair<LowerTriangularSquareMatrix<IK>, UpperTriangularSquareMatrix<IK>>
choleskyDecomposition(const d_matrix<IK>& P);

// LDL decomposition: A = L × D × L*
template<typename IK>
struct LDL_Decomposition_t {
    LowerTriangularSquareMatrix<IK> L;
    std::vector<IK> D;
    UpperTriangularSquareMatrix<IK> U;
};

template<typename IK>
LDL_Decomposition_t<IK> LDL_Decomposition(const d_matrix<IK>& P);

// Gram-Schmidt orthogonalization
template<typename IK>
std::vector<d_vector<IK>> gram_schmidt(
    const std::vector<d_vector<IK>>& A,
    real eps = 1e-6
);

template<typename IK>
std::vector<d_vector<IK>> complete_gram_schmidt(
    const std::vector<d_vector<IK>>& A,
    real eps = 1e-6
);

// QR decomposition: A = Q × R
template<typename IK>
struct QR_Decomposition_t {
    d_matrix<IK> Q, R;
};

template<typename IK>
QR_Decomposition_t<IK> QR_Decomposition(const d_matrix<IK>& A);

// LQ decomposition: A = L × Q
template<typename IK>
struct LQ_Decomposition_t {
    d_matrix<IK> L, Q;
};

template<typename IK>
LQ_Decomposition_t<IK> LQ_Decomposition(const d_matrix<IK>& A);

// Schur decomposition: A = Q × T × Q^T
template<typename IK>
struct SchurDecomposition_t {
    d_matrix<IK> P, T;
};

template<typename IK>
SchurDecomposition_t<IK> QR_Algorithm(d_matrix<IK> A, int iter);

// SVD: A = U × Σ × V^T
template<typename IK>
struct SVD_t {
    d_matrix<IK> U, D, V;
};

template<typename IK>
SVD_t<IK> SVD(d_matrix<IK> A, int iter);

// Polar decomposition: A = U × P
template<typename IK>
struct PolarDecomposition_t {
    d_matrix<IK> U, P;
};

template<typename IK>
PolarDecomposition_t<IK> polarDecomposition(const d_matrix<IK>& A, int iter);

// Pseudoinverse
template<typename IK>
d_matrix<IK> pinv(const d_matrix<IK>& A, int iter, real eps = 1e-7);

// Matrix square root
template<typename IK>
d_matrix<IK> sqrt(const d_matrix<IK>& A, int iter);

// Matrix functions
template<typename IK, typename Function>
d_matrix<IK> symmetric_matrix_function(
    const d_matrix<IK>& A,
    const Function& F,
    int iter
);

} // namespace cp::linalg
```

### Polynomial API

```cpp
namespace cp::linalg {

// Characteristic polynomial
template<typename R>
polynomial<R> faddev_lerrier_characteristic_polynomial(const d_matrix<R>& A);

template<typename R, int n>
polynomial<R> faddev_lerrier_characteristic_polynomial(const s_matrix<R, n, n>& A);

template<typename R>
polynomial<R> interpolation_characteristic_polynomial(d_matrix<R> M);

template<typename R, int n>
polynomial<R> interpolation_characteristic_polynomial(s_matrix<R, n, n> M);

// Minimal polynomial
template<typename IK>
polynomial<IK> minimal_polynomial(const d_matrix<IK>& T, const d_vector<IK>& u);

template<typename IK>
polynomial<IK> minimal_polynomial(const d_matrix<IK>& T, const std::vector<d_vector<IK>>& U);

template<typename IK>
polynomial<IK> minimal_polynomial(const d_matrix<IK>& T);

} // namespace cp::linalg
```

### Special Matrices API

```cpp
namespace cp::linalg {

// Band matrix
template<typename R>
struct BandMatrix {
    std::vector<std::vector<R>> M;
    int n;  // Size
    int r;  // Bandwidth

    BandMatrix(std::vector<std::vector<R>> _M);

    R operator()(int i, int j) const;
    auto& operator[](int i);
    R& at(int i, int j);
    const R& at(int i, int j) const;

    std::vector<R> solve(std::vector<R> u) const;
    R det() const;
};

// Triangular matrices
template<typename R, bool isUpper>
struct TriangularSquareMatrix {
    std::vector<std::vector<R>> M;
    int n;

    TriangularSquareMatrix(R k, int n);
    TriangularSquareMatrix(std::vector<std::vector<R>> _M);

    R operator()(int i, int j) const;
    auto& operator[](int i);
    R& at(int i, int j);
    const R& at(int i, int j) const;

    TriangularSquareMatrix<R, !isUpper> T() const;
    TriangularSquareMatrix<R, !isUpper> H() const;
    TriangularSquareMatrix inv() const;
    TriangularSquareMatrix pinv() const;

    std::vector<R> solve(std::vector<R> u) const;
    R det() const;
};

template<typename R>
using UpperTriangularSquareMatrix = TriangularSquareMatrix<R, true>;

template<typename R>
using LowerTriangularSquareMatrix = TriangularSquareMatrix<R, false>;

} // namespace cp::linalg
```

---

## Usage Examples

### Example 1: Basic Vector Operations

```cpp
#include "linear_algebra/vector.h"
using namespace cp::linalg;

int main() {
    // Static extent vector
    s_vector<int, 5> v1({1, 2, 3, 4, 5});
    s_vector<int, 5> v2({2, 3, 4, 5, 6});

    // Addition
    auto v3 = v1 + v2;  // {3, 5, 7, 9, 11}

    // Scalar multiplication
    auto v4 = 3 * v1;  // {3, 6, 9, 12, 15}

    // In-place operations
    v1 += v2;
    v1 *= 2;

    // Dynamic extent vector
    d_vector<double> v5(10, size_tag);
    for (int i = 0; i < 10; i++)
        v5[i] = i * 0.5;

    return 0;
}
```

### Example 2: Basic Matrix Operations

```cpp
#include "linear_algebra/matrix.h"
using namespace cp::linalg;

int main() {
    // Static extent matrix
    s_matrix<int, 3, 3> A({{1, 2, 3},
                           {4, 5, 6},
                           {7, 8, 9}});

    // Matrix operations
    auto At = A.T();           // Transpose
    int tr = A.tr();           // Trace = 15

    // Matrix multiplication
    s_matrix<int, 3, 3> B({{1, 0, 0},
                           {0, 1, 0},
                           {0, 0, 1}});
    auto C = A * B;  // C = A

    // Matrix-vector multiplication
    s_vector<int, 3> v({1, 2, 3});
    auto w = A * v;  // {14, 32, 50}

    return 0;
}
```

### Example 3: Solving Linear Systems

```cpp
#include "linear_algebra/matrix.h"
#include "nt/modular_arithmetic.h"
using namespace cp::linalg;

constexpr int64_t MOD = 1000000007;
using IF = cyclic<MOD>;

int main() {
    // Solve Ax = b in modular arithmetic
    s_matrix<IF, 3, 3> A({{1, 2, 3},
                          {2, 5, 3},
                          {1, 0, 8}});

    s_vector<IF, 3> b({1, 2, 3});

    // Solve (invertible system)
    auto x = A.solve(b, true);

    if (x.has_value()) {
        // Verify: A * x = b
        auto check = A * (*x);
        assert(check == b);
    }

    return 0;
}
```

### Example 4: Matrix Inverse and Determinant

```cpp
#include "linear_algebra/matrix.h"
#include "nt/modular_arithmetic.h"
using namespace cp::linalg;

constexpr int64_t MOD = 1000000007;
using IF = cyclic<MOD>;

int main() {
    // Create a matrix
    s_matrix<IF, 4, 4> A({{1, 2, 3, 4},
                          {3, 1, 2, 4},
                          {1, 4, 3, 1},
                          {5, 3, 1, 2}});

    // Compute determinant
    IF det = A.det();  // 35
    std::cout << "det(A) = " << det << "\n";

    // Compute inverse (if det != 0)
    if (det != IF(0)) {
        auto Ainv = A.inv();

        // Verify: A * A^(-1) = I
        auto I = A * Ainv;
        auto expected = s_matrix<IF, 4, 4>::eye();
        assert(I == expected);
    }

    return 0;
}
```

### Example 5: Rank and Nullity

```cpp
#include "linear_algebra/matrix.h"
using namespace cp::linalg;

int main() {
    // Create a rank-deficient matrix
    d_matrix<double> A({{1, 2, 3},
                        {2, 4, 6},
                        {3, 6, 9}});

    // All rows are linearly dependent
    size_t r = A.rank();      // 1
    size_t n = A.nullity();   // 2

    std::cout << "rank = " << r << ", nullity = " << n << "\n";

    // Verify rank-nullity theorem: rank + nullity = cols
    assert(r + n == A.cols());

    // Compute null space basis
    auto null_basis = A.null_basis();
    std::cout << "Null space dimension: " << null_basis.rows() << "\n";

    return 0;
}
```

### Example 6: Cholesky Decomposition

```cpp
#include "linear_algebra/decomposition.h"
using namespace cp::linalg;

int main() {
    // Create a positive-definite matrix
    d_matrix<double> A({{4, 2, 1},
                        {2, 5, 3},
                        {1, 3, 6}});

    // Cholesky decomposition: A = L × L^T
    auto [L, U] = choleskyDecomposition(A);

    // Verify
    auto reconstructed = d_matrix<double>(L.M) * d_matrix<double>(U.M);

    // Check if A ≈ L × L^T
    for (int i = 0; i < A.rows(); i++)
        for (int j = 0; j < A.cols(); j++)
            assert(std::abs(A[i][j] - reconstructed[i][j]) < 1e-9);

    return 0;
}
```

### Example 7: QR Decomposition

```cpp
#include "linear_algebra/decomposition.h"
using namespace cp::linalg;

int main() {
    // Create a matrix
    d_matrix<double> A({{1, 2, 3},
                        {4, 5, 6},
                        {7, 8, 9},
                        {10, 11, 12}});

    // QR decomposition: A = Q × R
    auto [Q, R] = QR_Decomposition(A);

    // Q is orthogonal: Q^T × Q = I
    auto QtQ = Q.T() * Q;

    // Verify A = Q × R
    auto reconstructed = Q * R;

    for (int i = 0; i < A.rows(); i++)
        for (int j = 0; j < A.cols(); j++)
            assert(std::abs(A[i][j] - reconstructed[i][j]) < 1e-6);

    return 0;
}
```

### Example 8: Singular Value Decomposition

```cpp
#include "linear_algebra/decomposition.h"
using namespace cp::linalg;

int main() {
    // Create a matrix
    d_matrix<double> A({{1, 2},
                        {3, 4},
                        {5, 6}});

    // SVD: A = U × Σ × V^T
    int iterations = 100;
    auto [U, Sigma, V] = SVD(A, iterations);

    // Singular values are on diagonal of Sigma
    std::cout << "Singular values:\n";
    for (int i = 0; i < std::min(Sigma.rows(), Sigma.cols()); i++)
        std::cout << Sigma[i][i] << " ";
    std::cout << "\n";

    // Verify A ≈ U × Σ × V^T
    auto reconstructed = U * Sigma * V.T();

    return 0;
}
```

### Example 9: Pseudoinverse

```cpp
#include "linear_algebra/decomposition.h"
using namespace cp::linalg;

int main() {
    // Create a non-square matrix
    d_matrix<double> A({{1, 2},
                        {3, 4},
                        {5, 6}});

    // Compute pseudoinverse A^+
    int iterations = 100;
    auto Aplus = pinv(A, iterations);

    // Verify properties:
    // 1. A × A^+ × A = A
    auto check1 = A * Aplus * A;

    // 2. A^+ × A × A^+ = A^+
    auto check2 = Aplus * A * Aplus;

    std::cout << "Pseudoinverse computed successfully\n";

    return 0;
}
```

### Example 10: Characteristic Polynomial

```cpp
#include "linear_algebra/special_polynomials.h"
using namespace cp::linalg;

int main() {
    // Create a matrix
    d_matrix<double> A({{1, 2},
                        {3, 4}});

    // Compute characteristic polynomial: det(A - λI)
    auto p1 = faddev_lerrier_characteristic_polynomial(A);
    auto p2 = interpolation_characteristic_polynomial(A);

    std::cout << "Characteristic polynomial (Faddeev-Leverrier): ";
    for (int i = 0; i <= p1.degree(); i++)
        std::cout << p1[i] << " ";
    std::cout << "\n";

    // For this matrix: p(λ) = λ² - 5λ - 2
    // Eigenvalues are roots of p(λ)

    return 0;
}
```

### Example 11: Minimal Polynomial

```cpp
#include "linear_algebra/special_polynomials.h"
using namespace cp::linalg;

int main() {
    // Create a matrix with repeated eigenvalues
    d_matrix<double> A({{2, 1, 0},
                        {0, 2, 1},
                        {0, 0, 2}});

    // Compute minimal polynomial
    auto mu = minimal_polynomial(A);

    std::cout << "Minimal polynomial: ";
    for (int i = 0; i <= mu.degree(); i++)
        std::cout << mu[i] << " ";
    std::cout << "\n";

    // For Jordan block with eigenvalue 2:
    // mu(λ) = (λ - 2)³

    return 0;
}
```

### Example 12: Gram-Schmidt Orthogonalization

```cpp
#include "linear_algebra/decomposition.h"
using namespace cp::linalg;

int main() {
    // Create vectors
    std::vector<d_vector<double>> vectors;
    vectors.push_back(d_vector<double>({1, 1, 0}));
    vectors.push_back(d_vector<double>({1, 0, 1}));
    vectors.push_back(d_vector<double>({0, 1, 1}));

    // Orthogonalize
    auto orthonormal = gram_schmidt(vectors);

    // Verify orthonormality
    for (int i = 0; i < orthonormal.size(); i++) {
        for (int j = 0; j < orthonormal.size(); j++) {
            double dot = 0;
            for (int k = 0; k < orthonormal[i].dim(); k++)
                dot += orthonormal[i][k] * orthonormal[j][k];

            if (i == j)
                assert(std::abs(dot - 1.0) < 1e-9);  // ||v_i|| = 1
            else
                assert(std::abs(dot) < 1e-9);  // v_i ⊥ v_j
        }
    }

    return 0;
}
```

### Example 13: Band Matrix Solver

```cpp
#include "linear_algebra/special_matrices.h"
using namespace cp::linalg;

int main() {
    // Create a tridiagonal matrix (bandwidth = 1)
    int n = 5;
    std::vector<std::vector<double>> M(n, std::vector<double>(3));

    for (int i = 0; i < n; i++) {
        if (i > 0) M[i][0] = 1.0;      // Lower diagonal
        M[i][1] = 4.0;                  // Main diagonal
        if (i < n - 1) M[i][2] = 1.0;  // Upper diagonal
    }

    BandMatrix<double> A(M);

    // Solve Ax = b
    std::vector<double> b = {1, 2, 3, 4, 5};
    auto x = A.solve(b);

    // O(n × r²) = O(n) for tridiagonal
    std::cout << "Solution: ";
    for (auto xi : x)
        std::cout << xi << " ";
    std::cout << "\n";

    return 0;
}
```

### Example 14: Matrix Power

```cpp
#include "linear_algebra/matrix.h"
#include "nt/modular_arithmetic.h"
using namespace cp::linalg;

constexpr int64_t MOD = 1000000007;
using IF = cyclic<MOD>;

int main() {
    // Compute A^n efficiently using binary exponentiation
    s_matrix<IF, 2, 2> A({{1, 1},
                          {1, 0}});

    // Fibonacci matrix
    int n = 100;
    auto An = pow(A, n);

    // An[0][1] contains F_n (nth Fibonacci number)
    std::cout << "F_" << n << " = " << An[0][1] << "\n";

    return 0;
}
```

### Example 15: Complex Matrices

```cpp
#include "linear_algebra/matrix.h"
#include <complex>
using namespace cp::linalg;

int main() {
    using C = std::complex<double>;

    // Create a complex matrix
    d_matrix<C> A({{C(1, 0), C(0, 1)},
                   {C(0, -1), C(1, 0)}});

    // Hermitian conjugate
    auto Ah = A.H();

    // Check if A is Hermitian: A = A*
    bool is_hermitian = (A == Ah);

    // Matrix operations work the same way
    auto det = A.det();
    auto inv = A.inv();

    std::cout << "det(A) = " << det << "\n";

    return 0;
}
```

---

## Complexity Analysis

### Vector Operations

| Operation | Time Complexity | Space Complexity |
|-----------|----------------|------------------|
| Constructor | O(n) | O(n) |
| Element access | O(1) | O(1) |
| Addition/Subtraction | O(n) | O(n) new / O(1) in-place |
| Scalar multiplication | O(n) | O(n) new / O(1) in-place |
| Inner product | O(n) | O(1) |
| Norm | O(n) | O(1) |

### Matrix Operations

| Operation | n×m Matrix | Square Matrix (n×n) | Notes |
|-----------|-----------|-------------------|-------|
| Constructor | O(n × m) | O(n²) | |
| Element access | O(1) | O(1) | |
| Addition/Subtraction | O(n × m) | O(n²) | |
| Scalar multiplication | O(n × m) | O(n²) | |
| **Matrix multiplication** | **O(n × p × m)** | **O(n³)** | (n×p) × (p×m) |
| Matrix-vector mult | O(n × m) | O(n²) | |
| Transpose | O(n × m) | O(n²) | |
| Trace | O(min(n,m)) | O(n) | |
| **Determinant** | N/A | **O(n³)** | Via Gaussian elim |
| **Inverse** | N/A | **O(n³)** | Via solve |
| **Rank** | **O(n² × m)** | **O(n³)** | Via REF |
| Nullity | O(n² × m) | O(n³) | Same as rank |

### Linear System Solving

| System Type | Time Complexity | Space Complexity | Notes |
|-------------|----------------|------------------|-------|
| General (n×m) | O(n² × m + r × m × p) | O((n+m) × (m+p)) | r = rank, p = # RHS |
| Square invertible | O(n² × (m + p)) | O(n × (m + p)) | |
| Triangular | O(n²) | O(n) | Per RHS |
| Band (bandwidth r) | O(n × r²) | O(n × r) | Much faster |

### Matrix Decompositions

| Decomposition | Time Complexity | Space Complexity | Iterations |
|---------------|----------------|------------------|------------|
| **Cholesky** | **O(n³)** | O(n²) | N/A |
| **LDL** | **O(n³)** | O(n²) | N/A |
| **QR** | **O(n² × m)** | O(n × m) | For n×m matrix |
| **SVD** | **O(iter × n² × m)** | O(n × m) | Iterative |
| **Schur** | **O(iter × n³)** | O(n²) | Iterative |
| Gram-Schmidt | O(m × n²) | O(m × n) | m vectors in ℝⁿ |

### Polynomial Computations

| Operation | Time Complexity | Space Complexity | Notes |
|-----------|----------------|------------------|-------|
| Characteristic poly (Faddeev) | O(n⁴) | O(n²) | |
| Characteristic poly (Interpolation) | O(n⁴) | O(n²) | n+1 determinants |
| Minimal polynomial | O(n⁴) | O(n²) | |

### Space Complexity Summary

**Static Extent (Compile-time):**
- Vector: O(n) on stack
- Matrix: O(n × m) on stack
- Very fast for small sizes (n, m ≤ 100)

**Dynamic Extent (Runtime):**
- Vector: O(n) on heap
- Matrix: O(n × m) on heap
- Necessary for large or variable sizes

---

## Supported Types

### Basic Numeric Types

The module works with any type satisfying the `ring` concept:

```cpp
// Integer types
s_vector<int, 5> v1;
s_vector<long long, 10> v2;

// Floating-point types
d_matrix<float> m1;
d_matrix<double> m2;
d_matrix<long double> m3;

// Complex numbers
#include <complex>
d_matrix<std::complex<double>> m4;
```

### Modular Arithmetic

Full support for modular arithmetic (essential for competitive programming):

```cpp
#include "nt/modular_arithmetic.h"

constexpr int64_t MOD = 1000000007;
using IF = cyclic<MOD>;

// Matrix operations in Z/pZ
s_matrix<IF, 100, 100> A;
auto det = A.det();      // Computed mod p
auto inv = A.inv();      // Modular inverse
```

**Advantages:**
- No overflow issues
- Exact arithmetic
- Required for many competitive programming problems

### Custom Ring Types

Any type implementing ring operations:

```cpp
struct MyRing {
    int value;

    MyRing operator+(const MyRing& other) const {
        return {value + other.value};
    }

    MyRing operator-(const MyRing& other) const {
        return {value - other.value};
    }

    MyRing operator*(const MyRing& other) const {
        return {value * other.value};
    }

    MyRing operator-() const {
        return {-value};
    }

    // ... other required operations
};

d_matrix<MyRing> M;
```

### Precision Handling

**For floating-point types:**

```cpp
// Epsilon for zero checking
inline bool is_zero(double x, double eps = 1e-9) {
    return std::abs(x) < eps;
}

// Used throughout the library
if (is_zero(matrix[i][j])) {
    // Treat as zero
}
```

**Recommendations:**
- Use `double` for most applications (sufficient precision)
- Use `long double` for very sensitive computations
- Use exact types (int, modular) when possible
- Be aware of numerical stability issues with large matrices

### Type Requirements

**For Ring Operations:**
- Addition: `R operator+(const R&, const R&)`
- Subtraction: `R operator-(const R&, const R&)`
- Multiplication: `R operator*(const R&, const R&)`
- Negation: `R operator-(const R&)`
- Zero element: `R{}` or `R()`

**For Field Operations (division required):**
- Division: `R operator/(const R&, const R&)`
- Multiplicative inverse

**Additional Requirements:**
- Equality: `bool operator==(const R&, const R&)` for comparisons
- Copy/move semantics: standard C++ requirements
- Default constructor: `R()`

---

## Performance Tips

### 1. Use Static Extent When Possible

```cpp
// Fast (stack allocation)
s_matrix<double, 100, 100> A;

// Slower (heap allocation)
d_matrix<double> B(100, 100, size_tag);
```

### 2. Mark Invertible Systems

```cpp
// Faster path (O(n²×p) instead of O(n²×m + r×m×p))
auto x = A.solve(b, true);  // true = invertible
```

### 3. Avoid Unnecessary Copies

```cpp
// Bad: creates copies
auto C = A * B * D;

// Better: use references
auto AB = A * B;
auto C = AB * D;

// Best: in-place when possible
A *= B;
```

### 4. Use Specialized Matrices

```cpp
// For tridiagonal systems
BandMatrix<double> A(M);  // O(n) solve
auto x = A.solve(b);

// For triangular systems
UpperTriangularSquareMatrix<double> U;
auto x = U.solve(b);  // O(n²) solve
```

### 5. Choose Appropriate Decomposition

- **Cholesky**: Fastest for positive-definite Hermitian matrices
- **LDL**: No square roots, good for general Hermitian
- **QR**: Numerical stability, orthogonal matrices
- **SVD**: Most general, but slowest (iterative)

### 6. Adjust SVD/Schur Iterations

```cpp
// Fewer iterations for approximate results
auto [U, D, V] = SVD(A, 10);  // Fast but less accurate

// More iterations for precision
auto [U, D, V] = SVD(A, 1000);  // Slow but accurate
```

### 7. Batch Operations

```cpp
// Solve AX = B for multiple RHS at once
matrix<double> B(n, p, size_tag);  // p right-hand sides
auto X = A.solve(B);  // More efficient than p separate solves
```

---

## Conclusion

The Linear Algebra module provides a comprehensive, efficient, and type-safe implementation of vectors, matrices, and linear algebra operations. Key highlights:

- **Generic Design**: Works with integers, floats, complex numbers, modular arithmetic
- **Performance**: Static extent for compile-time optimization, dynamic extent for flexibility
- **Complete API**: All essential operations from basic arithmetic to advanced decompositions
- **Numerical Algorithms**: Gaussian elimination, Cholesky, QR, SVD, and more
- **Competitive Programming**: Modular arithmetic support, exact computations
- **Well-Tested**: Comprehensive test suite with multiple numeric types

The module is designed for both educational purposes and practical competitive programming applications, providing a solid foundation for algorithm implementations requiring linear algebra.

For more examples and detailed usage, see the test files in `/workspace/tests/linear_algebra/`.
