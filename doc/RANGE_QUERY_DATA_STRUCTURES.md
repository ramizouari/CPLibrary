# Range Query Data Structures Documentation

This document provides comprehensive documentation for the range query data structures implemented in CPLibrary. These data structures are optimized for efficiently answering queries on ranges/intervals of arrays and support various binary operations.

## Table of Contents

- [Overview](#overview)
- [Segment Tree](#segment-tree)
- [Fenwick Tree (Binary Indexed Tree)](#fenwick-tree-binary-indexed-tree)
- [Sparse Table](#sparse-table)
- [Prefix Array](#prefix-array)
- [Li-Chao Tree](#li-chao-tree)
- [Comparison Table](#comparison-table)
- [Usage Patterns](#usage-patterns)
- [Advanced Topics](#advanced-topics)

## Overview

### What are Range Query Data Structures?

Range query data structures efficiently solve problems of the form:
- **Query**: Compute an aggregate function (sum, min, max, gcd, etc.) over a range `[l, r]` of an array
- **Update**: Modify elements in the array

These structures optimize the trade-off between query time and update time compared to naive approaches.

### Key Concepts

#### Binary Operations

All range query structures in this library are parameterized by a **binary operation** `O`, which must define:

```cpp
template<typename T>
struct binary_operation {
    using type = T;

    // The operation itself
    T reduce(const T& a, const T& b) const;

    // Neutral element (identity)
    static T neutral;
};
```

**Requirements:**
- **Associativity**: `F(a, F(b, c)) = F(F(a, b), c)`
- **Neutral element**: `F(a, neutral) = F(neutral, a) = a`

#### Common Operations

```cpp
// Addition (sum queries)
template<typename T>
struct plus_t : binary_operation<T> {
    T reduce(const T& a, const T& b) const { return a + b; }
    static inline T neutral = T(0);
};

// Minimum (range minimum query)
template<typename T>
struct min_t : binary_operation<T> {
    T reduce(const T& a, const T& b) const { return std::min(a, b); }
    static inline T neutral = std::numeric_limits<T>::max();
};

// Maximum (range maximum query)
template<typename T>
struct max_t : binary_operation<T> {
    T reduce(const T& a, const T& b) const { return std::max(a, b); }
    static inline T neutral = std::numeric_limits<T>::lowest();
};

// Multiplication
template<typename T>
struct multiplies_t : binary_operation<T> {
    T reduce(const T& a, const T& b) const { return a * b; }
    static inline T neutral = T(1);
};

// XOR
template<typename T>
struct xor_t : binary_operation<T> {
    T reduce(const T& a, const T& b) const { return a ^ b; }
    static inline T neutral = T(0);
};

// GCD
template<typename T>
struct gcd_t : binary_operation<T> {
    T reduce(const T& a, const T& b) const { return std::gcd(a, b); }
    static inline T neutral = T(0);
};
```

#### Invertible Operations

Some data structures (like Fenwick Tree) require **invertible operations** - operations with an inverse:

```cpp
template<typename T>
struct invertible_operation {
    // Compute inverse of a
    T inv(const T& a) const;
};

// Additive inverse: -a
template<typename T>
struct additive_inverse_t : invertible_operation<T> {
    T inv(const T& a) const { return -a; }
};

// XOR is self-inverse: a ^ a = 0
template<typename T>
struct involution_inverse_t : invertible_operation<T> {
    T inv(const T& a) const { return a; }
};
```

### Structure Variants

Each data structure comes in two variants:

1. **Fixed** (`cp::data_structures::fixed`): Static, compile-time parameterized binary operation
2. **Dynamic** (`cp::data_structures::dynamic`): Runtime-selectable binary operation via polymorphism

---

## Segment Tree

### Purpose

Segment trees answer range queries over **any associative binary operation** and support point updates efficiently.

### Implementation Details

**Location**: `include/data_structures/fixed/segment_tree.h`

#### Data Structure

```cpp
template<typename O>
struct segment_tree {
    using R = typename O::type;
    using type = R;

    std::vector<std::vector<R>> S;  // Tree levels
    std::vector<R> A;                // Original array
    int n, h;                        // Size and height

    // Construction, query, and update methods
};
```

#### Internal Representation

The segment tree stores data in a hierarchical structure:
- `S[0]`: Root level (1 element - aggregate of entire array)
- `S[h-1]`: Leaf level (n elements - original array values)
- `S[i]`: Level i contains `2^i` nodes

**Tree Structure:**
```
Level 0:              [sum(0-7)]
                      /        \
Level 1:      [sum(0-3)]      [sum(4-7)]
              /      \          /      \
Level 2:  [sum(0-1)] [sum(2-3)] [sum(4-5)] [sum(6-7)]
          /   \      /   \      /   \      /   \
Level 3: [a0] [a1] [a2] [a3] [a4] [a5] [a6] [a7]
```

#### Algorithm: Build

```cpp
void build() {
    // Copy elements to leaf level
    for(int i = 0; i < n; i++)
        S[h-1][i] = A[i];

    // Build tree bottom-up
    for(int level = h-2; level >= 0; level--) {
        for(int k = 0; k < (1 << level); k++) {
            S[level][k] = F(S[level+1][2*k], S[level+1][2*k+1]);
        }
    }
}
```

**Time Complexity**: O(n)

#### Algorithm: Point Update

```cpp
void update(int i, R u) {
    A[i] = u;
    S[h-1][i] = u;  // Update leaf

    int m = h-2;
    i /= 2;

    // Propagate upward
    while(m >= 0) {
        S[m][i] = F(S[m+1][2*i], S[m+1][2*i+1]);
        m--;
        i /= 2;
    }
}
```

**Time Complexity**: O(log n)

#### Algorithm: Range Query

```cpp
R query(int l, int r, int a, int b, int depth) {
    if(l >= r) return O::neutral;

    // Exact match at this level
    if(l == a && r == b)
        return S[depth][l >> (h-1-depth)];

    int mid = (a + b) / 2;

    // Completely in left or right subtree
    if(mid > r)
        return query(l, r, a, mid, depth+1);
    else if(mid < l)
        return query(l, r, mid, b, depth+1);

    // Split query
    else
        return F(query(l, mid, a, mid, depth+1),
                query(mid, r, mid, b, depth+1));
}
```

**Time Complexity**: O(log n)

### Time and Space Complexity

| Operation | Time | Space |
|-----------|------|-------|
| Construction | O(n) | O(n) |
| Point Update | O(log n) | - |
| Range Query | O(log n) | - |
| Memory | - | O(2n) = O(n) |

### API Reference

#### Fixed Variant

```cpp
namespace cp::data_structures::fixed {

template<typename O>
struct segment_tree {
    using R = typename O::type;

    // Constructor: Build from array
    segment_tree(const std::vector<R>& _A);

    // Update element at index i to value u
    void update(int i, R u);

    // Query range [l, r) with bounds checking
    R query(int l, int r);

private:
    std::vector<std::vector<R>> S;  // Tree levels
    std::vector<R> A;                // Array
    int n, h;                        // Size, height
    inline static O F = O();         // Operation instance
};

}
```

#### Dynamic Variant

```cpp
namespace cp::data_structures::dynamic {

template<typename R>
struct segment_tree {
    // Constructor: Build from array and operation
    segment_tree(const std::vector<R>& _A,
                 std::shared_ptr<binary_operation<R>> _F);

    // Same interface as fixed variant
    void update(int i, R u);
    R query(int l, int r);

private:
    std::vector<std::vector<R>> S;
    std::vector<R> A;
    int n, h;
    binary_operation_ptr<R> F;  // Runtime polymorphic operation
};

}
```

### Usage Examples

#### Example 1: Range Sum Queries

```cpp
#include "data_structures/range_queries.h"
#include <iostream>
#include <vector>

int main() {
    using namespace cp::data_structures::fixed;

    std::vector<int> data = {1, 3, 5, 7, 9, 11, 13, 15};

    // Create segment tree for sum queries
    segment_tree<plus_t<int>> seg(data);

    // Query sum of range [2, 6) -> 5+7+9+11 = 32
    int sum = seg.query(2, 6);
    std::cout << "Sum [2, 6): " << sum << std::endl;  // Output: 32

    // Update element at index 3 to 10
    seg.update(3, 10);

    // Query again [2, 6) -> 5+10+9+11 = 35
    sum = seg.query(2, 6);
    std::cout << "Sum after update: " << sum << std::endl;  // Output: 35

    return 0;
}
```

#### Example 2: Range Minimum Queries

```cpp
#include "data_structures/range_queries.h"
#include <iostream>
#include <vector>
#include <limits>

int main() {
    using namespace cp::data_structures::fixed;

    // Initialize neutral element for min
    min_t<int>::neutral = std::numeric_limits<int>::max();

    std::vector<int> data = {5, 2, 8, 1, 9, 3, 7, 4};
    segment_tree<min_t<int>> seg(data);

    // Query minimum in range [1, 5) -> min(2, 8, 1, 9) = 1
    int minimum = seg.query(1, 5);
    std::cout << "Min [1, 5): " << minimum << std::endl;  // Output: 1

    // Update: change element at index 3 from 1 to 15
    seg.update(3, 15);

    // Query again [1, 5) -> min(2, 8, 15, 9) = 2
    minimum = seg.query(1, 5);
    std::cout << "Min after update: " << minimum << std::endl;  // Output: 2

    return 0;
}
```

#### Example 3: Range GCD Queries

```cpp
#include "data_structures/range_queries.h"
#include <iostream>
#include <vector>
#include <numeric>

// Define GCD operation
template<typename T>
struct gcd_t {
    using type = T;
    T reduce(const T& a, const T& b) const {
        return std::gcd(a, b);
    }
    static inline T neutral = T(0);
};

int main() {
    using namespace cp::data_structures::fixed;

    std::vector<int> data = {12, 18, 24, 36, 48};
    segment_tree<gcd_t<int>> seg(data);

    // Query GCD of range [1, 4) -> gcd(18, 24, 36) = 6
    int result = seg.query(1, 4);
    std::cout << "GCD [1, 4): " << result << std::endl;  // Output: 6

    return 0;
}
```

#### Example 4: Matrix Range Queries

```cpp
#include "data_structures/range_queries.h"
#include "linear_algebra/matrix.h"
#include <iostream>

int main() {
    using namespace cp;
    using namespace cp::data_structures::fixed;

    using Matrix = s_matrix<int, 2, 2>;

    // Create array of matrices
    std::vector<Matrix> matrices;
    matrices.push_back({{1, 0}, {0, 1}});  // Identity
    matrices.push_back({{2, 1}, {1, 1}});
    matrices.push_back({{1, 1}, {0, 1}});

    // Segment tree for matrix multiplication
    segment_tree<multiplies_t<Matrix>> seg(matrices);

    // Query product of matrices [0, 2)
    Matrix result = seg.query(0, 2);
    std::cout << "Product [0, 2):" << std::endl;
    std::cout << result(0,0) << " " << result(0,1) << std::endl;
    std::cout << result(1,0) << " " << result(1,1) << std::endl;

    return 0;
}
```

### Special Features

1. **Generic over Binary Operations**: Works with any associative operation
2. **Hierarchical Storage**: Level-based storage for efficient traversal
3. **Automatic Size Adjustment**: Pads array to power of 2 internally
4. **Type Safety**: Compile-time type checking with templates

---

## Fenwick Tree (Binary Indexed Tree)

### Purpose

Fenwick trees (Binary Indexed Trees) efficiently compute prefix sums and support updates for **group operations** (operations with inverses). They are more space-efficient than segment trees for these specific operations.

### Implementation Details

**Location**: `include/data_structures/fixed/fenwick_tree.h`

#### Data Structure

```cpp
template<typename O>
struct fenwick_tree {
    using T = typename O::type;
    int n;
    std::vector<T> bit;  // Binary indexed tree array
    inline static O F = O();

    // Construction, query, and update methods
};
```

#### Internal Representation

Fenwick tree stores partial sums in an array where `bit[i]` contains the sum of a specific range determined by the binary representation of `i`.

**Key Insight**: Each index `i` is responsible for elements in range `[i - LSB(i) + 1, i]`, where `LSB(i)` is the least significant bit.

**Example for n=8:**
```
Index:    1    2    3    4    5    6    7    8
bit[i]:  a[0] a[0] a[2] a[0] a[4] a[4] a[6] a[0]
         +a[1]      +a[1]      +a[5]      +a[1]
                    +a[2]               +a[2]
                    +a[3]               +a[3]
                                        +a[4]
                                        +a[5]
                                        +a[6]
                                        +a[7]
```

#### Algorithm: Point Update (Add)

```cpp
void add(int x, T delta) {
    // Propagate update upward through tree
    for (int i = x; i < n; i = i | (i + 1))
        bit[i] = F(bit[i], delta);
}
```

**Time Complexity**: O(log n)

**How it works**: The expression `i | (i + 1)` adds the least significant bit, moving to the next node that depends on index `i`.

#### Algorithm: Prefix Sum

```cpp
T sum(int x) {
    if(x < 0) return O::neutral;

    T ret = O::neutral;
    // Traverse down the tree
    for (int i = x; i >= 0; i = (i & (i + 1)) - 1)
        ret = F(ret, bit[i]);
    return ret;
}
```

**Time Complexity**: O(log n)

**How it works**: The expression `(i & (i + 1)) - 1` removes the least significant bit, moving to the previous independent node.

#### Algorithm: Range Query

```cpp
T query(int a, int b) {
    // sum[a, b] = sum[0, b] - sum[0, a-1]
    return F(F.inv(sum(a-1)), sum(b));
}
```

**Time Complexity**: O(log n)

**Requirements**: Operation must be invertible!

### Time and Space Complexity

| Operation | Time | Space |
|-----------|------|-------|
| Construction | O(n log n) | O(n) |
| Point Update (add) | O(log n) | - |
| Point Update (set) | O(log n) | - |
| Prefix Query | O(log n) | - |
| Range Query | O(log n) | - |
| Memory | - | O(n) |

### API Reference

```cpp
namespace cp::data_structures::fixed {

template<typename O>
struct fenwick_tree {
    using T = typename O::type;

    // Constructor: Initialize with size
    fenwick_tree(int _n);

    // Constructor: Build from array
    fenwick_tree(const std::vector<T>& X);

    // Compute prefix sum [0, x]
    T sum(int x);

    // Compute range sum [a, b]
    T query(int a, int b);

    // Alias for query
    T sum(int a, int b);

    // Add delta to element at position x
    void add(int x, T delta);

    // Set element at position x to delta
    void update(int x, T delta);

private:
    int n;
    std::vector<T> bit;
    inline static O F = O();
};

}
```

### Usage Examples

#### Example 1: Basic Range Sum

```cpp
#include "data_structures/range_queries.h"
#include <iostream>
#include <vector>

int main() {
    using namespace cp::data_structures::fixed;

    int n = 10;
    fenwick_tree<plus_t<int>> fen(n);

    // Build array: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    for (int i = 0; i < n; i++) {
        fen.add(i, i + 1);
    }

    // Query range sum [2, 7] -> 3+4+5+6+7+8 = 33
    int sum = fen.query(2, 7);
    std::cout << "Sum [2, 7]: " << sum << std::endl;  // Output: 33

    // Add 5 to element at index 4
    fen.add(4, 5);  // Now array[4] = 5 + 5 = 10

    // Query again [2, 7] -> 3+4+10+6+7+8 = 38
    sum = fen.query(2, 7);
    std::cout << "Sum after update: " << sum << std::endl;  // Output: 38

    return 0;
}
```

#### Example 2: Prefix Sums

```cpp
#include "data_structures/range_queries.h"
#include <iostream>
#include <vector>

int main() {
    using namespace cp::data_structures::fixed;

    std::vector<int> data = {1, 2, 3, 4, 5};
    fenwick_tree<plus_t<int>> fen(data);

    // Print all prefix sums
    for (int i = 0; i < data.size(); i++) {
        std::cout << "sum[0, " << i << "] = "
                  << fen.sum(i) << std::endl;
    }
    // Output:
    // sum[0, 0] = 1
    // sum[0, 1] = 3
    // sum[0, 2] = 6
    // sum[0, 3] = 10
    // sum[0, 4] = 15

    return 0;
}
```

#### Example 3: XOR Queries

```cpp
#include "data_structures/range_queries.h"
#include <iostream>
#include <vector>

int main() {
    using namespace cp::data_structures::fixed;

    std::vector<int> data = {1, 3, 5, 7, 9};
    fenwick_tree<xor_t<int>> fen(data);

    // Query XOR of range [1, 4]
    // XOR(3, 5, 7, 9) = 3 ^ 5 ^ 7 ^ 9 = 12
    int result = fen.query(1, 4);
    std::cout << "XOR [1, 4]: " << result << std::endl;

    // XOR is self-inverse: a ^ a = 0
    // So Fenwick tree works perfectly for XOR

    return 0;
}
```

#### Example 4: Frequency Counting

```cpp
#include "data_structures/range_queries.h"
#include <iostream>

int main() {
    using namespace cp::data_structures::fixed;

    // Track frequency of elements
    const int MAX_VAL = 100;
    fenwick_tree<plus_t<int>> freq(MAX_VAL + 1);

    // Insert elements
    std::vector<int> elements = {5, 10, 5, 20, 10, 5, 15};
    for (int x : elements) {
        freq.add(x, 1);  // Increment frequency
    }

    // Count elements in range [0, 15]
    int count = freq.sum(15);
    std::cout << "Count of elements <= 15: " << count << std::endl;

    // Count elements in range [5, 15]
    count = freq.query(5, 15);
    std::cout << "Count in [5, 15]: " << count << std::endl;

    return 0;
}
```

### Special Features

1. **Space Efficient**: Uses exactly n elements (vs 2n for segment tree)
2. **Cache Friendly**: Better cache locality than segment tree
3. **Simple Implementation**: Elegant bit manipulation
4. **Invertible Operations Only**: Requires group structure (with inverse)

### When to Use Fenwick Tree vs Segment Tree

**Use Fenwick Tree when:**
- Operation is invertible (sum, XOR, etc.)
- Need space efficiency
- Only point updates needed
- Simpler implementation preferred

**Use Segment Tree when:**
- Operation is not invertible (min, max, GCD)
- Need range updates (lazy propagation)
- More flexibility needed

---

## Sparse Table

### Purpose

Sparse tables answer **immutable** range queries (no updates) for **idempotent** operations in O(1) time per query after O(n log n) preprocessing.

### Implementation Details

**Location**: `include/data_structures/fixed/sparse_array.h`

#### Data Structure

```cpp
template<typename O>
struct sparse_array {
    using T = typename O::type;
    using type = T;

    int n, h;
    std::vector<std::vector<T>> S;  // Sparse table

    // Construction and query methods
};
```

#### Internal Representation

Sparse table stores precomputed answers for all ranges of length 2^k:
- `S[k][i]` = result of operation on range `[i, i + 2^k)`

**Example for n=8:**
```
S[0][i]: [a0] [a1] [a2] [a3] [a4] [a5] [a6] [a7]  (2^0 = 1 element)
S[1][i]: [a0,a1] [a1,a2] [a2,a3] ...              (2^1 = 2 elements)
S[2][i]: [a0-a3] [a1-a4] [a2-a5] ...              (2^2 = 4 elements)
S[3][i]: [a0-a7] [a1-a8]                          (2^3 = 8 elements)
```

#### Algorithm: Build

```cpp
sparse_array(const std::vector<T>& A) : n(bit_ceil(A.size())),
                                         h(bit_log(n)),
                                         S(h+1) {
    // Allocate levels
    int r = 1;
    for(int i = h; i >= 0; i--, r *= 2)
        S[i].resize(n - r + 1, O::neutral);

    // Copy base values
    for(int i = 0; i < A.size(); i++)
        S[h][i] = A[i];

    // Build table bottom-up
    r = 1;
    for(int i = h-1; i >= 0; i--, r *= 2) {
        for(int j = 0; j <= n - 2*r; j++)
            S[i][j] = F(S[i+1][j], S[i+1][j+r]);
    }
}
```

**Time Complexity**: O(n log n)

#### Algorithm: Range Query

```cpp
T query(int l, int r) const {
    if(l >= r) return O::neutral;

    int d = r - l;             // Range length
    int s = bit_floor(d);      // Largest power of 2 <= d
    int b = bit_log(s);        // log2(s)

    // Combine two overlapping ranges
    return F(S[h-b][l], S[h-b][r-s]);
}
```

**Time Complexity**: O(1)

**Key Insight**: For idempotent operations (where `F(x, x) = x`), we can use two overlapping ranges of length `2^k` that cover `[l, r]`.

### Time and Space Complexity

| Operation | Time | Space |
|-----------|------|-------|
| Construction | O(n log n) | O(n log n) |
| Range Query | O(1) | - |
| Update | Not supported | - |
| Memory | - | O(n log n) |

### API Reference

```cpp
namespace cp::data_structures::fixed {

template<typename O>
struct sparse_array {
    using T = typename O::type;

    // Constructor: Build from array
    sparse_array(const std::vector<T>& A);

    // Query range [l, r) in O(1)
    T query(int l, int r) const;

private:
    int n, h;
    std::vector<std::vector<T>> S;
    inline static O F = O();
};

}
```

### Usage Examples

#### Example 1: Range Minimum Query (RMQ)

```cpp
#include "data_structures/range_queries.h"
#include <iostream>
#include <vector>
#include <limits>

int main() {
    using namespace cp::data_structures::fixed;

    // Initialize neutral element
    min_t<int>::neutral = std::numeric_limits<int>::max();

    std::vector<int> data = {7, 2, 3, 0, 5, 10, 3, 12, 18};
    sparse_array<min_t<int>> st(data);

    // O(1) range minimum queries
    std::cout << "Min [0, 4]: " << st.query(0, 4) << std::endl;  // 0
    std::cout << "Min [2, 7]: " << st.query(2, 7) << std::endl;  // 0
    std::cout << "Min [5, 9]: " << st.query(5, 9) << std::endl;  // 3

    // Note: Updates not supported!

    return 0;
}
```

#### Example 2: Range GCD Query

```cpp
#include "data_structures/range_queries.h"
#include <iostream>
#include <vector>
#include <numeric>

// GCD is idempotent: gcd(x, x) = x
template<typename T>
struct gcd_t {
    using type = T;
    T reduce(const T& a, const T& b) const {
        return std::gcd(a, b);
    }
    static inline T neutral = T(0);
};

int main() {
    using namespace cp::data_structures::fixed;

    std::vector<int> data = {12, 18, 24, 36, 48, 60};
    sparse_array<gcd_t<int>> st(data);

    // O(1) GCD queries
    std::cout << "GCD [0, 3]: " << st.query(0, 3) << std::endl;  // 6
    std::cout << "GCD [2, 5]: " << st.query(2, 5) << std::endl;  // 12
    std::cout << "GCD [0, 6]: " << st.query(0, 6) << std::endl;  // 6

    return 0;
}
```

#### Example 3: Range Maximum Query

```cpp
#include "data_structures/range_queries.h"
#include <iostream>
#include <vector>
#include <limits>

int main() {
    using namespace cp::data_structures::fixed;

    max_t<int>::neutral = std::numeric_limits<int>::lowest();

    std::vector<int> data = {3, 1, 4, 1, 5, 9, 2, 6, 5};
    sparse_array<max_t<int>> st(data);

    // O(1) range maximum queries
    std::cout << "Max [0, 5]: " << st.query(0, 5) << std::endl;  // 5
    std::cout << "Max [3, 8]: " << st.query(3, 8) << std::endl;  // 9
    std::cout << "Max [2, 6]: " << st.query(2, 6) << std::endl;  // 9

    return 0;
}
```

### Special Features

1. **O(1) Query Time**: Fastest query among all structures
2. **Idempotent Operations**: Works for min, max, GCD, LCM, AND, OR
3. **Static Data**: No update support (immutable)
4. **Overlapping Ranges**: Uses idempotency to overlap ranges for O(1) queries

### Idempotency Requirement

**Idempotent**: `F(x, x) = x`

**Idempotent operations:**
- min, max (idempotent)
- GCD, LCM (idempotent)
- Bitwise AND, OR (idempotent)

**Non-idempotent operations:**
- Sum: sum(x, x) = 2x ≠ x
- Product: product(x, x) = x² ≠ x
- XOR: xor(x, x) = 0 ≠ x

---

## Prefix Array

### Purpose

Prefix arrays compute prefix sums for **group operations** and answer range queries. Updates are expensive (O(n)), making this structure useful when queries dominate or data is mostly static.

### Implementation Details

**Location**: `include/data_structures/fixed/prefix_array.h`

#### Data Structure

```cpp
template<typename O>
struct prefix_array {
    using R = typename O::type;
    using type = typename O::type;

    std::vector<R> A;  // Original array
    std::vector<R> P;  // Prefix sums

    // Construction, query, and update methods
};
```

#### Internal Representation

```
Array A:    [a0, a1, a2, a3, a4]
Prefix P:   [e,  a0, a0+a1, a0+a1+a2, a0+a1+a2+a3, a0+a1+a2+a3+a4]
             ^
             neutral element
```

Where `P[i] = F(a[0], a[1], ..., a[i-1])`

#### Algorithm: Build

```cpp
prefix_array(const std::vector<R>& _A) : A(_A), P(_A.size()+1) {
    P[0] = O::neutral;
    for(int i = 0; i < A.size(); i++)
        P[i+1] = F(P[i], A[i]);
}
```

**Time Complexity**: O(n)

#### Algorithm: Range Query

```cpp
R query(int l, int r) {
    // sum[l, r) = P[r] - P[l]
    return F(F.inv(P[l]), P[r]);
}
```

**Time Complexity**: O(1)

#### Algorithm: Update

```cpp
void update(int i, R u) {
    A[i] = u;
    // Rebuild all prefixes from i onward
    for(int j = i+1; j < P.size(); j++)
        P[j] = F(P[j-1], A[j-1]);
}
```

**Time Complexity**: O(n)

### Time and Space Complexity

| Operation | Time | Space |
|-----------|------|-------|
| Construction | O(n) | O(n) |
| Range Query | O(1) | - |
| Update | O(n) | - |
| Memory | - | O(n) |

### API Reference

```cpp
namespace cp::data_structures::fixed {

template<typename O>
struct prefix_array {
    using R = typename O::type;

    // Constructor: Build from array
    prefix_array(const std::vector<R>& _A);

    // Query range [l, r) in O(1)
    R query(int l, int r);

    // Update element at index i (slow: O(n))
    void update(int i, R u);

private:
    std::vector<R> A;  // Original array
    std::vector<R> P;  // Prefix array
    inline static O F = O();
};

}
```

### Usage Examples

#### Example 1: Range Sum with Rare Updates

```cpp
#include "data_structures/range_queries.h"
#include <iostream>
#include <vector>

int main() {
    using namespace cp::data_structures::fixed;

    std::vector<int> data = {1, 2, 3, 4, 5};
    prefix_array<plus_t<int>> prefix(data);

    // Many O(1) queries
    std::cout << "Sum [0, 3): " << prefix.query(0, 3) << std::endl;  // 6
    std::cout << "Sum [2, 5): " << prefix.query(2, 5) << std::endl;  // 12
    std::cout << "Sum [1, 4): " << prefix.query(1, 4) << std::endl;  // 9

    // Rare O(n) update
    prefix.update(2, 10);  // Change A[2] from 3 to 10

    std::cout << "Sum [2, 5): " << prefix.query(2, 5) << std::endl;  // 19

    return 0;
}
```

### When to Use

**Use Prefix Array when:**
- Queries vastly outnumber updates
- Updates are very rare or data is static after initialization
- Want simplest implementation
- O(1) query time is critical

**Don't use when:**
- Frequent updates (use Fenwick or Segment tree)
- Need better than O(n) update time

---

## Li-Chao Tree

### Purpose

Li-Chao trees maintain a set of linear functions and answer queries for the minimum (or maximum) value at a given point among all functions. Useful for convex hull trick and dynamic programming optimizations.

### Implementation Details

**Location**: `include/data_structures/li_chao.h`

#### Data Structure

```cpp
template<std::integral I, typename T>
struct li_chao_tree {
    using fn_ptr = std::shared_ptr<function<I,T>>;

    std::vector<std::vector<fn_ptr>> S;  // Tree levels
    integer l;
    natural n, h;

    // Methods for adding functions and querying
};
```

**Note**: This structure is currently under development (incomplete implementation in the codebase).

#### Concept

Li-Chao tree stores functions (typically linear) at tree nodes and efficiently determines which function gives the minimum/maximum value at any query point.

**Key Operations:**
1. **Add function**: Insert a linear function `f(x) = ax + b`
2. **Query point**: Find minimum value at point `x` among all functions

#### Algorithm Overview

For each node covering range `[l, r]`:
1. Store the "dominant" function at the midpoint
2. Recursively store subdominant functions in child nodes

**Time Complexity:**
- Add function: O(log n)
- Query point: O(log n)

### Usage Context

**Common Applications:**
- Convex Hull Trick optimization
- Dynamic Programming with slope optimization
- Maintaining lower/upper envelope of lines
- Competitive programming DP optimizations

**Example Use Case:**
```
Problem: Given n lines, find minimum value at various x-coordinates

Solution:
1. Add all lines to Li-Chao tree
2. Query minimum at each x in O(log n)
```

### API Reference (Planned)

```cpp
namespace cp::data_structures {

template<std::integral I, typename T>
struct li_chao_tree {
    // Constructor: Initialize with coordinate range
    li_chao_tree(I l, I r);

    // Add a function to the tree
    void add_function(fn_ptr fn);

    // Query minimum value at point x
    T minimum(I x);

private:
    std::vector<std::vector<fn_ptr>> S;
    integer l;
    natural n, h;
};

}
```

---

## Comparison Table

### Feature Comparison

| Data Structure | Update | Query | Space | Operations | Notes |
|----------------|--------|-------|-------|------------|-------|
| **Segment Tree** | O(log n) | O(log n) | O(n) | Associative | Most versatile |
| **Fenwick Tree** | O(log n) | O(log n) | O(n) | Group (invertible) | Space efficient |
| **Sparse Table** | No updates | O(1) | O(n log n) | Idempotent | Fastest queries |
| **Prefix Array** | O(n) | O(1) | O(n) | Group | Simple, query-heavy |
| **Li-Chao Tree** | O(log n) | O(log n) | O(n) | Min/Max of functions | Specialized |

### Operation Support

| Operation | Segment Tree | Fenwick Tree | Sparse Table | Prefix Array |
|-----------|--------------|--------------|--------------|--------------|
| Sum | ✅ | ✅ | ❌* | ✅ |
| Min/Max | ✅ | ❌ | ✅ | ❌ |
| GCD | ✅ | ❌ | ✅ | ❌ |
| XOR | ✅ | ✅ | ❌* | ✅ |
| Product | ✅ | ❌ | ❌ | ❌ |
| Custom Associative | ✅ | ❌ | ✅** | ❌ |

*Not idempotent
**If idempotent

### Decision Guide

```
Need updates?
├─ Yes
│  ├─ Operation invertible? (has inverse)
│  │  ├─ Yes → Fenwick Tree (space efficient)
│  │  └─ No → Segment Tree
│  └─ Operation not invertible → Segment Tree
│
└─ No (static data)
   ├─ Operation idempotent?
   │  ├─ Yes → Sparse Table (O(1) query)
   │  └─ No → Prefix Array or Segment Tree
   └─ Rare updates → Prefix Array (simplest)
```

---

## Usage Patterns

### Pattern 1: Choosing the Right Operation Type

```cpp
// For sum queries (invertible)
using SumTree = fenwick_tree<plus_t<int>>;

// For min/max queries (idempotent, no updates)
using RMQ = sparse_array<min_t<int>>;

// For GCD queries (idempotent with updates)
using GCDTree = segment_tree<gcd_t<int>>;

// For product queries (not invertible)
using ProdTree = segment_tree<multiplies_t<int>>;
```

### Pattern 2: Custom Operations

```cpp
// Define custom operation
template<typename T>
struct custom_op {
    using type = T;

    T reduce(const T& a, const T& b) const {
        // Custom logic
        return /* ... */;
    }

    static inline T neutral = /* ... */;
};

// Use with segment tree
segment_tree<custom_op<MyType>> tree(data);
```

### Pattern 3: Dynamic vs Fixed

```cpp
// Fixed: Compile-time operation
fixed::segment_tree<plus_t<int>> fixed_tree(data);

// Dynamic: Runtime operation
auto op = std::make_shared<monoid_plus_t<int>>();
dynamic::segment_tree<int> dynamic_tree(data, op);

// Fixed is faster (inline, no virtual calls)
// Dynamic is more flexible (can change operation)
```

### Pattern 4: Modular Arithmetic

```cpp
#include "nt/modular_arithmetic.h"

constexpr int MOD = 1000000007;
using mint = cyclic<MOD>;

// Sum queries with modular arithmetic
std::vector<mint> data = {1, 2, 3, 4, 5};
segment_tree<plus_t<mint>> tree(data);

mint sum = tree.query(0, 5);  // Automatically mod MOD
```

### Pattern 5: Matrix Queries

```cpp
#include "linear_algebra/matrix.h"

using Matrix = s_matrix<int, 2, 2>;

// Matrix multiplication queries
std::vector<Matrix> matrices = /* ... */;
segment_tree<multiplies_t<Matrix>> tree(matrices);

// Query product of matrices [l, r)
Matrix product = tree.query(l, r);
```

---

## Advanced Topics

### 1. Lazy Propagation (Segment Tree)

For **range updates** efficiently, segment trees can use lazy propagation:

**Concept**: Delay updates until necessary
- Store pending updates at nodes
- Propagate lazily when visiting nodes

**Time Complexity**:
- Range update: O(log n)
- Range query: O(log n)

**Note**: Not yet implemented in current codebase, but standard extension.

### 2. Persistent Data Structures

**Concept**: Maintain history of all versions
- Each update creates a new version
- Old versions remain accessible
- Uses path copying

**Space**: O(log n) per update
**Use case**: Time-travel queries, functional programming

### 3. 2D Range Queries

**Concept**: Extend to 2D arrays
- 2D Segment Tree: Query rectangles
- 2D Fenwick Tree: Query rectangles (if invertible)

**Time Complexity**: O(log² n) queries and updates

### 4. Range Update Range Query (RURQ)

**Segment Tree with Lazy Propagation:**
- Update range: add value to all elements in [l, r]
- Query range: get sum/min/max of [l, r]

**Fenwick Tree with Difference Array:**
- Use difference array technique
- Two Fenwick trees for range updates

### 5. Dynamic Segment Tree (Coordinate Compression)

**Concept**: Allocate nodes on-demand
- Handle large coordinate ranges (up to 10⁹)
- Only create nodes as needed

**Use case**: Sparse data, large coordinate space

### 6. Heavy-Light Decomposition (Tree Queries)

**Concept**: Decompose tree into paths
- Use segment trees on paths
- Answer path queries on trees

**Location**: `include/graph/tree/range_queries.h`
**Documented separately**: See tree algorithms documentation

---

## Best Practices

### 1. Initialize Neutral Elements

```cpp
// For min queries
min_t<int>::neutral = std::numeric_limits<int>::max();

// For max queries
max_t<int>::neutral = std::numeric_limits<int>::lowest();
```

### 2. Handle Empty Ranges

```cpp
// Always check before querying
if (l < r) {
    result = tree.query(l, r);
} else {
    result = O::neutral;  // Empty range
}
```

### 3. Choose Appropriate Structure

```cpp
// Many updates → Fenwick or Segment tree
// No updates → Sparse table (if idempotent)
// Rare updates, many queries → Prefix array
```

### 4. Type Safety

```cpp
// Use strong types to prevent errors
using Sum = segment_tree<plus_t<int>>;
using Min = segment_tree<min_t<int>>;

Sum sum_tree(data);
Min min_tree(data);

// Can't accidentally mix operations
```

### 5. Testing

```cpp
// Test against naive implementation
for (int i = 0; i < Q; i++) {
    int l = random(0, n);
    int r = random(l, n);

    // Fast query
    auto fast_result = tree.query(l, r);

    // Naive verification
    auto naive_result = naive_query(data, l, r);

    assert(fast_result == naive_result);
}
```

---

## Performance Considerations

### Cache Locality

**Segment Tree**: Level-based storage improves cache performance
**Fenwick Tree**: Linear array, excellent cache locality
**Sparse Table**: Larger memory footprint, more cache misses

### Memory Usage

```
For n = 10⁶:
- Segment Tree: ~8 MB (2n integers)
- Fenwick Tree: ~4 MB (n integers)  ← Most space efficient
- Sparse Table: ~80 MB (n log n integers)
```

### Compilation

```bash
# Optimize for performance
g++ -std=c++23 -O3 -march=native -DNDEBUG \
    -I/path/to/CPLibrary/include program.cpp -o program
```

---

## Conclusion

CPLibrary provides a comprehensive suite of range query data structures optimized for competitive programming:

- **Segment Tree**: Most versatile, handles any associative operation
- **Fenwick Tree**: Space-efficient, for invertible operations
- **Sparse Table**: Lightning-fast O(1) queries for immutable idempotent operations
- **Prefix Array**: Simplest, for query-dominant scenarios
- **Li-Chao Tree**: Specialized for function optimization problems

Choose the appropriate structure based on your problem's requirements (update frequency, operation type, space constraints) and leverage the generic design to work with custom types and operations.

For more examples and applications, see:
- [EXAMPLES.md](../EXAMPLES.md) - Comprehensive usage examples
- [API.md](../API.md) - Complete API reference
- [tests/data_structures/](tests/data_structures/) - Unit tests and benchmarks

---

**Last Updated**: December 2025
**Author**: Documentation generated for CPLibrary
**License**: See repository for license information
