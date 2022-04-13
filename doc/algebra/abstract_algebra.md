# Abstract Algebra

```c++
#include "algebra/abstract_algebra.h"
```

## 1. Functions

This is a sub-module of algebra. It contains functions that works for many algebraic structures:

### a. `pow`

#### Details

Calculates the 'multiplicative' power of an element `a`

#### Requirements

- `a` is an element of a defined multiplicative monoid $\mathcal{M}$
- `a` is an element of a multiplicative monoid that can be well-founded by additional parameters

### b. `gcd`

#### Details

Calculates the greatest common divisor of two elements `a,b` belonging to an euclidean domain $\mathtt{E}$

#### Requirements

- `a,b` are both elements of the same euclidean domain $\mathtt{E}$
- $\mathtt{E}$ has a total pre-order $<$ compatible with the divisibility relation
- $\mathtt{E}$ supports eulidean division with an operator $/$

### c. `egcd`

#### Details

Calculates the extended eulidean algorithm of two elements `a,b` belonging to an euclidean domain $\mathtt{E}$

#### Requirements

- `a,b` are both elements of the same euclidean domain $\mathtt{E}$
- $\mathtt{E}$ has a total pre-order $<$ compatible with the divisibility relation
- $\mathtt{E}$ supports eulidean division with an operator $/$

### d. `bezout`

#### Details

Calculates the bezout coefficients of two elements `a,b` belonging to an euclidean domain $\mathtt{E}$

#### Requirements

- `a,b` are both elements of the same euclidean domain $\mathtt{E}$
- $\mathtt{E}$ has a total pre-order $<$ compatible with the divisibility relation
- $\mathtt{E}$ supports eulidean division with an operator $/$

## 2. Concepts

This module defines the following concepts

* Monoid
* Group
* Ring
* Euclidean Domain
* Field

