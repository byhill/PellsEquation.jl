# PellsEquation

[![Build Status](https://ci.appveyor.com/api/projects/status/github/byhill/PellsEquation.jl?svg=true)](https://ci.appveyor.com/project/byhill/PellsEquation-jl)

## Overview

This package gives provides support for finding solutions to the [generalized](https://en.wikipedia.org/wiki/Pell's_equation#Generalized_Pell's_equation) [Pell's equation](https://en.wikipedia.org/wiki/Pell's_equation).
More specifically, give an integer $N$ and an integer $D \geq 0$,
this package will find all nonnegative integer solutions $(x, y)$ to the diophantine equation
$x^2 - Dy^2 = N$.

To do this, you can call `pellseqn(D, N=1)`.
This will return a (possibly infinite) iterator whose eltype is `Tuple{BigInt,BigInt}`.

## Examples

```julia
julia> Pkg.add("PellsEquation")

julia> using PellsEquation, .Iterators

# Find the first five solutions to x^2 - 2y^2 = 1.
julia> collect(take(pellseqn(2), 5))
5-element Vector{Tuple{BigInt, BigInt}}:
 (1, 0)
 (3, 2)
 (17, 12)
 (99, 70)
 (577, 408)

# Find the first five solutions to x^2 - 2y^2 = 1,
# this time explicitly stating the value of N.
julia> collect(take(pellseqn(2, 1), 5))
5-element Vector{Tuple{BigInt, BigInt}}:
 (1, 0)
 (3, 2)
 (17, 12)
 (99, 70)
 (577, 408)

# Find the first five solutions of x^2 - 5y^2 = 1.
julia> collect(take(pellseqn(5, -1), 5))
5-element Vector{Tuple{BigInt, BigInt}}:
 (2, 1)
 (38, 17)
 (682, 305)
 (12238, 5473)
 (219602, 98209)

# Find the solutions to x^2 - 4y^2 = 25.
# Here, there are only 2 solutions in total,
# so it stops iterating once the two solutions are found.
julia> collect(take(pellseqn(4, 25), 5))
2-element Vector{Tuple{BigInt, BigInt}}:
 (5, 0)
 (13, 6)

# Find the solutions to x^2 - 4y^2 = 10.
# There are no solutions so it returns and empty iterator.
julia> collect(take(pellseqn(4, 10), 5))
Tuple{BigInt, BigInt}[]

# Find the first five solutions to x^2 - 921y^2 = -12.
# Notice how the solutions can grow extremely quickly,
# even with very modest sized inputs.
# This is why we always use BigInts,
# even if the inputs themselves are not BigInts.
julia> collect(take(pellseqn(921, -12), 5))
5-element Vector{Tuple{BigInt, BigInt}}:
 (123013602, 4053436)
 (620494807415610916348542, 20445999054572056617484)
 (3129847429634131057287425640460782483138, 103131979224431202543397867947637928044)
 (15787311699815721226111492312624904230758974206686324318, 520209607285982995690281237198588485553526346187087196)
 (79633022474924155281204153475574155565247711565894067613899610557324322, 2623997304693724359853375113257533978570443685138259398969894101570076)
```


## Continued fractions of quadratic numbers

Since it comes essentially for free, 
another use of this package is to find the
convergents for the continued fraction of quadratic numbers,
ie., numbers of the form $(\sqrt{D} + P) / Q$,
where $D, P, Q$ are integers.

You can do this by calling
`continued_fraction(D, P=0, Q=1)`.
This returns an iterator whose eltypes are 3-Tuples of the form `(ai::BigInt, Pi::BigInt, Qi::BigInt)`,
where `ai` is the `i`'th coefficent of the continued fraction of $(\sqrt{D} + P) / Q$,
and `Pi / Qi` is its `i`'th convergent.

This can be useful if you want to approximate a floating point quadratic number
with a rational number so you can avoid floating-point arithmetic.

### Examples

```julia
# Find the first twelve convergents of (sqrt(14) - 9) / -13.
julia> collect(Iterators.take(continued_fraction(14, -9, -13), 12))
12-element Vector{Tuple{Int64, BigInt, BigInt}}:
 (0, 0, 1)
 (2, 1, 2)
 (2, 2, 5)
 (8, 17, 42)
 (1, 19, 47)
 (1, 36, 89)
 (18, 667, 1649)
 (1, 703, 1738)
 (12, 9103, 22505)
 (1, 9806, 24243)
 (18, 185611, 458879)
 (1, 195417, 483122)

julia> 667//1649 == 0 + (1//(2 + 1//(2 + 1//(8 + 1//(1 + 1//(1 + 1//18))))))
true

julia> (sqrt(14) - 9) / -13 ≈ 195417 / 483122
true
```
