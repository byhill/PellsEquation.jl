module PellsEquation

using ModularSquareRoots: sqrtmod
using Primes: eachfactor, divisors
using .Iterators

export pellseqn


"""
    issquare(n::Integer)
Returns `true` if and only if `n` is a perfect square
"""
@inline issquare(n::Integer) = isqrt(n)^2 == n


"""
    squarefactors(n::Integer)

Returns a list of positive integers d such that d^2 divides n.
"""
function squarefactors(n::T) where {T<:Integer}
    factors = [one(T)]
    for (p, e) in eachfactor(n)
        x = [p^i for i in 1:e÷2]
        L = length(factors)
        for d in factors[1:L], n in x
            push!(factors, d * n)
        end
    end
    sort!(factors)

    return factors
end


# --------------------------------------------------------------------
# PQa algorithm

struct PQA{T}
    D::T
    P::T
    Q::T
end


PQa(D::T, P::Integer, Q::Integer) where {T<:Integer} = PQa(D, convert(T, P), convert(T, Q))
function PQa(D::T, P::T, Q::T) where {T<:Integer}
    iszero(Q) && throw(DomainError(Q, "Q must be nonzero in the PQa algorithm"))
    isless(0, D) || throw(DomainError(D, "D must be postive in the PQa algorithm"))
    issquare(D) && throw(DomainError(D, "D must not be square in the PQa algorithm"))
    mod(powermod(P, 2, Q) - D, Q) == 0 || throw(DomainError(P, "Require P² ≡ D (mod Q) in the PQa algorithm"))

    return PQA{T}(D, P, Q)
end


function Base.iterate(it::PQA{T}) where {T<:Integer}
    (A₋₂, A₋₁) = (zero(T), one(T))
    (B₋₂, B₋₁) = (one(T), zero(T))
    (G₋₂, G₋₁) = (-it.P, it.Q)

    a₀ = (it.P + isqrt(it.D)) ÷ it.Q
    A₀ = a₀ * A₋₁ + A₋₂
    B₀ = a₀ * B₋₁ + B₋₂
    G₀ = a₀ * G₋₁ + G₋₂

    ret = (i=0, P=it.P, Q=it.Q, a=a₀, A=A₀, B=B₀, G=G₀)
    state = (0, it.P, it.Q, a₀, (A₀, A₋₁), (B₀, B₋₁), (G₀, G₋₁))
    return (ret, state)
end


function Base.iterate(it::PQA{T}, state) where {T<:Integer}
    (i, Pᵢ₋₁, Qᵢ₋₁, aᵢ₋₁, (Aᵢ₋₁, Aᵢ₋₂), (Bᵢ₋₁, Bᵢ₋₂), (Gᵢ₋₁, Gᵢ₋₂)) = state
    D = it.D

    i += 1
    Pᵢ = aᵢ₋₁ * Qᵢ₋₁ - Pᵢ₋₁
    Qᵢ = (D - Pᵢ^2) ÷ Qᵢ₋₁
    aᵢ = (Pᵢ + isqrt(D)) ÷ Qᵢ
    Aᵢ = aᵢ * Aᵢ₋₁ + Aᵢ₋₂
    Bᵢ = aᵢ * Bᵢ₋₁ + Bᵢ₋₂
    Gᵢ = aᵢ * Gᵢ₋₁ + Gᵢ₋₂

    ret = (i=i, P=Pᵢ, Q=Qᵢ, a=aᵢ, A=Aᵢ, B=Bᵢ, G=Gᵢ)
    state = (i, Pᵢ, Qᵢ, aᵢ, (Aᵢ, Aᵢ₋₁), (Bᵢ, Bᵢ₋₁), (Gᵢ, Gᵢ₋₁))
    return (ret, state)
end

Base.IteratorSize(::Type{PQA{T}}) where {T} = Base.IsInfinite()
eltype(::Type{PQA{T}}) where {T} = NamedTuple{(:i, :P, :Q, :a, :A, :B, :G),Tuple{Int64,Vararg{T,6}}}


# --------------------------------------------------------------------
# Pell's Equation:  x^2 - Dy^2 = 1


"""
    pellsfundamentalsoln1(D::Integer)

Finds the minimal positive solution (x₀, y₀) to x^2 - Dy^2 = ±1.
Returns (l, x₀, y₀) where l is the number of iterations of PQa
it took to find x₀, y₀.
"""
function fundamental_soln(D::T) where {T<:Integer}
    a₀ = isqrt(D)
    (G, B) = (zero(T), zero(T))
    for pqa in PQa(D, 0, 1)
        if pqa.a == 2a₀
            (x₀, y₀) = (G, B)
            return (pqa.i, x₀, y₀)
        end
        (G, B) = (pqa.G, pqa.B)
    end
end


function pellseqn1(D::T) where {T<:Integer}
    (l, x₀, y₀) = fundamental_soln(D)
    if isodd(l)  # Then the fundamental solution we found is for x^2 - Dy^2 = -1
        (x₀, y₀) = (x₀ * x₀ + D * y₀ * y₀, 2x₀ * y₀)
    end

    queue = Tuple{T,T}[(1, 0)]
    return PellsEqn{T}(D, x₀, y₀, queue)
end


# --------------------------------------------------------------------
# Negative Pell's Equation:  x^2 - Dy^2 = -1


function pellseqn_neg1(D::T) where {T<:Integer}
    (l, x₀, y₀) = fundamental_soln(D)

    # Has no solutions, return empty iterator
    iseven(l) && return PellsEqn{T}(D, 0, 0, Tuple{T,T}[])

    queue = Tuple{T,T}[(x₀, y₀)]
    (u, v) = (x₀ * x₀ + D * y₀ * y₀, 2x₀ * y₀)
    return PellsEqn{T}(D, u, v, queue)
end


# --------------------------------------------------------------------
# Generalized Pell's Equation:  x^2 - Dy^2 = N


function fundamental_soln(D::T, P::T, Q::T) where {T<:Integer}
    (P_reduced, Q_reduced) = (zero(T), zero(T))
    (G, B) = (abs(Q), zero(T))

    for pqa in PQa(D, P, Q)
        abs(pqa.Q) == 1 && return (G, B)
        (pqa.P, pqa.Q) == (P_reduced, Q_reduced) && return nothing

        (G, B) = (pqa.G, pqa.B)
        if (P_reduced, Q_reduced) == (0, 0)
            ξ = (pqa.P + isqrt(D)) ÷ pqa.Q
            ξ̄ = (pqa.P - isqrt(D)) ÷ pqa.Q
            if ξ ≥ 1 && ξ̄ == 0
                (P_reduced, Q_reduced) = (pqa.P, pqa.Q)
            end
        end
    end
end


function fundamental_solns(D::T, N::T) where {T<:Integer}
    (t, u) = isempty(pellseqn_neg1(D)) ? (zero(T), zero(T)) : first(pellseqn_neg1(D))
    solns = Tuple{T,T}[]

    for f in squarefactors(N)
        m = N ÷ f^2
        for z in sqrtmod(D, abs(m))
            soln = fundamental_soln(D, z, abs(m))
            isnothing(soln) && continue
            (r, s) = soln
            if r^2 - D * s^2 == m
                push!(solns, (f * r, f * s))
            elseif (t, u) != (0, 0)
                push!(solns, (f * (r * t + D * s * u), f * (r * u + s * t)))
            end
        end
    end
    unique!(solns)

    return solns
end


function _pellseqn(D::T, N::T) where {T<:Integer}
    i, t, u = fundamental_soln(D)
    if isodd(i)  # Then the fundamental solution we found is for x^2 - Dy^2 = -1
        t, u = t * t + D * u * u, 2t * u
    end

    queue = Tuple{T,T}[]
    for (x, y) in fundamental_solns(D, N)
        isless(x, 0) && ((x, y) = (-x, -y))
        isless(y, 0) && ((x, y) = (x * t + D * y * u, x * u + y * t))
        isless(x, 0) && ((x, y) = (-x, -y))
        push!(queue, (x, y))
    end
    sort!(queue)

    return PellsEqn{T}(D, t, u, queue)
end


# --------------------------------------------------------------------
# Iterator implementation

@enum SolutionType begin  # x^2 - Dy^2 = N, where...
    pellsolution  # D nonsquare, N nonzero
    zerocoeff  # D or N are zero and the other is square
    finite  # D > 0 square, or N is zero and D is nonsquare
end

struct PellsEqn{T}
    D::T
    u::T
    v::T
    queue::Vector{Tuple{T,T}}
    solntype::SolutionType
end

PellsEqn{T}(D, u, v, queue) where {T} = PellsEqn{T}(D, u, v, queue, pellsolution)


function Base.iterate(it::PellsEqn)
    state = it.queue
    isempty(state) && return nothing

    (x, y) = popfirst!(state)
    (D, u, v) = (it.D, it.u, it.v)
    if it.solntype == pellsolution
        (xₙ, yₙ) = (x * u + D * y * v, x * v + y * u)
        push!(state, (xₙ, yₙ))
    elseif it.solntype == zerocoeff
        (xₙ, yₙ) = (x + u, y + v)
        push!(state, (xₙ, yₙ))
    end

    return ((x, y), state)
end

function Base.iterate(it::PellsEqn, state)
    isempty(state) && return nothing

    (x, y) = popfirst!(state)
    (D, u, v) = (it.D, it.u, it.v)
    if it.solntype == pellsolution
        (xₙ, yₙ) = (x * u + D * y * v, x * v + y * u)
        push!(state, (xₙ, yₙ))
    elseif it.solntype == zerocoeff
        (xₙ, yₙ) = (x + u, y + v)
        push!(state, (xₙ, yₙ))
    end

    return ((x, y), state)
end

Base.IteratorSize(::Type{PellsEqn{T}}) where {T} = Base.SizeUnknown()
Base.eltype(::Type{PellsEqn{T}}) where {T} = Tuple{T,T}


# --------------------------------------------------------------------
# Public Interface


"""
    pells_eqn(D::Integer, N::Integer=1)

Returns an iterator that yields all nonnegative integer solutions `(x, y)`
of the equation

`x^2 - D⋅y^2 = N`.

The eltype of pells_eqn is always `Tuple{BigInt,BigInt}`.

!!! note
    There are often an infinite number of nonnegative solutions to Pell's equation,
    thus the iterator can have infinite length.

```
julia> using .Iterators

julia> collect(take(pellseqn(2), 5))
5-element Vector{Tuple{BigInt, BigInt}}:
 (3, 2)
 (17, 12)
 (99, 70)
 (577, 408)
 (3363, 2378)

julia> collect(take(pellseqn(2, 1), 5))
5-element Vector{Tuple{BigInt, BigInt}}:
 (3, 2)
 (17, 12)
 (99, 70)
 (577, 408)
 (3363, 2378)

julia> collect(take(pellseqn(5, -1), 5))
5-element Vector{Tuple{BigInt, BigInt}}:
 (2, 1)
 (38, 17)
 (682, 305)
 (12238, 5473)
 (219602, 98209)

julia> collect(take(pellseqn(4, 25), 5))
2-element Vector{Tuple{BigInt, BigInt}}:
 (5, 0)
 (13, 6)

julia> collect(take(pellseqn(4, 10), 5))
Tuple{BigInt, BigInt}[]

julia> collect(take(pellseqn(921, -12), 5))
5-element Vector{Tuple{BigInt, BigInt}}:
 (123013602, 4053436)
 (620494807415610916348542, 20445999054572056617484)
 (3129847429634131057287425640460782483138, 103131979224431202543397867947637928044)
 (15787311699815721226111492312624904230758974206686324318, 520209607285982995690281237198588485553526346187087196)
 (79633022474924155281204153475574155565247711565894067613899610557324322, 2623997304693724359853375113257533978570443685138259398969894101570076)
```
"""
pellseqn(D::Integer, N::Integer=1) = pellseqn(big(D), big(N))

function pellseqn(D::BigInt, N::BigInt=big(1))
    D < 0 && throw(DomainError(D, "D must be a nonnegative integer"))

    T = BigInt
    TT = Tuple{T,T}
    if iszero(D)
        # Has solutions iff N is square
        # Solutions are (sqrt(N), 0), (sqrt(N), 1), (sqrt(N), 2), ...
        if issquare(N)
            return PellsEqn(D, zero(T), one(T), [(isqrt(N), zero(T))], zerocoeff)
        else
            return PellsEqn(D, zero(T), one(T), TT[], finite)
        end

    elseif iszero(N)
        # If D is square,
        # solutions are (0, 0), (sqrt(D), 1), (2sqrt(D), 2), ...
        # else, the only solution is (0, 0)
        if issquare(D)
            return PellsEqn(D, isqrt(D), one(D), [(zero(T), zero(T))], zerocoeff)
        else
            return PellsEqn(D, isqrt(D), one(D), [(zero(T), zero(T))], finite)
        end

    elseif issquare(D)
        D = isqrt(D)
        # N = x^2 - D^2y^2 = (x - Dy)(x + Dy) = d1 * d2
        solns = TT[]
        for d1 in divisors(N)
            N < 0 && (d1 = -d1)
            d2 = N ÷ d1

            x = d1 + d2
            iseven(x) || continue
            x ≥ 0 || continue
            x >>= 1

            y = (d2 - d1)
            y ≥ 0 || continue
            iszero(mod(y, 2D)) || continue
            y ÷= 2D

            push!(solns, (x, y))
        end

        return PellsEqn(D, zero(T), zero(T), sort(solns), finite)

    elseif isone(N)
        return pellseqn1(D)
    elseif isone(-N)
        return pellseqn_neg1(D)
    else
        return _pellseqn(D, N)
    end
end


end
