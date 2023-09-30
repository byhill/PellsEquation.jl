module PellsEquation

using ModularSquareRoots: sqrtmod
using Primes: factor

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
    for (p, e) in factor(n)
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

    a₀ = floor(T, (it.P + sqrt(it.D)) / it.Q)
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
    aᵢ = floor(T, (Pᵢ + sqrt(D)) / Qᵢ)
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

    queue = Tuple{T,T}[(x₀, y₀)]
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
            ξ = (pqa.P + sqrt(D)) / pqa.Q
            ξ̄ = (pqa.P - sqrt(D)) / pqa.Q
            if ξ > 1 && (-1 < ξ̄ < 0)
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
    (t, u) = first(pellseqn1(D))

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

struct PellsEqn{T}
    D::T
    u::T
    v::T
    queue::Vector{Tuple{T,T}}
end

function Base.iterate(it::PellsEqn)
    state = it.queue
    isempty(state) && return nothing

    (D, u, v) = (it.D, it.u, it.v)
    (x, y) = popfirst!(state)
    (xₙ, yₙ) = (x * u + D * y * v, x * v + y * u)
    push!(state, (xₙ, yₙ))
    return ((x, y), state)
end

function Base.iterate(it::PellsEqn, state)
    (D, u, v) = (it.D, it.u, it.v)
    (x, y) = popfirst!(state)
    (xₙ, yₙ) = (x * u + D * y * v, x * v + y * u)
    push!(state, (xₙ, yₙ))
    return ((x, y), state)
end

Base.IteratorSize(::Type{PellsEqn{T}}) where {T} = Base.SizeUnknown()
eltype(::Type{PellsEqn{T}}) where {T} = Tuple{T,T}


# --------------------------------------------------------------------
# Public Interface


"""
    pells_eqn([T,], D::Integer, N::Integer=1)

Returns an iterator that iterates over all nonnegative integer solutions `(x, y)`
of the equation

`x^2 - D⋅y^2 = N`.
"""
pellseqn(D::Integer, N::Integer=1) = pellseqn(BigInt, D, N)
pellseqn(::Type{T}, D::Integer, N::Integer=1) where {T<:Integer} = pellseqn(T, convert(T, D), convert(T, N))

function pellseqn(::Type{T}, D::T, N::T=one(T)) where {T<:Integer}
    isless(0, D) || throw(DomainError(D, "D must be a positive integer"))
    iszero(N) && throw(DomainError(N, "N must be a nonzero integer"))
    issquare(D) && throw(DomainError(D, "D must not be a perfect square"))

    isone(N) && return pellseqn1(D)
    isone(-N) && return pellseqn_neg1(D)
    return _pellseqn(D, N)
end


end
