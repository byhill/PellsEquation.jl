module PellsEquation

using Base.Iterators: countfrom
using Primes: factor

using("sqrtmod.jl")

export pells_eqn, pells_variant

"""
    issquare(n::Integer)
Returns `true` if and only if `n` is a perfect square
"""
@inline issquare(n::Integer) = isqrt(n)^2 == n


"""
    squarefactors(n::Integer)

Returns a list of positive integers f such that f^2 divides n.
"""
function squarefactors(n::Integer)
    return sort(foldl((a, (p, e)) -> vcat((a * [p^(i) for i in 0:e÷2]')...), factor(n), init=[one(typeof(n))]))
end


### PQa algorithm

struct PQA{T}
    D::T
    P::T
    Q::T
end

PQa(D::Integer, P::Integer, Q::Integer) = PQa(BigInt, D, P, Q)

function PQa(T::Type, D::Integer, P::Integer, Q::Integer)
    Q == 0 && throw(DomainError(Q, "Q must be nonzero"))
    D > 0 || throw(DomainError(D, "D must be postive"))
    issquare(D) && throw(DomainError(D, "D must not be square"))
    mod(powermod(P, 2, Q) - D, Q) == 0 || throw(DomainError(P, "Require P^2 ≡ D (mod Q)"))

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


# --------------------------------------------------------------------
# Pell's Equation:  x^2 - Dy^2 = 1


"""
    pellsfundamentalsoln1(D::Integer)

Finds the minimal positive solution (x₀, y₀) to x^2 - Dy^2 = ±1.
Returns (l, x₀, y₀) where l is the number of iterations of PQa
it took to find x₀, y₀.
"""
function pells_fundamental_soln_1(D::Integer)
    a₀ = isqrt(D)
    (G, B) = (0, 0)
    for pqa in PQa(D, 0, 1)
        if pqa.a == 2a₀
            (x₀, y₀) = (G, B)
            return (pqa.i, x₀, y₀)
        end
        (G, B) = (pqa.G, pqa.B)
    end
end


struct PellsEqn_1{T}
    D::T
    x₀::T
    y₀::T
end

"""
    pellsequation([T::Type=BigInt,] D::Integer)
    pellsequation([T::Type=BigInt,] D::Integer, x₀::Integer, y₀::Integer)


Iterate through all positive integer solutions ``(x, y)`` of the Diophantine equation
``x^2 - D⋅y^2 = 1``.

If a solution ``(x₀, y₀)`` is already known, supplying x₀ and y₀ as arguments
will allow the iterator to skip the costly process of finding
the fundamental solution to the Diophantine equation,
and solutions will be generated upward starting from ``(x₀, y₀)``.
"""
pells_eqn_1(D::Integer) = pells_eqn_1(BigInt, D)

function pells_eqn_1(::Type{T}, D::Integer) where {T<:Integer}
    D > 0 || throw(DomainError(D, "D must be a positive integer"))
    issquare(D) && throw(DomainError(D, "D must not be a perfect square"))

    (l, x₀, y₀) = pells_fundamental_soln_1(D)
    if isodd(l)
        (x₀, y₀) = (x₀ * x₀ + D * y₀ * y₀, 2x₀ * y₀)
    end
    return PellsEqn_1{T}(D, x₀, y₀)
end

pells_eqn_1(D::Integer, x₀::Integer, y₀::Integer) = pells_eqn_1(BigInt, D, x₀, y₀)

function pells_eqn_1(::Type{T}, D::Integer, x₀::Integer, y₀::Integer) where {T<:Integer}
    D > 0 || throw(DomainError(D, "D must be a positive integer"))
    issquare(D) && throw(DomainError(D, "D must not be a perfect square"))
    x₀^2 - D * y₀^2 == 1 || throw(ArgumentError("x₀^2 - D * y₀^2 ≠ 1"))

    return PellsEqn_1{T}(D, x₀, y₀)
end


function Base.iterate(it::PellsEqn_1, state=(it.x₀, it.y₀))
    (xₙ, yₙ) = state
    (D, x₀, y₀) = (it.D, it.x₀, it.y₀)
    return (state, (x₀ * xₙ + y₀ * yₙ * D, x₀ * yₙ + y₀ * xₙ))
end

Base.IteratorSize(::Type{PellsEqn_1{T}}) where {T} = Base.IsInfinite()
eltype(::Type{PellsEqn_1{T}}) where {T} = Tuple{T,T}


# --------------------------------------------------------------------
# Negative Pell's Equation:  x^2 - Dy^2 = -1

struct PellsEqn_Neg1{T}
    D::T
    x₀::T
    y₀::T
    u::T
    v::T
end

"""
    negative_pellsequation([T::Type=BigInt,] D::Integer)
    negative_pellsequation([T::Type=BigInt,] D::Integer, x₀::Integer, y₀::Integer)

Iterate through all positive integer solutions ``(x, y)`` of the Diophantine equation
``x^2 - D⋅y^2 = -1``.

If a solution ``(x₀, y₀)`` is already known, supplying x₀ and y₀ as arguments
will allow the iterator to skip the costly process of finding
the fundamental solution to the Diophantine equation,
and solutions will be generated upward starting from ``(x₀, y₀)``.
"""
pells_eqn_neg1(D::Integer) = pells_eqn_neg1(BigInt, D)

function pells_eqn_neg1(::Type{T}, D::Integer) where {T<:Integer}
    D > 0 || throw(DomainError(D, "D must be a positive integer"))
    issquare(D) && throw(DomainError(D, "D must not be a perfect square"))

    (l, x₀, y₀) = pells_fundamental_soln_1(D)
    if isodd(l)
        (u, v) = (x₀^2 + D * y₀^2, 2x₀ * y₀)
        return PellsEqn_Neg1{T}(D, x₀, y₀, u, v)
    end
    return PellsEqn_Neg1{T}(D, 0, 0, 0, 0)  # Has no solutions, return empty iterator.
end

pells_eqn_neg1(D::Integer, x₀::Integer, y₀::Integer) = pells_eqn_neg1(BigInt, D, x₀, y₀)

function pells_eqn_neg1(::Type{T}, D::Integer, x₀::Integer, y₀::Integer) where {T<:Integer}
    D > 0 || throw(DomainError(D, "D must be a positive integer"))
    x₀^2 - D * y₀^2 == -1 || throw(ArgumentError("$x₀^2 - $D * $y₀^2 ≠ -1"))
    return PellsEqn_Neg1{T}(D, x₀, y₀, x₀^2 + D * y₀^2, 2x₀ * y₀)
end


function Base.iterate(it::PellsEqn_Neg1)
    (x₀, y₀) = (it.x₀, it.y₀)
    (x₀, y₀) == (0, 0) && return nothing
    (D, u, v) = (it.D, it.u, it.v)

    return ((x₀, y₀), (x₀ * u + D * y₀ * v, x₀ * v + y₀ * u))
end

function Base.iterate(it::PellsEqn_Neg1, state)
    (xₙ, yₙ) = state
    (D, u, v) = (it.D, it.u, it.v)

    return ((xₙ, yₙ), (xₙ * u + D * yₙ * v, xₙ * v + yₙ * u))
end

Base.IteratorSize(::Type{PellsEqn_Neg1{T}}) where {T} = Base.SizeUnknown()
eltype(::Type{PellsEqn_Neg1{T}}) where {T} = Tuple{T,T}



# --------------------------------------------------------------------
# Negative Pell's Equation:  x^2 - Dy^2 = -1

struct PellsEqnGeneral{T}
    D::T
    u::T
    v::T
    queue::Vector{Tuple{T,T}}
end


function fundamentalsoln(::Type{T}, D::Integer, P::Integer, Q::Integer) where {T<:Integer}
    (P_reduced, Q_reduced) = (zero(T), zero(T))
    (G, B) = (abs(Q), zero(T))

    for pqa in PQa(T, D, P, Q)
        abs(pqa.Q) == 1 && return (G, B)
        (pqa.P, pqa.Q) == (P_reduced, Q_reduced) && return nothing

        (G, B) = (pqa.G, pqa.B)
        if (P_reduced, Q_reduced) == (0, 0)
            ζ = (pqa.P + sqrt(D)) / pqa.Q
            ζ̄ = (pqa.P - sqrt(D)) / pqa.Q
            if ζ > 1 && (-1 < ζ̄ < 0)
                (P_reduced, Q_reduced) = (pqa.P, pqa.Q)
            end
        end
    end
end


function fundamental_solns(::Type{T}, D::Integer, N::Integer) where {T<:Integer}
    (t, u) = isempty(pells_eqn_neg1(T, D)) ? (zero(T), zero(T)) : first(pells_eqn_neg1(T, D))
    solns = Tuple{T,T}[]

    for f in squarefactors(N)
        m = N ÷ f^2
        for z in modularsquareroots(D, abs(m))
            soln = fundamentalsoln(T, D, z, abs(m))
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

pellseqn_general(D::Integer, N::Integer) = pellseqn_general(BigInt, D, N)

function pellseqn_general(::Type{T}, D::Integer, N::Integer) where {T<:Integer}
    (t, u) = first(pells_eqn_1(D))

    solns = fundamental_solns(T, D, N)

    queue = Tuple{T,T}[]
    for (x, y) in solns
        x < 0 && ((x, y) = (-x, -y))
        y < 0 && ((x, y) = (x * t + D * y * u, x * u + y * t))
        x < 0 && ((x, y) = (-x, -y))
        push!(queue, (x, y))
    end
    sort!(queue)

    return PellsEqnGeneral{BigInt}(D, t, u, queue)
end


function Base.iterate(it::PellsEqnGeneral)
    state = it.queue
    isempty(state) && return nothing

    (D, u, v) = (it.D, it.u, it.v)
    (x, y) = popfirst!(state)
    (xₙ, yₙ) = (x * u + D * y * v, x * v + y * u)
    push!(state, (xₙ, yₙ))
    return ((x, y), state)
end


function Base.iterate(it::PellsEqnGeneral, state)
    (D, u, v) = (it.D, it.u, it.v)
    (x, y) = popfirst!(state)
    (xₙ, yₙ) = (x * u + D * y * v, x * v + y * u)
    push!(state, (xₙ, yₙ))
    return ((x, y), state)
end

Base.IteratorSize(::Type{PellsEqnGeneral{T}}) where {T} = Base.SizeUnknown()
eltype(::Type{PellsEqnGeneral{T}}) where {T} = Tuple{T,T}


# --------------------------------------------------------------------
# Public Interface


"""
    pells_eqn([T,], D::Integer, N::Integer=1)

Returns an iterator that iterates over all nonnegative integer solutions `(x, y)`
of the equation

`x^2 - D⋅y^2 = N`.
"""
pells_eqn(D::Integer, N::Integer=BigInt(1)) = pells_eqn(BigInt, D, N)

function pells_eqn(::Type{T}, D::Integer, N::Integer=T(1)) where {T<:Integer}
    D > 0 || throw(DomainError(D, "D must be a positive integer"))
    N == 0 && throw(DomainError(N, "N must be a nonzero integer"))
    issquare(D) && throw(DomainError(D, "D must not be a perfect square"))

    N == 1 && return pells_eqn_1(T, D)
    N == -1 && return pells_eqn_neg1(T, D)
    return pellseqn_general(T, D, N)
end


"""
    pells_variant([T,] a::Integer, b::Integer, c::Integer)

Returns an iterator that iterates over all nonegative integer solutions `(x, y)`
of the equation

`a⋅x^2 - b⋅y^2 = c'.
"""
pells_variant(a::Integer, b::Integer, c::Integer) = pells_variant(BigInt, a, b, c)

function pells_variant(::Type{T}, a::Integer, b::Integer, c::Integer) where {T<:Integer}
    c == 0 && throw(DomainError(c, "c must be nonzero"))
    a > 0 || throw(DomainError(a, "a must be a positive integer"))
    b > 0 || throw(DomainError(b, "b must be a positive integer"))
    D = a * b
    N = a * c
    issquare(D) && throw(DomainError(a * b, "a⋅b must not be a perfect square"))

    (u, v) = first(pells_eqn_1(T, D))
    any(soln -> mod(soln[1], a) == 0, fundamental_solns(T, D, N)) || return ()
    return Iterators.map(
        soln -> ((soln[1] ÷ a) * u + b * soln[2] * v, soln[2] * u + soln[1] * v),
        Iterators.filter(soln -> mod(soln[1], a) == 0, pells_eqn(T, D, N))
    )
end

end  # module PellsEquation
