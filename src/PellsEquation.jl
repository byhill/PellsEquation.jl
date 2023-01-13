module PellsEquation

export PQa

export pellsequation, negative_pellsequation

"""
    issquare(n::Integer)
Returns `true` if and only if `n` is a perfect square
"""
@inline issquare(n::Integer) = isqrt(n)^2 == n


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

Base.IteratorSize(::Type{PellsEqn{T}}) where {T} = Base.IsInfinite()


# --------------------------------------------------------------------
# Pell's Equation:  x^2 - Dy^2 = 1

struct PellsEqn{T}
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
pellsequation(D::Integer) = pellsequation(BigInt, D)

function pellsequation(::Type{T}, D::Integer) where {T<:Integer}
    D > 0 || throw(DomainError(D, "D must be a positive integer"))
    issquare(D) && throw(DomainError(D, "D must not be a perfect square"))

    pqaI = Iterators.Stateful(PQa(D, 0, 1))
    a₀ = isqrt(D)
    l = 0
    for pqa in pqaI
        pqa.a == 2a₀ || continue
        l = pqa.i
        break
    end

    for _ in 1:l-2
        popfirst!(pqaI)
    end
    pqa = popfirst!(pqaI)

    x₀ = pqa.G
    y₀ = pqa.B
    return PellsEqn{T}(D, x₀, y₀)
end

pellsequation(D::Integer, x₀::Integer, y₀::Integer) = pellsequation(BigInt, D, x₀, y₀)

function pellsequation(::Type{T}, D::Integer, x₀::Integer, y₀::Integer) where {T<:Integer}
    D > 0 || throw(DomainError(D, "D must be a positive integer"))
    issquare(D) && throw(DomainError(D, "D must not be a perfect square"))
    x₀^2 - D * y₀^2 == 1 || throw(ArgumentError("x₀^2 - D * y₀^2 ≠ 1"))
    return PellsEqn{T}(D, x₀, y₀)
end

function Base.iterate(it::PellsEqn, state=(it.x₀, it.y₀))
    (xₙ, yₙ) = state
    (D, x₀, y₀) = (it.D, it.x₀, it.y₀)
    return (state, (x₀ * xₙ + y₀ * yₙ * D, x₀ * yₙ + y₀ * xₙ))
end

Base.IteratorSize(::Type{PellsEqn{T}}) where {T} = Base.SizeUnknown()
eltype(::Type{PellsEqn{T}}) where {T} = Tuple{T,T}


# --------------------------------------------------------------------
# Negative Pell's Equation:  x^2 - Dy^2 = -1

struct NegPellsEqn{T}
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
negative_pellsequation(D::Integer) = negative_pellsequation(BigInt, D)

function negative_pellsequation(::Type{T}, D::Integer) where {T<:Integer}
    D > 0 || throw(DomainError(D, "D must be a positive integer"))
    issquare(D) && throw(DomainError(D, "D must not be a perfect square"))

    a₀ = isqrt(D)
    (G, B) = (0, 0)
    for pqa in PQa(D, 0, 1)
        if pqa.a == 2a₀
            (x₀, y₀) = (G, B)
            (u, v) = (x₀^2 + D * y₀^2, 2x₀ * y₀)
            return isodd(pqa.i) ? NegPellsEqn{T}(D, x₀, y₀, u, v) : NegPellsEqn{T}(D, 0, 0, 0, 0)
        end
        (G, B) = (pqa.G, pqa.B)
    end

    return nothing
end

negative_pellsequation(D::Integer, x₀::Integer, y₀::Integer) = negative_pellsequation(BigInt, D, x₀, y₀)

function negative_pellsequation(::Type{T}, D::Integer, x₀::Integer, y₀::Integer) where {T<:Integer}
    D > 0 || throw(DomainError(D, "D must be a positive integer"))
    x₀^2 - D * y₀^2 == -1 || throw(ArgumentError("$x₀^2 - $D * $y₀^2 ≠ -1"))
    return NegPellsEqn{T}(D, x₀, y₀, x₀^2 + D * y₀^2, 2x₀ * y₀)
end


function Base.iterate(it::NegPellsEqn)
    (x₀, y₀) = (it.x₀, it.y₀)
    (x₀, y₀) == (0, 0) && return nothing
    (D, u, v) = (it.D, it.u, it.v)

    return ((x₀, y₀), (x₀ * u + D * y₀ * v, x₀ * v + y₀ * u))
end

function Base.iterate(it::NegPellsEqn, state)
    (xₙ, yₙ) = state
    (D, u, v) = (it.D, it.u, it.v)

    return ((xₙ, yₙ), (xₙ * u + D * yₙ * v, xₙ * v + yₙ * u))
end

Base.IteratorSize(::Type{NegPellsEqn{T}}) where {T} = Base.SizeUnknown()
eltype(::Type{NegPellsEqn{T}}) where {T} = Tuple{T,T}

end  # module PellsEquation
