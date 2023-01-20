function modularsquareroots(n::Integer, m::Integer)
    solns = []
    n = mod(n, m)
    for x = 0:m÷2
        powermod(x, 2, m) == n && (push!(solns, x), push!(solns, -x))
    end
    return solns
end