#%%
using IterTools
using Combinatorics
using ProgressMeter
#%%

function identity_relation(n::Int64; start=0)#the identity element
    return [(j, j) for j in start:(n+start-1)]
end
function identity_relation(vector::Vector)
    return [(vj, vj) for vj in vector]
end
function lattice_point(n::Int64; start=0::Int64)
    X =start:(n+start-1)
    return [p for p in product(X,X)]
end

function direct_prod(X,Y)
    return [(x,y) for x in X, y in Y]
end

function composition_of_relation(R::Vector, S::Vector)
    first_type = Base.tuple_type_head(eltype(R))
    second_type = Base.tuple_type_head(Base.tuple_type_tail(eltype(S)))
    
    _ans = Vector{Tuple{first_type, second_type}}()
    for x in R, y in S
        if (x[2] == y[1])# (a, b) ∈ R and (b, c) ∈ S
            push!(_ans, (x[1], y[2]))
        end
    end
    return _ans
end
#%%

function identity_relation_set(n::Int64; start=0)#the identity element
    return Set(identity_relation(n; start))
end
function lattice_point_set(n::Int64; start=0::Int64)#faster
    return Set(lattice_point(n; start))
end

function direct_prod_set(X,Y)
    return Set(direct_prod(X,Y))
end
function composition(R::Vector{Tuple{Int64, Int64}},S::Vector{Tuple{Int64, Int64}})
    _ans = Set{Tuple{Int64, Int64}}()
    for x in R, y in S
        if (x[2] == y[1])
            push!(_ans,(x[1], y[2]))
        end
    end
    return _ans
end
#%%
function FindInvertible(n::Int64)
    if n ∉ 1:3
        throw(DomainError(n, "n must be in {1, 2, 3}"))
    end
    X2 = lattice_point_set(n; start=0)
    PX2 = powerset(collect(X2))
    itr = 0# the number of the invertibles
    for R in PX2, S in PX2 
        if composition(R,S) == identity_relation_set(n)
            @show R, S
            itr += 1
        end
    end
    return itr
end

function invertible_test()
    @time for n in 1:3
        @show n
        FindInvertible(n)
    end
    return 0
end
#%%
function isreflexive(R; identity)
    issubset(identity, Set(R))# true if reflexive
end
function istransitive(R)
    issubset(composition(R,R),R)# true if transitive
end
#%%
function transpose_rel(R)
    Set([(x[2],x[1]) for x in R])
end
function issymmetric(R)
    issetequal(transpose_rel(R),R)# true if symmetric
end
#%%
function FindEquivalence(n)
    X2 = lattice_point_set(n; start=0)
    PX2 = powerset(collect(X2))
    itr = 0
    I = identity_relation_set(n; start=0)
    for R in PX2
        if (isreflexive(R; identity=I) && istransitive(R) && issymmetric(R))
            @show R
            itr += 1
        end
    end
    return itr
end
# @time FindEquivalence(4)

#%%
# 2

# (R, S) = ([(0, 0), (1, 1)], [(0, 0), (1, 1)])
# (R, S) = ([(0, 1), (1, 0)], [(0, 1), (1, 0)])

# 3

# (R, S) = ([(0, 0), (1, 2), (2, 1)], [(0, 0), (1, 2), (2, 1)])
# (R, S) = ([(0, 0), (1, 1), (2, 2)], [(0, 0), (1, 1), (2, 2)])
# (R, S) = ([(1, 2), (2, 0), (0, 1)], [(0, 2), (2, 1), (1, 0)])
# (R, S) = ([(2, 0), (1, 1), (0, 2)], [(2, 0), (1, 1), (0, 2)])
# (R, S) = ([(0, 1), (2, 2), (1, 0)], [(0, 1), (2, 2), (1, 0)])
# (R, S) = ([(0, 2), (2, 1), (1, 0)], [(1, 2), (2, 0), (0, 1)])
#   0.250575 seconds (2.89 M allocations: 342.603 MiB, 20.54% gc time)