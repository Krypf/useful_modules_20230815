# multi-valued logic (three)
# %%

# Constants for the three truth values
const T_num = 1.0
const F_num = 0.0
const U_num = 0.5

# Define a custom type for three-valued logic
struct ThreeValuedLogic
    value::Symbol
    number::Float64
    # Constructor with constraint
    function ThreeValuedLogic(value::Symbol, number::Float64)
        if number ∈ (F_num, U_num, T_num)
            return new(value, number)
        else
            throw(ArgumentError("The 'number' $number is invalid"))
        end
    end
end

const _true = ThreeValuedLogic(:True, T_num)
const _false = ThreeValuedLogic(:False, F_num)
const _unknown = ThreeValuedLogic(:Unknown, U_num)

#%%

function land(P::Float64, Q::Float64)# Logical and function
    min(P, Q)
end
function lor(P::Float64, Q::Float64)# Logical or function
    max(P, Q)
end

function simpleNegation(P::Float64)
    # Tabel
    # T | F
    # U | U
    # F | T
    T_num - P
end
function totalNegThree(x::Float64)
    # Tabel
    # T | F
    # U | F
    # F | T
    (x == F_num) * T_num + (x == U_num) * F_num + (x == T_num) * F_num
end
function partialNegThree(x::Float64)
    # Tabel
    # T | F
    # U | T
    # F | T
    (x == F_num) * T_num + (x == U_num) * T_num + (x == T_num) * F_num
end
#%%
function totalNegation(P::Float64)
    condition = (P >= U_num)
    if condition
        return F_num
    else
        return T_num
    end
end
function totalNegation2(P::Float64)
    condition = (P >= U_num)
    return (condition) * F_num + !condition * T_num
end

#%%

function partialNegation(P)
    condition = (P <= U_num)
    if condition
        return T_num
    else
        return F_num
    end
end
#%%

NAND(P,Q) = simpleNegation(land(P,Q))
NOR(P,Q) = simpleNegation(lor(P,Q))

function simpleImplication(P,Q)
    lor(simpleNegation(P),Q)
end
function simpleImplication2(P,Q)
    simpleNegation(
        land((P),simpleNegation(Q))
    )
end
function Implication(P,Q;NegationQ,Negation)
    Negation(
        land((P),NegationQ(Q))
    )
end
strictImplication(P,Q) = Implication(P,Q;
NegationQ=partialNegation,
Negation=totalNegation,
)
looseImplication(P,Q) = Implication(P,Q;
NegationQ=totalNegation,
Negation=partialNegation,
)

function lorDeMorgan(P,Q)
    simpleNegation(
        land(
            simpleNegation(P),
            simpleNegation(Q)
        )
    )
end
function totalAffirmation(P)
    simpleNegation(partialNegation(P))
end
function partialAffirmation(P)
    simpleNegation(totalNegation(P))
end

function looseContradiction(P)
    land.(partialAffirmation.(P),partialNegation.(P))
end
function looseExcludedMiddle(P)
    lorDeMorgan.(totalNegation.(P),totalAffirmation.(P))
end# simpleNegation(looseContradiction())

function ImplyLukasiewicz(P,Q)
    min(T_num, T_num - P + Q)
end
ImplyInvLukasiewicz(P,Q) = ImplyLukasiewicz(Q,P)

function lorNaiveLukasiewicz(P,Q)
    simpleImplication(simpleImplication(P,Q),Q)
    # (A -> B) -> B
    # (-A | B) -> B
    # -(-A | B) | B
    # (A & -B) | B
    # (A | B) & (-B | B)
    # (A | B)
end
function lorLukasiewicz(P,Q)
    ImplyLukasiewicz(ImplyLukasiewicz(P,Q),Q)
end# 20221220

function equivalence(P,Q;imply,land=land)
    land(imply(P,Q),imply(Q,P))
end

simpleEquiv(P,Q) = equivalence(P,Q;imply=simpleImplication)

EquivLukasiewicz(P,Q) = equivalence(P,Q;imply=ImplyLukasiewicz)

# logical not in Lukasiewicz logic is equivalent to the simple negation (not).
lnotLukasiewicz(P) = simpleImplication(P,land(P, totalNegation(P)))#20221220

function modusPonens(P,Q;
    AND=land,ImplyPQ=simpleImplication,Therefore=simpleImplication)
    x = AND(ImplyPQ(P,Q),P)
    return Therefore(x,Q)
end
# modusPonens.(P,Q;ImplyPQ=ImplyLukasiewicz)
# modusPonens.(P,Q;Therefore=ImplyLukasiewicz)
# modusPonens.(P,Q;ImplyPQ=ImplyLukasiewicz,Therefore=ImplyLukasiewicz)

# %%
function matrixof(f,P,Q)# matrix of a two-variable logical function
    N = size(P)[1]
    _input = Vector{Tuple{Float64, Float64}}(undef,N*N)
    _input = reshape(_input,(N,N))# matrix of tuples of an input
    _ans = zeros(N,N)
    for j in 1:N
        for k in 1:N
            _input[j,k] = (P[j],Q[k])
            _ans[j,k] = f(P[j],Q[k])
        end
    end
    # @show _input
    return _ans
end
function tuplesof(f,P,Q)# matrix of a two-variable logical function
    N = size(P)[1]
    _input_output = Vector{Tuple{Float64, Float64, Float64}}(undef,N*N)
    _input_output = reshape(_input_output,(N,N))# matrix of tuples of an input
    for j in 1:N
        for k in 1:N
            _input_output[j,k] = (P[j],Q[k],f(P[j],Q[k]))
        end
    end
    show(stdout, "text/plain", _input_output)
    println()
    return _input_output
end
#%%
function genOneVarAll(N)
    OneVarAll = [Vector{Float64}(undef,N) for j in 1:N^N]
    _index = 1  
    for j in 1:N, k in 1:N, l in 1:N
        OneVarAll[_index] = [P[j], P[k], P[l]]
        _index += 1
    end
    return OneVarAll
end

#%%
# ALL UNKNOWN
U = Iterators.repeated(U_num,3)

n2 = 2
P = 0:(1/n2):1
Q = 0:(1/n2):1

f00u(P) = land(U_num, totalAffirmation(P))
simpleContradiction(P) = land(P, simpleNegation(P))
LogicalFunctionsThree = [
    land.(P, totalNegation.(P)),
    f00u.(P),
    totalAffirmation.(P),
    simpleContradiction.(P),
    land.(U,P),
    identity.(P),
    looseContradiction.(P),#[0,1,0]
    lorDeMorgan.(land.(U,P),looseContradiction.(P)),
    partialAffirmation.(P),
    f00u.(simpleNegation.(P)),# land.(U, totalNegation.(P)),
    lorDeMorgan.(land.(U, totalAffirmation.(P)),land.(U, totalNegation.(P))),#[0,0,0.5] or [0.5,0,0]
    lorDeMorgan.(totalAffirmation.(P),land.(U, totalNegation.(P))),#[0,0,1] or [0.5,0,0]
    land.(U, simpleNegation.(P)),
    collect(U),
    lorDeMorgan.(U,P),
    lorDeMorgan.(land.(U, totalNegation.(P)),looseContradiction.(P)),#[0.5,1,0]=[0.5,0,0] or [0,1,0]
    lorDeMorgan.(U,looseContradiction.(P)),#[0.5,1,0.5]
    lorDeMorgan.(U,partialAffirmation.(P)),

    totalNegation.(P),
    lorDeMorgan.(totalNegation.(P),land.(U, totalAffirmation.(P))),# [1,0,0.5]
    looseExcludedMiddle.(P),
    simpleNegation.(P),
    lorDeMorgan.(U,simpleNegation.(P)),
    lorDeMorgan.(P,simpleNegation.(P)),
    partialNegation.(P),
    lorDeMorgan.(U,partialNegation.(P)),
    lorDeMorgan.(P,partialNegation.(P)),
]
#%%
if abspath(PROGRAM_FILE) == @__FILE__
    # Test code
    @show genOneVarAll(3) == LogicalFunctionsThree
    # double negation
    @show partialNegation.(partialNegation.(P)) == totalAffirmation.(P)
    @show totalNegation.(partialNegation.(P)) == totalAffirmation.(P)
    @show totalNegation.(totalNegation.(P)) == partialAffirmation.(P)
    @show partialNegation.(totalNegation.(P)) == partialAffirmation.(P)

    @show totalNegation.(simpleNegation.(P)) == totalAffirmation.(P)
    @show partialNegation.(simpleNegation.(P)) == partialAffirmation.(P)

    @show simpleNegation.(totalNegation.(simpleNegation.(P))) == partialNegation.(P)
    @show partialNegation.(P) == partialNegThree.(P) == partialNegThreeByTripleNeg.(P)
end

#%% two-variable functions
# basis
# 反転を考えると右下4分の1を考えるだけで済む
# Given the inversion, we only need to consider the lower right quarter.
f000000001(P,Q) = land(totalAffirmation(P),totalAffirmation(Q))
f000000010(P,Q) = land(totalAffirmation(P),looseContradiction(Q))
f000010000(P,Q) = land(looseContradiction(P),looseContradiction(Q))# for Lukasiewicz imply

# inversion across P = Q
f000001000(P,Q) = f000000010(Q,P)
# inversion with respect to Q
f000000100(P,Q) = f000000001(P,simpleNegation(Q))
# other two examples
f000100000(P,Q) = f000000010(simpleNegation(Q),P)
f100000000(P,Q) = f000000001(simpleNegation(P),simpleNegation(Q))

# u basis
f00000000u(P,Q) = land(f00u(P),f00u(Q))
f0000000u0(P,Q) = land(f00u(P),simpleContradiction(Q))
f0000u0000(P,Q) = land(simpleContradiction(P),simpleContradiction(Q))

f010111010(P,Q) = lor(looseContradiction(P),looseContradiction(Q))
f111u110u1(P,Q) = lor(f000010000(P,Q),simpleImplication(P,Q))# ImplyLukasiewicz

println("MultiValuedLogic.jl imported")