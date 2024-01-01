# multi-valued logic (three)
# %%
function land(P,Q)
    min(P,Q)
end
function lor(P,Q)
    max(P,Q)
end

function simpleNegation(P)
    1 - P
end
function totalNegThree(x)
    (x == 0.0) * 1 + (x == 0.5) * 0 + (x == 1.0) * 0
end
function partialNegThree(x)
    (x == 0.0) * 1 + (x == 0.5) * 1 + (x == 1.0) * 0
end
function partialNegThreeByTripleNeg(P)
    simpleNegation(totalNegation(simpleNegation(P)))
end

totalNegation(P) = totalNegThree(P)
partialNegation(P) = partialNegThree(P)
# function totalNegation(P)
#     if P >= .5
#         return 0
#     else
#         return 1
#     end
# end
# function partialNegation(P)
#     if P <= .5
#         return 1
#     else
#         return 0
#     end
# end

NAND(P,Q) = simpleNegation(land(P,Q))
NOR(P,Q) = simpleNegation(lor(P,Q))
NAND.(P,P)

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
    min(1, 1-P+Q)
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

simpleEquiv(P,Q) = equivalence.(P,Q;imply=simpleImplication)

EquivLukasiewicz(P,Q) = equivalence.(P,Q;imply=ImplyLukasiewicz)

lnotLukasiewicz(P) = simpleImplication.(P,land.(P, totalNegation.(P)))#20221220

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
U = Iterators.repeated(0.5,3)

n2 = 2
P = 0:(1/n2):1
Q = 0:(1/n2):1

f00u(P) = land(0.5, totalAffirmation(P))
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
    #Test code
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

#%%two-variable functions
# basis
f000000001(P,Q) = land(totalAffirmation(P),totalAffirmation(Q))
f000000010(P,Q) = land(totalAffirmation(P),looseContradiction(Q))
f000010000(P,Q) = land(looseContradiction(P),looseContradiction(Q))#for Lukasiewicz imply
f00000000u(P,Q) = land(f00u(P),f00u(Q))
f0000000u0(P,Q) = land(f00u(P),simpleContradiction(Q))
f0000u0000(P,Q) = land(simpleContradiction(P),simpleContradiction(Q))

f010111010(P,Q) = lor(looseContradiction(P),looseContradiction(Q))
f111u110u1(P,Q) = lor(f000010000(P,Q),simpleImplication(P,Q))

f000001000(P,Q) = f000000010(Q,P)
f000000100(P,Q) = f000000001(P,simpleNegation(Q))
f000100000(P,Q) = f000000010(simpleNegation(Q),P)
f100000000(P,Q) = f000000001(simpleNegation(P),simpleNegation(Q))
# f0(P,Q) =  land(f00u(P),simpleContradiction(Q))
matrixof(f0,P,Q)
