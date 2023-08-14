#%% numerical integral
function trapezoidal_small(fc, x, j)
    delta_x = x[j+1] - x[j]
    mean_value = (fc(x[j+1]) + fc(x[j])) / 2
    return mean_value * delta_x
end
function trapezoidal(fc, x)
    # one dimensional trapezoidal rule
    # fc: function
    # x : interval (Vector)
    _ans = 0.0
    M = length(x)# M = N + 1
    for j in 1:M-1
        _ans += trapezoidal_small(fc, x, j)
    end
    return _ans
end
function trapezoidal_uniform(fc, (a,b), N)
    # one dimensional trapezoidal rule
    # on a uniform grid
    # fc: function
    # (a,b) : interval (tuple)
    M = N + 1
    delta_x = (b - a) / N
    x = range(a, b; length=M)
    _ans = (fc(a) + fc(b)) / 2
    for j in 2:M-1# not N-1
        _ans += fc(x[j])
    end
    _ans *= delta_x
    return _ans
end
#%%
function simpson_small(fc, x, J)
    return fc(x[J-1]) + 4 * fc(x[J]) + fc(x[J+1])
end
function simpsons_rule(fc, (a,b), n)
    # one dimensional Simpson's rule
    # on a uniform grid
    # fc: function
    # (a,b) : interval (tuple)
    N = 2n
    M = N + 1
    delta_x = (b - a) / N
    x = range(a, b; length=M)
    _ans = 0.0
    for j in 1:n
        # print(2j, " ",2j+1, " ")
        _ans += simpson_small(fc, x, 2j)
    end
    _ans *= delta_x / 3
    return _ans
end
#%%
using Combinatorics
using OffsetArrays
using LinearAlgebra

function nevilles_algorithm(T, h, j, k)
    numer = T[j,k-1] - T[j-1,k-1]
    denom = (h[j-k] / h[j])^2 - 1
    term = numer / denom
    # T[j,k] = 
    return T[j,k-1] + term
end
function lower_triangle(L)
    # without 1st column
    [(j, k ) for (j, k ) in Base.product(L,L) if j >= k >= 1]
end

function rombergs_method(fc, (a,b), N0, m)
    # N: interval width
    # m: size of a matrix - 1
    indices = 0:m
    m_plus_1 = m + 1
    T = OffsetArray(LowerTriangular(zeros(m_plus_1, m_plus_1)), indices, indices)
    h0 = (b - a) / N0
    h = [h0 / 2^j for j in indices]; h = OffsetArray(h, indices);
    N = [Int(N0 * 2^j) for j in indices]; N = OffsetArray(N, indices);

    for j in indices
        # println(N[j])
        T[j,0] = trapezoidal_uniform(fc, (a,b), Int(N[j]))
    end
    for (j,k) in lower_triangle(indices)
        T[j,k] = nevilles_algorithm(T, h, j, k)
    end
    println(T[m,m])
    return T
end

#%% distributions
using SpecialFunctions
using Distributions

function gaussian(x; mu=0, sigma=1)
    a = 1 / (sqrt(2pi) * sigma)
    return a * exp(- (x - mu)^2 / (2 * sigma^2 ))
end
function gaussian(x; kwargs=(mu=0, sigma=1))
    mu, sigma = kwargs[:mu], kwargs[:sigma]
    a = 1 / (sqrt(2pi) * sigma)
    return a * exp(- (x - mu)^2 / (2 * sigma^2 ))
end

function get_random_gaussian(size; mu=0, sigma=1)
    d = Normal(mu, sigma)
    return rand(d, size)
end
function get_random_gaussian(size; kwargs::NamedTuple)
    mu, sigma = kwargs[:mu], kwargs[:sigma]
    d = Normal(mu, sigma)
    return rand(d, size)
end
function prod_gaussian(x::Vector{Float64}; kwargs::NamedTuple)
    prod(gaussian.(x; kwargs=kwargs))
end
function uniform_dist(x::Vector{Float64}; kwargs::NamedTuple)
    n = length(x)
    1 / 2^n
end
function uniform_rand(n; kwargs = nothing)
    c = 1
    return c * (2 * rand(n) .- 1)
end

#%%
function num_integral_montecarlo(n, N, fc, prob, rand_vars; kwargs_fc, kwargs_prob, kwargs_rand)
    #=
    numerical 
    x is a vector of random variables
    a function is defined on [0,c]^n
    The distribution of x is rand_vars
    The theoretical formula is prob
    =#
    _ans = 0.0
    for _ in 1:N
        x = rand_vars(n; kwargs = kwargs_rand)# n is the number of random variables
        # @show x
        _ans += fc(x; kwargs = kwargs_fc) / prob(x; kwargs = kwargs_prob)
    end
    _ans /= N
    return _ans
end
#%%
using Statistics
function iterate_montecarlo(num_iteration)
    _ans = [num_integral_montecarlo(n, N, fc, prob, rand_vars; kwargs_fc = nothing, kwargs_prob=kwargs_normal, kwargs_rand = kwargs_normal) for _ in 1:num_iteration]
    _mean = mean(_ans)
    all_variance = mean(_ans.^2) - _mean^2
    # the square root must be operated on the standard deviation
    _stdev = sqrt(all_variance / (N - 1))

    return _mean, _stdev, (_mean - 2 * _stdev, _mean + 2 * _stdev)
end
#%%
volume_circle(n, R=1) = pi^(n / 2) * R^n / gamma(n / 2 + 1)# from 03test.jl
#%%
volume_circle(4 - 0.1im)
volume_circle(4)
volume_circle(4 + 0.1im)