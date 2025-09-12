export LambertProblem, LambertSolution, solve, init, solve!, select_lambert_algorithm

abstract type AbstractAstroProblem <: SciMLBase.AbstractSciMLProblem end
abstract type AbstractAstroSolution end

"""
    LambertProblem{TT, MT}

Represents a Lambert problem with gravitational parameter, position vectors, and time of flight.
Inherits from AbstractSciMLProblem to enable use with EnsembleProblem and other SciML interfaces.

When solved without specifying an algorithm, the system will automatically select the optimal
solver based on problem characteristics (transfer angle, number of revolutions, etc.).

# Fields
- `μ::MT`: Gravitational parameter (GM of attracting body)
- `r1::Vector`: Initial position vector
- `r2::Vector`: Final position vector  
- `tof::TT`: Time of flight between r1 and r2

# Examples
```julia
# Automatic algorithm selection
problem = LambertProblem(μ, r1, r2, tof)
solution = solve(problem)  # Heuristic selects best algorithm

# Multi-revolution with automatic selection
solution = solve(problem, M=2, prograde=true)  # Will select Arora or Izzo

# Manual algorithm selection
solution = solve(problem, IzzoSolver())
```
"""
struct LambertProblem{TT<:Number,RT1<:Number,RT2<:Number,MT<:Number} <: AbstractAstroProblem
    μ::MT
    r1::Vector{RT1}
    r2::Vector{RT2}
    tof::TT
end

function LambertProblem(μ::Number, r1::AbstractVector, r2::AbstractVector, tof::Number)
    return LambertProblem(μ, r1[1:3], r2[1:3], tof)
end

function LambertProblem(
    μ::Number,
    coord1::AstroCoords.AstroCoord,
    coord2::AstroCoords.AstroCoord,
    tof::Number,
)
    r1 = Cartesian(coord1, μ)[1:3]
    r2 = Cartesian(coord2, μ)[1:3]
    return LambertProblem(μ, r1, r2, tof)
end

"""
    LambertSolution{VT1, VT2}

Solution to a Lambert problem.

# Fields
- `v1::VT1`: Initial velocity vector (preserves input vector type)
- `v2::VT2`: Final velocity vector (preserves input vector type)
- `numiter::Int`: Number of iterations used
- `retcode::Symbol`: Return code (:SUCCESS, :MAXIMUM_ITERATIONS, etc.)
"""
struct LambertSolution{VT1,VT2} <: AbstractAstroSolution
    v1::VT1
    v2::VT2
    numiter::Int
    retcode::Symbol
end

function SciMLBase.remake(
    prob::LambertProblem;
    μ = nothing,
    r1 = nothing,
    r2 = nothing,
    tof = nothing,
)
    return LambertProblem(
        μ === nothing ? prob.μ : μ,
        r1 === nothing ? prob.r1 : r1,
        r2 === nothing ? prob.r2 : r2,
        tof === nothing ? prob.tof : tof,
    )
end

# Iterator struct for SciMLBase init/solve! interface
struct LambertIterator{P<:LambertProblem,A<:AbstractLambertSolver,K<:NamedTuple}
    prob::P
    alg::A
    kwargs::K
end

# SciMLBase init function - initializes the iterator
function SciMLBase.init(prob::LambertProblem, alg::AbstractLambertSolver; kwargs...)
    return LambertIterator(prob, alg, values(kwargs))
end

# Default algorithm version of init with heuristic selection
function SciMLBase.init(
    prob::LambertProblem;
    alg = nothing,
    M::Int = 0,
    prograde::Bool = true,
    kwargs...,
)
    if alg === nothing
        # Use heuristic algorithm selection based on problem characteristics
        alg = select_lambert_algorithm(prob, M, prograde)
    end
    return SciMLBase.init(prob, alg; kwargs...)
end

# SciMLBase solve! function - performs the computation
function SciMLBase.solve!(iterator::LambertIterator)
    # Apply kwargs to the algorithm if they correspond to solver parameters
    alg_with_kwargs = _apply_kwargs_to_solver(iterator.alg, iterator.kwargs)

    # Call the specific solver method (existing SciMLBase.solve methods)
    return SciMLBase.solve(iterator.prob, alg_with_kwargs)
end

# Convenience solve method following SciMLBase pattern
function solve(prob::LambertProblem, alg::AbstractLambertSolver; kwargs...)
    return SciMLBase.solve!(SciMLBase.init(prob, alg; kwargs...))
end

# Convenience solve method with heuristic algorithm selection
function solve(
    prob::LambertProblem;
    alg = nothing,
    M::Int = 0,
    prograde::Bool = true,
    kwargs...,
)
    return SciMLBase.solve!(
        SciMLBase.init(prob; alg = alg, M = M, prograde = prograde, kwargs...),
    )
end

# Helper function to apply kwargs to solver parameters
function _apply_kwargs_to_solver(alg::AbstractLambertSolver, kwargs::NamedTuple)
    if isempty(kwargs)
        return alg
    end

    # Get the current algorithm parameters
    alg_params = Dict()
    for field in fieldnames(typeof(alg))
        alg_params[field] = getfield(alg, field)
    end

    # Update with provided kwargs that match solver parameters
    for (key, value) in pairs(kwargs)
        if haskey(alg_params, key)
            alg_params[key] = value
        end
    end

    # Reconstruct the algorithm with updated parameters
    return typeof(alg)(; alg_params...)
end
