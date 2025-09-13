module Lambert

using AstroCoords
using LinearAlgebra
using Parameters
using Roots
using SciMLBase
using StaticArraysCore
import SciMLBase: solve, remake, init, solve!

# Abstract solver type
abstract type AbstractLambertSolver <: SciMLBase.AbstractSciMLAlgorithm end

include("lambert_problem.jl")
include("utils.jl")

# Solvers
include("arora_solver.jl")
include("avanzini_solver.jl")
include("battin_solver.jl")
include("gauss_solver.jl")
include("gooding_solver.jl")
include("izzo_solver.jl")
include("vallado_solver.jl")

end
