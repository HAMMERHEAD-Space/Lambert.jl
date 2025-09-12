using Test
using AstroCoords
using Lambert
using LinearAlgebra
using SciMLBase
using Random
using Statistics

using AllocCheck
using Aqua
using JET

# List of all solvers to test
const ALL_SOLVERS = [
    GoodingSolver(),
    IzzoSolver(),
    BattinSolver(),
    GaussSolver(),
    ValladoSolver(),
    AroraSolver(),
    AvanziniSolver(),
]

# Robust solvers (excluding those with known accuracy issues)
const ROBUST_SOLVERS = [
    GoodingSolver(),
    IzzoSolver(),
    BattinSolver(),
    ValladoSolver(),
]

const MULTIREV_SOLVERS = [
    IzzoSolver(),
    BattinSolver(),
    GaussSolver(),
    AroraSolver(),
]

@testset "AstroProblemsLambert.jl Tests" begin
    # Run existing Lambert solver tests
    include("lambert_solvers_tests.jl")
    
    # Run new EnsembleProblem integration tests  
    include("ensemble_tests.jl")
    
    # Run performance and allocation tests
    #include("test_performance.jl")
end


@testset "Aqua Tests" begin
    Aqua.test_all(AstroProblemsLambert;
        ambiguities=(recursive = false),
        deps_compat=(check_extras = false),
    )
end

@testset "JET Testing" begin
    rep = JET.test_package(
        AstroProblemsLambert; toplevel_logger=nothing, target_modules=(@__MODULE__,)
    )
end
