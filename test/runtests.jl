using Test
using AstroCoords
using Lambert
using LinearAlgebra
using SciMLBase
using StaticArraysCore
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
    RussellSolver(),
]

# Robust solvers (excluding those with known accuracy issues)
const ROBUST_SOLVERS =
    [GoodingSolver(), IzzoSolver(), BattinSolver(), ValladoSolver(), RussellSolver()]

const MULTIREV_SOLVERS =
    [IzzoSolver(), BattinSolver(), GaussSolver(), AroraSolver(), RussellSolver()]

@testset "AstroProblemsLambert.jl Tests" begin
    # Run existing Lambert solver tests
    include("lambert_solvers_tests.jl")

    # Run new EnsembleProblem integration tests  
    include("ensemble_tests.jl")

    # Run porkchop grid tests (core functionality, no plotting dependencies)
    include("test_porkchop_grid.jl")

    # Run porkchop plot tests (requires Plots.jl extension)
    include("test_porkchop_plot.jl")

    # Run performance and allocation tests
    include("test_performance.jl")
end


# Differentiability tests — gated behind LAMBERT_TEST_DIFF environment variable.
#   "true" / "all"       → all backends (ForwardDiff, Enzyme, Mooncake, PolyesterForwardDiff, Zygote)
#   "ForwardDiff"        → ForwardDiff only (fast smoke-test)
#   comma-separated list → specific backends
#   unset / "false"      → skip
const _DIFF_ENV = get(ENV, "LAMBERT_TEST_DIFF", "false")

if _DIFF_ENV ∉ ("false", "")
    using DifferentiationInterface
    using FiniteDiff
    using ForwardDiff  # always loaded — used as reference for analytical Jacobian tests

    _run_all = _DIFF_ENV ∈ ("true", "all")
    _requested = _run_all ? Set{String}() : Set(strip.(split(_DIFF_ENV, ",")))
    _need(name) = _run_all || name ∈ _requested

    _backend_list = Tuple{String,Any}[]

    if _need("ForwardDiff")
        push!(_backend_list, ("ForwardDiff", DifferentiationInterface.AutoForwardDiff()))
    end
    if _need("Enzyme")
        using Enzyme
        push!(
            _backend_list,
            (
                "Enzyme",
                DifferentiationInterface.AutoEnzyme(;
                    mode = Enzyme.set_runtime_activity(Enzyme.Forward),
                ),
            ),
        )
    end
    if _need("Mooncake")
        using Mooncake
        push!(
            _backend_list,
            ("Mooncake", DifferentiationInterface.AutoMooncake(; config = nothing)),
        )
    end
    if _need("PolyesterForwardDiff")
        using PolyesterForwardDiff
        push!(
            _backend_list,
            ("PolyesterForwardDiff", DifferentiationInterface.AutoPolyesterForwardDiff()),
        )
    end
    if _need("Zygote")
        using Zygote
        push!(_backend_list, ("Zygote", DifferentiationInterface.AutoZygote()))
    end

    if isempty(_backend_list)
        error(
            "LAMBERT_TEST_DIFF=\"$_DIFF_ENV\" did not match any backend. " *
            "Valid names: ForwardDiff, Enzyme, Mooncake, PolyesterForwardDiff, Zygote",
        )
    end

    const _BACKENDS = Tuple(_backend_list)

    @info "Running differentiability tests with backends: $(join([b[1] for b in _BACKENDS], ", "))"

    @testset "Differentiability" begin
        include("differentiability/test_differentiability.jl")
    end
else
    @info "Skipping differentiability tests (set LAMBERT_TEST_DIFF to enable)"
end

@testset "Aqua Tests" begin
    Aqua.test_all(
        Lambert;
        ambiguities = (recursive = false),
        deps_compat = (check_extras = false),
    )
end

@testset "JET Testing" begin
    rep = JET.test_package(
        Lambert;
        toplevel_logger = nothing,
        target_modules = (@__MODULE__,),
    )
end
