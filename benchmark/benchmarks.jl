using BenchmarkTools
using Lambert
using LinearAlgebra
using SciMLBase
using StaticArraysCore

const SUITE = BenchmarkGroup()

# ---------------------------------------------------------------------------
# Shared test geometry (Vector for direct calls, LambertProblem converts internally)
# ---------------------------------------------------------------------------

const μ_earth = 3.986004418e5  # km³/s²

# Vallado Example 5.7 — moderate transfer angle (~50°), works for all solvers
const r1_vallado = [15945.34, 0.0, 0.0]
const r2_vallado = [12214.83899, 10249.46731, 0.0]
const tof_vallado = 76.0 * 60.0

# Curtis Example 5.2 — large transfer angle (~160°), robust solvers only
const r1_curtis = [5000.0, 10000.0, 2100.0]
const r2_curtis = [-14600.0, 2500.0, 7000.0]
const tof_curtis = 3600.0

# GMAT Hyperbolic — high-energy escape trajectory, robust solvers only
const r1_hyper = [7100.0, 200.0, 1300.0]
const r2_hyper = [-38113.5870, 67274.1946, 29309.5799]
const tof_hyper = 12000.0

# Near-180° transfer — challenging geometry, robust solvers only
const r1_180 = [7000.0, 0.0, 0.0]
const r2_180 = [-6900.0, 500.0, 0.0]
const tof_180 = 5400.0

# Multi-revolution geometry
const r1_multi = r1_vallado
const r2_multi = r2_vallado
const tof_multi_m1 = tof_vallado * 5.0
const tof_multi_m2 = tof_vallado * 7.0

# ---------------------------------------------------------------------------
# Zero-revolution benchmarks — Vallado case (all solvers)
# ---------------------------------------------------------------------------

SUITE["zero_rev"] = BenchmarkGroup(["single_rev", "solver"])

const ALL_SOLVER_PAIRS = [
    ("Gooding", GoodingSolver()),
    ("Izzo", IzzoSolver()),
    ("Battin", BattinSolver()),
    ("Gauss", GaussSolver()),
    ("Vallado", ValladoSolver()),
    ("Arora", AroraSolver()),
    ("Avanzini", AvanziniSolver()),
    ("Russell", RussellSolver()),
]

for (name, solver) in ALL_SOLVER_PAIRS
    SUITE["zero_rev"][name] = @benchmarkable solve(prob, $solver) setup =
        (prob = LambertProblem($μ_earth, $r1_vallado, $r2_vallado, $tof_vallado))
end

# ---------------------------------------------------------------------------
# Challenging geometry benchmarks — robust solvers only
# ---------------------------------------------------------------------------

SUITE["challenging"] = BenchmarkGroup(["challenging", "solver"])

const ROBUST_SOLVER_PAIRS = [
    ("Gooding", GoodingSolver()),
    ("Izzo", IzzoSolver()),
    ("Battin", BattinSolver()),
    ("Vallado", ValladoSolver()),
    ("Russell", RussellSolver()),
]

for (name, solver) in ROBUST_SOLVER_PAIRS
    SUITE["challenging"][name] = BenchmarkGroup()

    SUITE["challenging"][name]["curtis"] = @benchmarkable solve(prob, $solver) setup =
        (prob = LambertProblem($μ_earth, $r1_curtis, $r2_curtis, $tof_curtis))

    SUITE["challenging"][name]["hyperbolic"] = @benchmarkable solve(prob, $solver) setup =
        (prob = LambertProblem($μ_earth, $r1_hyper, $r2_hyper, $tof_hyper))

    SUITE["challenging"][name]["near_180"] = @benchmarkable solve(prob, $solver) setup =
        (prob = LambertProblem($μ_earth, $r1_180, $r2_180, $tof_180))
end

# ---------------------------------------------------------------------------
# Multi-revolution benchmarks — capable solvers only
# ---------------------------------------------------------------------------

SUITE["multi_rev"] = BenchmarkGroup(["multi_rev", "solver"])

# Solvers with low_path parameter
const MULTIREV_WITH_LOWPATH = [("Izzo", IzzoSolver), ("Russell", RussellSolver)]

for (name, SolverType) in MULTIREV_WITH_LOWPATH
    SUITE["multi_rev"][name] = BenchmarkGroup()

    SUITE["multi_rev"][name]["M1_low"] = @benchmarkable solve(prob, solver) setup = (
        prob = LambertProblem($μ_earth, $r1_multi, $r2_multi, $tof_multi_m1);
        solver = $SolverType(M = 1, low_path = true)
    )

    SUITE["multi_rev"][name]["M2_low"] = @benchmarkable solve(prob, solver) setup = (
        prob = LambertProblem($μ_earth, $r1_multi, $r2_multi, $tof_multi_m2);
        solver = $SolverType(M = 2, low_path = true)
    )
end

# Solvers without low_path parameter
const MULTIREV_NO_LOWPATH =
    [("Battin", BattinSolver), ("Gauss", GaussSolver), ("Arora", AroraSolver)]

for (name, SolverType) in MULTIREV_NO_LOWPATH
    SUITE["multi_rev"][name] = BenchmarkGroup()

    SUITE["multi_rev"][name]["M1"] = @benchmarkable solve(prob, solver) setup = (
        prob = LambertProblem($μ_earth, $r1_multi, $r2_multi, $tof_multi_m1);
        solver = $SolverType(M = 1)
    )

    SUITE["multi_rev"][name]["M2"] = @benchmarkable solve(prob, solver) setup = (
        prob = LambertProblem($μ_earth, $r1_multi, $r2_multi, $tof_multi_m2);
        solver = $SolverType(M = 2)
    )
end

# ---------------------------------------------------------------------------
# Retrograde benchmarks — robust solvers
# ---------------------------------------------------------------------------

SUITE["retrograde"] = BenchmarkGroup(["retrograde", "solver"])

const RETROGRADE_SOLVER_TYPES = [
    ("Gooding", GoodingSolver),
    ("Izzo", IzzoSolver),
    ("Battin", BattinSolver),
    ("Vallado", ValladoSolver),
    ("Russell", RussellSolver),
]

for (name, SolverType) in RETROGRADE_SOLVER_TYPES
    SUITE["retrograde"][name] = @benchmarkable solve(prob, solver) setup = (
        prob = LambertProblem($μ_earth, $r1_hyper, $r2_hyper, $tof_hyper);
        solver = $SolverType(prograde = false)
    )
end

# ---------------------------------------------------------------------------
# Russell solver — convergence order comparison
# ---------------------------------------------------------------------------

SUITE["russell_order"] = BenchmarkGroup(["russell", "convergence"])

for order = 1:3
    label = order == 1 ? "Newton" : order == 2 ? "Halley" : "3rd_order"
    solver = RussellSolver(order = order)

    SUITE["russell_order"]["$(label)_vallado"] = @benchmarkable solve(prob, $solver) setup =
        (prob = LambertProblem($μ_earth, $r1_vallado, $r2_vallado, $tof_vallado))

    SUITE["russell_order"]["$(label)_curtis"] = @benchmarkable solve(prob, $solver) setup =
        (prob = LambertProblem($μ_earth, $r1_curtis, $r2_curtis, $tof_curtis))
end

# ---------------------------------------------------------------------------
# Direct function call benchmarks (no SciMLBase overhead)
# ---------------------------------------------------------------------------

SUITE["direct_call"] = BenchmarkGroup(["direct", "raw"])

SUITE["direct_call"]["gooding1990"] =
    @benchmarkable gooding1990($μ_earth, $r1_vallado, $r2_vallado, $tof_vallado)

SUITE["direct_call"]["izzo2015"] =
    @benchmarkable izzo2015($μ_earth, $r1_vallado, $r2_vallado, $tof_vallado)

SUITE["direct_call"]["battin1984"] =
    @benchmarkable battin1984($μ_earth, $r1_vallado, $r2_vallado, $tof_vallado)

SUITE["direct_call"]["gauss1809"] =
    @benchmarkable gauss1809($μ_earth, $r1_vallado, $r2_vallado, $tof_vallado)

SUITE["direct_call"]["vallado2013"] =
    @benchmarkable vallado2013($μ_earth, $r1_vallado, $r2_vallado, $tof_vallado)

SUITE["direct_call"]["arora2013"] =
    @benchmarkable arora2013($μ_earth, $r1_vallado, $r2_vallado, $tof_vallado)

SUITE["direct_call"]["avanzini2008"] =
    @benchmarkable avanzini2008($μ_earth, $r1_vallado, $r2_vallado, $tof_vallado)

SUITE["direct_call"]["russell2021"] =
    @benchmarkable russell2021($μ_earth, $r1_vallado, $r2_vallado, $tof_vallado)

# ---------------------------------------------------------------------------
# Porkchop grid benchmarks
# ---------------------------------------------------------------------------

SUITE["porkchop"] = BenchmarkGroup(["porkchop", "grid"])

# Circular coplanar orbits (simplified Earth–Mars analogy)
const R_inner = 1.0e8   # km  (≈ 1 AU simplified)
const R_outer = 1.524e8  # km  (≈ Mars orbit)
const n_inner = sqrt(μ_earth / R_inner^3)
const n_outer = sqrt(μ_earth / R_outer^3)

function _inner_state(t)
    θ = n_inner * t
    sθ, cθ = sincos(θ)
    v_mag = R_inner * n_inner
    return SVector{6}(R_inner * cθ, R_inner * sθ, 0.0, -v_mag * sθ, v_mag * cθ, 0.0)
end

function _outer_state(t)
    θ = n_outer * t
    sθ, cθ = sincos(θ)
    v_mag = R_outer * n_outer
    return SVector{6}(R_outer * cθ, R_outer * sθ, 0.0, -v_mag * sθ, v_mag * cθ, 0.0)
end

# Hohmann-like synodic period for grid sizing
const T_synodic = 2π / abs(n_inner - n_outer)
const dep_range = range(0.0, T_synodic, length = 20)
const arr_range = range(0.3 * T_synodic, 1.5 * T_synodic, length = 20)

# 20×20 grid — small (400 Lambert solves)
SUITE["porkchop"]["20x20_Izzo"] = @benchmarkable porkchop_grid(
    $μ_earth,
    $_inner_state,
    $_outer_state,
    $dep_range,
    $arr_range;
    solver = IzzoSolver(),
)

SUITE["porkchop"]["20x20_Russell"] = @benchmarkable porkchop_grid(
    $μ_earth,
    $_inner_state,
    $_outer_state,
    $dep_range,
    $arr_range;
    solver = RussellSolver(),
)

SUITE["porkchop"]["20x20_Gooding"] = @benchmarkable porkchop_grid(
    $μ_earth,
    $_inner_state,
    $_outer_state,
    $dep_range,
    $arr_range;
    solver = GoodingSolver(),
)

# 50×50 grid — medium (2500 Lambert solves)
const dep_range_50 = range(0.0, T_synodic, length = 50)
const arr_range_50 = range(0.3 * T_synodic, 1.5 * T_synodic, length = 50)

SUITE["porkchop"]["50x50_Izzo"] = @benchmarkable porkchop_grid(
    $μ_earth,
    $_inner_state,
    $_outer_state,
    $dep_range_50,
    $arr_range_50;
    solver = IzzoSolver(),
)

SUITE["porkchop"]["50x50_Russell"] = @benchmarkable porkchop_grid(
    $μ_earth,
    $_inner_state,
    $_outer_state,
    $dep_range_50,
    $arr_range_50;
    solver = RussellSolver(),
)

SUITE["porkchop"]["50x50_Gooding"] = @benchmarkable porkchop_grid(
    $μ_earth,
    $_inner_state,
    $_outer_state,
    $dep_range_50,
    $arr_range_50;
    solver = GoodingSolver(),
)

# ---------------------------------------------------------------------------
# Tune / load parameters
# ---------------------------------------------------------------------------

paramspath = joinpath(dirname(@__FILE__), "params.json")

if isfile(paramspath)
    loadparams!(SUITE, BenchmarkTools.load(paramspath)[1], :evals)
else
    tune!(SUITE)
    BenchmarkTools.save(paramspath, BenchmarkTools.params(SUITE))
end
