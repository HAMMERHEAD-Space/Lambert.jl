# Analytical Jacobian tests for Lambert.jl
#
# Validates lambert_jacobian() against ForwardDiff (which differentiates through
# the actual solver code via Dual numbers, NOT through ChainRulesCore rules).
#
# This file runs BEFORE ChainRulesCore is loaded, so there is no extension
# intercepting the solve() call — ForwardDiff sees the real computation.
#
# Expects: DifferentiationInterface, FiniteDiff, ForwardDiff — already loaded

using DifferentiationInterface: jacobian, AutoFiniteDiff, AutoForwardDiff

const _REF_BACKEND = AutoForwardDiff()

# Test problem: Earth orbit transfer
const _DIFF_MU = 3.986004418e5
const _DIFF_R1 = SVector{3}(15945.34, 0.0, 0.0)
const _DIFF_R2 = SVector{3}(12214.83899, 10249.46731, 0.0)
const _DIFF_TOF = 76.0 * 60.0

"""
Map input vector x = [μ, r1..., r2..., tof] to output [v1..., v2...] for a given solver.
"""
function _lambert_map(x::AbstractVector{T}, solver) where {T}
    prob = LambertProblem(
        x[1],
        SVector{3,T}(x[2], x[3], x[4]),
        SVector{3,T}(x[5], x[6], x[7]),
        x[8],
    )
    sol = Lambert.solve(prob, solver)
    return Vector(vcat(sol.v1, sol.v2))
end

const _DIFF_X0 = Vector(vcat(_DIFF_MU, _DIFF_R1, _DIFF_R2, _DIFF_TOF))

# Solvers whose internal accuracy matches the STM's Kepler propagation closely
const _STM_TIGHT_SOLVERS = [
    ("Izzo", IzzoSolver()),
    ("Gooding", GoodingSolver()),
    ("Russell", RussellSolver()),
    ("Arora", AroraSolver()),
]

# Solvers with lower internal accuracy (Battin, Gauss use different parameterizations)
const _STM_LOOSE_SOLVERS = [("Battin", BattinSolver()), ("Gauss", GaussSolver())]

# Solvers that use Roots.jl internally (not AD-compatible)
const _STM_ROOTFINDER_SOLVERS =
    [("Avanzini", AvanziniSolver()), ("Vallado", ValladoSolver())]

function _assemble_jacobian(J_an)
    return hcat(
        vcat(J_an.dv1_dmu, J_an.dv2_dmu),
        vcat(J_an.dv1_dr1, J_an.dv2_dr1),
        vcat(J_an.dv1_dr2, J_an.dv2_dr2),
        vcat(J_an.dv1_dtof, J_an.dv2_dtof),
    )
end

for (name, solver) in _STM_TIGHT_SOLVERS
    @testset "$name" begin
        prob = LambertProblem(_DIFF_MU, _DIFF_R1, _DIFF_R2, _DIFF_TOF)
        sol = Lambert.solve(prob, solver)
        J_an = _assemble_jacobian(lambert_jacobian(prob, sol))
        J_ref = jacobian(x -> _lambert_map(x, solver), _REF_BACKEND, _DIFF_X0)
        @test J_an ≈ J_ref rtol = 1e-6 atol = 1e-14
    end
end
for (name, solver) in _STM_LOOSE_SOLVERS
    @testset "$name" begin
        prob = LambertProblem(_DIFF_MU, _DIFF_R1, _DIFF_R2, _DIFF_TOF)
        sol = Lambert.solve(prob, solver)
        J_an = _assemble_jacobian(lambert_jacobian(prob, sol))
        J_ref = jacobian(x -> _lambert_map(x, solver), _REF_BACKEND, _DIFF_X0)
        @test J_an ≈ J_ref rtol = 1e-2 atol = 1e-10
    end
end
# Avanzini/Vallado use Roots.jl (not AD-compatible, FiniteDiff unreliable through
# bisection). Cross-check: their STM should match a high-accuracy solver's STM.
J_ref_izzo = let
    prob = LambertProblem(_DIFF_MU, _DIFF_R1, _DIFF_R2, _DIFF_TOF)
    sol = Lambert.solve(prob, IzzoSolver())
    _assemble_jacobian(lambert_jacobian(prob, sol))
end
for (name, solver) in _STM_ROOTFINDER_SOLVERS
    @testset "$name" begin
        prob = LambertProblem(_DIFF_MU, _DIFF_R1, _DIFF_R2, _DIFF_TOF)
        sol = Lambert.solve(prob, solver)
        J_an = _assemble_jacobian(lambert_jacobian(prob, sol))
        @test all(isfinite, J_an)
        @test J_an ≈ J_ref_izzo rtol = 1e-2 atol = 1e-10
    end
end
