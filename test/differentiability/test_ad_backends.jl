# AD backend tests for Lambert.jl
#
# Validates that each AD backend produces Jacobians matching ForwardDiff
# (the reference) when differentiating through Lambert.solve().
#
# ChainRulesCore is loaded BEFORE this file, so the LambertChainRulesCoreExt
# extension is active. Reverse-mode backends (Zygote, Mooncake) use the
# extension's rrule; ForwardDiff still uses Dual number propagation.
#
# Expects:
#   _BACKENDS — Tuple of (name::String, backend) pairs (defined in runtests.jl)
#   _lambert_map, _DIFF_X0, _DIFF_SOLVERS, _REF_BACKEND — from test_analytical_jacobian.jl

# AD-compatible solvers (excludes Avanzini/Vallado which use Roots.jl)
const _AD_SOLVERS = [
    ("Izzo", IzzoSolver()),
    ("Gooding", GoodingSolver()),
    ("Russell", RussellSolver()),
    ("Battin", BattinSolver()),
    ("Arora", AroraSolver()),
    ("Gauss", GaussSolver()),
]

for (bname, backend) in _BACKENDS
    @testset "$bname" begin
        for (sname, solver) in _AD_SOLVERS
            @testset "$sname" begin
                J_ref = jacobian(x -> _lambert_map(x, solver), _REF_BACKEND, _DIFF_X0)
                J_ad = jacobian(x -> _lambert_map(x, solver), backend, _DIFF_X0)
                @test J_ad ≈ J_ref rtol = 1e-6 atol = 1e-14
            end
        end
    end
end
