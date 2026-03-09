module LambertChainRulesCoreExt

using Lambert
using ChainRulesCore
using LinearAlgebra
using StaticArraysCore

import Lambert: LambertProblem, LambertSolution, AbstractLambertSolver, lambert_jacobian

"""
Forward-mode rule for `solve(prob::LambertProblem, solver::AbstractLambertSolver)`.

Propagates tangent vectors through the analytical Jacobian computed via the
two-body state transition matrix.
"""
function ChainRulesCore.frule(
    (_, Δprob, _),
    ::typeof(Lambert.solve),
    prob::LambertProblem,
    solver::AbstractLambertSolver;
    kwargs...,
)
    sol = Lambert.solve(prob, solver; kwargs...)
    J = lambert_jacobian(prob, sol)

    Δμ = Δprob isa ChainRulesCore.Tangent ? Δprob.μ : ZeroTangent()
    Δr1 = Δprob isa ChainRulesCore.Tangent ? Δprob.r1 : ZeroTangent()
    Δr2 = Δprob isa ChainRulesCore.Tangent ? Δprob.r2 : ZeroTangent()
    Δtof = Δprob isa ChainRulesCore.Tangent ? Δprob.tof : ZeroTangent()

    Δv1 = ZeroTangent()
    Δv2 = ZeroTangent()

    if !(Δμ isa AbstractZero)
        Δv1 = Δv1 + J.dv1_dmu * Δμ
        Δv2 = Δv2 + J.dv2_dmu * Δμ
    end
    if !(Δr1 isa AbstractZero)
        Δv1 = Δv1 + J.dv1_dr1 * SVector{3}(Δr1[1], Δr1[2], Δr1[3])
        Δv2 = Δv2 + J.dv2_dr1 * SVector{3}(Δr1[1], Δr1[2], Δr1[3])
    end
    if !(Δr2 isa AbstractZero)
        Δv1 = Δv1 + J.dv1_dr2 * SVector{3}(Δr2[1], Δr2[2], Δr2[3])
        Δv2 = Δv2 + J.dv2_dr2 * SVector{3}(Δr2[1], Δr2[2], Δr2[3])
    end
    if !(Δtof isa AbstractZero)
        Δv1 = Δv1 + J.dv1_dtof * Δtof
        Δv2 = Δv2 + J.dv2_dtof * Δtof
    end

    Δsol = Tangent{LambertSolution}(;
        v1 = Δv1,
        v2 = Δv2,
        numiter = NoTangent(),
        retcode = NoTangent(),
    )
    return sol, Δsol
end

"""
Reverse-mode rule for `solve(prob::LambertProblem, solver::AbstractLambertSolver)`.

The pullback computes `Jᵀ · d̄sol` efficiently using the analytical Jacobian.
"""
function ChainRulesCore.rrule(
    ::typeof(Lambert.solve),
    prob::LambertProblem,
    solver::AbstractLambertSolver;
    kwargs...,
)
    sol = Lambert.solve(prob, solver; kwargs...)
    J = lambert_jacobian(prob, sol)

    function solve_pullback(d̄sol)
        d̄v1 = d̄sol.v1
        d̄v2 = d̄sol.v2

        # Transpose multiplication: d̄input = Jᵀ [d̄v1; d̄v2]
        sv1 =
            d̄v1 isa AbstractZero ? SVector{3}(0.0, 0.0, 0.0) :
            SVector{3}(d̄v1[1], d̄v1[2], d̄v1[3])
        sv2 =
            d̄v2 isa AbstractZero ? SVector{3}(0.0, 0.0, 0.0) :
            SVector{3}(d̄v2[1], d̄v2[2], d̄v2[3])

        d̄μ = dot(J.dv1_dmu, sv1) + dot(J.dv2_dmu, sv2)
        d̄r1 = J.dv1_dr1' * sv1 + J.dv2_dr1' * sv2
        d̄r2 = J.dv1_dr2' * sv1 + J.dv2_dr2' * sv2
        d̄tof = dot(J.dv1_dtof, sv1) + dot(J.dv2_dtof, sv2)

        d̄prob = Tangent{LambertProblem}(; μ = d̄μ, r1 = d̄r1, r2 = d̄r2, tof = d̄tof)
        return NoTangent(), d̄prob, NoTangent()
    end

    return sol, solve_pullback
end

end # module
