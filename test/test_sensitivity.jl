@testset "Analytical Sensitivity" begin
    # Reference problem used throughout
    μ = 3.986004418e5
    r1 = [15945.34, 0.0, 0.0]
    r2 = [12214.83899, 10249.46731, 0.0]
    tof = 76.0 * 60

    prob = LambertProblem(μ, r1, r2, tof)
    sol = solve(prob, IzzoSolver())

    @testset "lambert_jacobian types and structure" begin
        jac = lambert_jacobian(prob, sol)

        @test jac.dv1_dr1 isa SMatrix{3,3}
        @test jac.dv1_dr2 isa SMatrix{3,3}
        @test jac.dv1_dtof isa SVector{3}
        @test jac.dv1_dmu isa SVector{3}
        @test jac.dv2_dr1 isa SMatrix{3,3}
        @test jac.dv2_dr2 isa SMatrix{3,3}
        @test jac.dv2_dtof isa SVector{3}
        @test jac.dv2_dmu isa SVector{3}
    end

    @testset "lambert_jacobian cross-solver consistency" begin
        jac_izzo = lambert_jacobian(prob, sol)

        for solver in [GoodingSolver(), RussellSolver(), McElreathSolver()]
            other_sol = solve(prob, solver)
            jac_other = lambert_jacobian(prob, other_sol)
            @test jac_other.dv1_dr1 ≈ jac_izzo.dv1_dr1 atol = 1e-4
            @test jac_other.dv1_dr2 ≈ jac_izzo.dv1_dr2 atol = 1e-4
            @test jac_other.dv2_dtof ≈ jac_izzo.dv2_dtof atol = 1e-4
            @test jac_other.dv1_dmu ≈ jac_izzo.dv1_dmu atol = 1e-6
        end
    end

    @testset "lambert_jacobian finite-difference validation" begin
        jac = lambert_jacobian(prob, sol)

        ε = 1.0
        for i = 1:3
            r1p = copy(r1)
            r1p[i] += ε
            sol_p = solve(LambertProblem(μ, r1p, r2, tof), IzzoSolver())
            fd_col = (sol_p.v1 - sol.v1) / ε
            @test fd_col ≈ jac.dv1_dr1[:, i] rtol = 1e-3
        end

        ε_tof = 0.1
        sol_tp = solve(LambertProblem(μ, r1, r2, tof + ε_tof), IzzoSolver())
        fd_tof = (sol_tp.v1 - sol.v1) / ε_tof
        @test fd_tof ≈ jac.dv1_dtof rtol = 1e-3
    end

    @testset "_twobody_stm symplectic property" begin
        r1_sv = SVector{3}(7000.0, 0.0, 0.0)
        v1_sv = SVector{3}(0.0, 7.5, 0.0)
        Φrr, Φrv, Φvr, Φvv = Lambert._twobody_stm(r1_sv, v1_sv, 3600.0, 398600.0)

        I3 = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
        @test Φrr * Φvv' - Φrv * Φvr' ≈ I3 atol = 1e-6
    end

    @testset "_solve_kepler_universal satisfies Kepler equation" begin
        r1_mag = 7000.0
        r1dv1 = 100.0
        α = 1e-4
        sqrtmu = √398600.0
        dt = 3600.0

        χ = Lambert._solve_kepler_universal(r1_mag, r1dv1, α, sqrtmu, dt)
        ψ = α * χ^2
        C2 = Lambert.c2(ψ)
        C3 = Lambert.c3(ψ)
        F =
            (r1dv1 / sqrtmu) * χ^2 * C2 + r1_mag * χ * (1.0 - ψ * C3) + χ^3 * C3 -
            sqrtmu * dt
        @test abs(F) < 1e-8
        @test χ > 0
    end

    @testset "_mu_orbit_sensitivity has expected scale" begin
        r1_sv = SVector{3}(7000.0, 0.0, 0.0)
        v1_sv = SVector{3}(0.0, 7.5, 0.0)
        Φrr, Φrv, Φvr, Φvv = Lambert._twobody_stm(r1_sv, v1_sv, 3600.0, 398600.0)
        dr2, dv2 = Lambert._mu_orbit_sensitivity(
            r1_sv,
            v1_sv,
            3600.0,
            398600.0,
            Φrr,
            Φrv,
            Φvr,
            Φvv,
        )

        @test norm(dr2) > 0
        @test norm(dv2) > 0
        @test norm(dr2) < 1.0
        @test norm(dv2) < 1e-3
    end
end
