@testset "Utility Functions" begin
    # ── Stumpff functions ─────────────────────────────────────────────────
    @testset "Stumpff c2 branches" begin
        @test Lambert.c2(1.0) ≈ (1 - cos(1.0)) / 1.0
        @test Lambert.c2(-1.0) ≈ (cosh(1.0) - 1) / 1.0
        @test Lambert.c2(1e-8) ≈ 0.5 atol = 1e-6
        @test Lambert.c2(-1e-8) ≈ 0.5 atol = 1e-6
        @test Lambert.c2(0.0) ≈ 0.5 atol = 1e-6
    end

    @testset "Stumpff c3 branches" begin
        @test Lambert.c3(1.0) ≈ (1.0 - sin(1.0)) / 1.0
        @test Lambert.c3(-1.0) ≈ (sinh(1.0) - 1.0) / 1.0
        @test Lambert.c3(1e-8) ≈ 1 / 6 atol = 1e-6
        @test Lambert.c3(-1e-8) ≈ 1 / 6 atol = 1e-6
        @test Lambert.c3(0.0) ≈ 1 / 6 atol = 1e-6
    end

    # ── Validation ────────────────────────────────────────────────────────
    @testset "assert_parameters_are_valid" begin
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 7000.0, 0.0]

        @test_throws AssertionError Lambert.assert_parameters_are_valid(
            -1.0,
            r1,
            r2,
            3600.0,
            0,
        )
        @test_throws AssertionError Lambert.assert_parameters_are_valid(
            0.0,
            r1,
            r2,
            3600.0,
            0,
        )
        @test_throws AssertionError Lambert.assert_parameters_are_valid(
            398600.0,
            [0.0, 0.0, 0.0],
            r2,
            3600.0,
            0,
        )
        @test_throws AssertionError Lambert.assert_parameters_are_valid(
            398600.0,
            r1,
            [0.0, 0.0, 0.0],
            3600.0,
            0,
        )
        @test_throws AssertionError Lambert.assert_parameters_are_valid(
            398600.0,
            r1,
            r2,
            -1.0,
            0,
        )
        @test_throws AssertionError Lambert.assert_parameters_are_valid(
            398600.0,
            r1,
            r2,
            0.0,
            0,
        )
        @test_throws AssertionError Lambert.assert_parameters_are_valid(
            398600.0,
            r1,
            r2,
            3600.0,
            -1,
        )
        @test_throws AssertionError Lambert.assert_parameters_are_valid(
            398600.0,
            r1,
            r1,
            3600.0,
            0,
        )
        Lambert.assert_parameters_are_valid(398600.0, r1, r2, 3600.0, 0)
    end

    @testset "assert_transfer_angle_not_zero" begin
        @test_throws AssertionError Lambert.assert_transfer_angle_not_zero(0.0)
        @test_throws AssertionError Lambert.assert_transfer_angle_not_zero(1e-15)
        Lambert.assert_transfer_angle_not_zero(0.1)
    end

    @testset "assert_transfer_angle_not_pi" begin
        @test_throws AssertionError Lambert.assert_transfer_angle_not_pi(π)
        Lambert.assert_transfer_angle_not_pi(3.0)
    end

    @testset "validate_multi_revolution" begin
        Lambert.validate_multi_revolution(0, "TestSolver")
        @test_throws ErrorException Lambert.validate_multi_revolution(1, "TestSolver")
    end

    # ── Geometry ──────────────────────────────────────────────────────────
    @testset "get_transfer_angle" begin
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 7000.0, 0.0]
        @test Lambert.get_transfer_angle(r1, r2, true) ≈ π / 2 atol = 1e-10
        @test Lambert.get_transfer_angle(r1, r2, false) ≈ 2π - π / 2 atol = 1e-10

        r1_col = [7000.0, 0.0, 0.0]
        r2_col = [14000.0, 0.0, 0.0]
        @test Lambert.get_transfer_angle(r1_col, r2_col, true) ≈ 0.0 atol = 1e-10

        r2_opp = [-7000.0, 0.0, 0.0]
        @test Lambert.get_transfer_angle(r1_col, r2_opp, true) ≈ π atol = 1e-10
    end

    @testset "get_orbit_normal_vector" begin
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 7000.0, 0.0]
        i_h_pro = Lambert.get_orbit_normal_vector(r1, r2, true)
        @test i_h_pro[3] > 0
        @test norm(i_h_pro) ≈ 1.0 atol = 1e-10

        i_h_retro = Lambert.get_orbit_normal_vector(r1, r2, false)
        @test i_h_retro[3] < 0
        @test norm(i_h_retro) ≈ 1.0 atol = 1e-10
    end

    @testset "lambert_geometry" begin
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 7000.0, 0.0]
        r1_n, r2_n, c_n, dθ = Lambert.lambert_geometry(r1, r2, true)
        @test r1_n ≈ 7000.0
        @test r2_n ≈ 7000.0
        @test c_n ≈ sqrt(7000^2 + 7000^2) atol = 1e-6
        @test dθ ≈ π / 2 atol = 1e-10
    end

    @testset "compute_unit_vector_geometry collinear fallback" begin
        r1 = [7000.0, 0.0, 0.0]
        r2 = [14000.0, 0.0, 0.0]
        i_r1, i_r2, i_h, cos_dθ, sin_dθ =
            Lambert.compute_unit_vector_geometry(r1, r2, 7000.0, 14000.0)
        @test cos_dθ ≈ 1.0 atol = 1e-10
        @test sin_dθ ≈ 0.0 atol = 1e-10
        @test i_h ≈ SVector{3}(0.0, 0.0, 1.0)
    end

    @testset "compute_unit_vector_geometry normal case" begin
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 7000.0, 0.0]
        i_r1, i_r2, i_h, cos_dθ, sin_dθ =
            Lambert.compute_unit_vector_geometry(r1, r2, 7000.0, 7000.0)
        @test cos_dθ ≈ 0.0 atol = 1e-10
        @test sin_dθ ≈ 1.0 atol = 1e-10
        @test norm(i_h) ≈ 1.0 atol = 1e-10
    end

    # ── Convergence & status helpers ──────────────────────────────────────
    @testset "check_collinear_vectors" begin
        @test Lambert.check_collinear_vectors(0.0) == :COLLINEAR_VECTORS
        @test Lambert.check_collinear_vectors(1e-12) == :COLLINEAR_VECTORS
        @test Lambert.check_collinear_vectors(1.0) == :OK
    end

    @testset "check_convergence scalar" begin
        @test Lambert.check_convergence(1.0, 1.0, 1e-10, 1e-10) == true
        @test Lambert.check_convergence(1.0, 2.0, 1e-10, 1e-10) == false
        @test Lambert.check_convergence(1.0, 1.0 + 1e-12, 1e-10, 1e-10) == true
    end

    @testset "check_convergence vector" begin
        v1 = [1.0, 2.0, 3.0]
        v2 = [1.0, 2.0, 3.0]
        @test Lambert.check_convergence(v1, v2, 1e-10, 1e-10) == true
        @test Lambert.check_convergence(v1, [1.0, 2.0, 4.0], 1e-10, 1e-10) == false
        @test Lambert.check_convergence(v1, [1.0 + 1e-12, 2.0, 3.0], 1e-10, 1e-10) == true
    end

    @testset "handle_max_iterations" begin
        @test Lambert.handle_max_iterations(35, 35) == :MAXIMUM_ITERATIONS
        @test Lambert.handle_max_iterations(36, 35) == :MAXIMUM_ITERATIONS
        @test Lambert.handle_max_iterations(34, 35) == :SUCCESS
        @test Lambert.handle_max_iterations(0, 35) == :SUCCESS
    end

    @testset "check_newton_derivative" begin
        @test Lambert.check_newton_derivative(1.0) == :OK
        @test Lambert.check_newton_derivative(0.0) == :DERIVATIVE_TOO_SMALL
        @test Lambert.check_newton_derivative(1e-15) == :DERIVATIVE_TOO_SMALL
    end

    # ── Numerical helpers ─────────────────────────────────────────────────
    @testset "robust_sqrt" begin
        @test Lambert.robust_sqrt(4.0) ≈ 2.0
        @test Lambert.robust_sqrt(0.0) ≈ 0.0
        @test Lambert.robust_sqrt(-1e-20) ≈ 0.0
        @test Lambert.robust_sqrt(1e-20) ≈ 0.0
    end

    @testset "robust_division" begin
        @test Lambert.robust_division(10.0, 2.0) ≈ 5.0
        @test Lambert.robust_division(10.0, 0.0) ≈ 0.0
        @test Lambert.robust_division(10.0, 1e-20) ≈ 0.0
        @test Lambert.robust_division(10.0, 0.0; fallback = 999.0) ≈ 999.0
    end

    @testset "apply_damping" begin
        @test Lambert.apply_damping(0.1, 1.0) ≈ 0.1
        delta_damped = Lambert.apply_damping(10.0, 1.0)
        @test abs(delta_damped) < abs(10.0)
        @test abs(delta_damped / 1.0) ≈ 0.5
    end

    @testset "validate_lagrange_coefficients" begin
        @test Lambert.validate_lagrange_coefficients(0.5, 100.0, 0.8) == :OK
        @test Lambert.validate_lagrange_coefficients(Inf, 100.0, 0.8) ==
              :INVALID_LAGRANGE_COEFFICIENTS
        @test Lambert.validate_lagrange_coefficients(0.5, NaN, 0.8) ==
              :INVALID_LAGRANGE_COEFFICIENTS
        @test Lambert.validate_lagrange_coefficients(0.5, 1e-20, 0.8) ==
              :INVALID_LAGRANGE_COEFFICIENTS
    end

    # ── Velocity reconstruction ───────────────────────────────────────────
    @testset "create_lambert_solution" begin
        v1 = SVector{3}(1.0, 2.0, 3.0)
        v2 = SVector{3}(4.0, 5.0, 6.0)
        sol = Lambert.create_lambert_solution(v1, v2, 10, :SUCCESS)
        @test sol isa LambertSolution
        @test sol.v1 == v1
        @test sol.v2 == v2
        @test sol.numiter == 10
        @test sol.retcode == :SUCCESS
    end

    @testset "compute_semiperimeter" begin
        @test Lambert.compute_semiperimeter(7000.0, 8000.0, 5000.0) ≈
              (7000.0 + 8000.0 + 5000.0) / 2.0
    end

    @testset "compute_unit_vectors" begin
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 7000.0, 0.0]
        i_r1, i_r2, i_h_unnorm = Lambert.compute_unit_vectors(r1, r2, 7000.0, 7000.0)
        @test i_r1 ≈ SVector{3}(1.0, 0.0, 0.0)
        @test i_r2 ≈ SVector{3}(0.0, 1.0, 0.0)
        @test norm(i_h_unnorm) > 0
    end

    @testset "compute_tangential_unit_vectors" begin
        i_h = SVector{3}(0.0, 0.0, 1.0)
        i_r1 = SVector{3}(1.0, 0.0, 0.0)
        i_r2 = SVector{3}(0.0, 1.0, 0.0)
        i_t1, i_t2 = Lambert.compute_tangential_unit_vectors(i_h, i_r1, i_r2)
        @test norm(i_t1) ≈ 1.0 atol = 1e-10
        @test norm(i_t2) ≈ 1.0 atol = 1e-10
        @test dot(i_t1, i_r1) ≈ 0.0 atol = 1e-10
        @test dot(i_t2, i_r2) ≈ 0.0 atol = 1e-10
    end

    @testset "reconstruct_velocities_from_components" begin
        i_r1 = SVector{3}(1.0, 0.0, 0.0)
        i_t1 = SVector{3}(0.0, 1.0, 0.0)
        i_r2 = SVector{3}(0.0, 1.0, 0.0)
        i_t2 = SVector{3}(-1.0, 0.0, 0.0)
        v1, v2 = Lambert.reconstruct_velocities_from_components(
            1.0,
            2.0,
            3.0,
            4.0,
            i_r1,
            i_t1,
            i_r2,
            i_t2,
        )
        @test v1 ≈ SVector{3}(1.0, 2.0, 0.0)
        @test v2 ≈ SVector{3}(-4.0, 3.0, 0.0)
    end

    @testset "reconstruct_velocities_fg" begin
        r1 = SVector{3}(7000.0, 0.0, 0.0)
        r2 = SVector{3}(0.0, 7000.0, 0.0)
        f = 0.5
        g = 100.0
        gdot = 0.8
        v1, v2 = Lambert.reconstruct_velocities_fg(f, g, gdot, r1, r2)
        @test v1 ≈ (r2 - f * r1) / g
        @test v2 ≈ (gdot * r2 - r1) / g
    end

    # ── Normalization ─────────────────────────────────────────────────────
    @testset "normalize_inputs and denormalize_velocities" begin
        μ = 398600.0
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 8000.0, 0.0]
        tof = 3600.0

        μ_hat, r1_hat, r2_hat, tof_hat, L_ref, T_ref =
            Lambert.normalize_inputs(μ, r1, r2, tof)
        @test μ_hat ≈ 1.0 atol = 1e-10
        @test norm(r1_hat) ≈ 1.0 atol = 1e-10
        @test L_ref ≈ norm(r1)

        v1_hat = [0.001, 0.002, 0.0]
        v2_hat = [-0.001, 0.001, 0.0]
        v1_dim, v2_dim = Lambert.denormalize_velocities(v1_hat, v2_hat, L_ref, T_ref)
        scale = L_ref / T_ref
        @test v1_dim ≈ v1_hat * scale
        @test v2_dim ≈ v2_hat * scale
    end

    # ── Miscellaneous ─────────────────────────────────────────────────────
    @testset "compute_parabolic_time_of_flight" begin
        tof_p = Lambert.compute_parabolic_time_of_flight(398600.0, 10000.0, 5000.0)
        @test tof_p > 0
        @test tof_p < 1e6
    end

    @testset "universal_anomaly_initial_guess" begin
        μ = 398600.0
        s = 10000.0
        χ_ell = Lambert.universal_anomaly_initial_guess(μ, 5000.0, s, 0, 3000.0)
        @test χ_ell > 0
        χ_hyp = Lambert.universal_anomaly_initial_guess(μ, 2000.0, s, 0, 3000.0)
        @test χ_hyp < 0
        χ_mr = Lambert.universal_anomaly_initial_guess(μ, 5000.0, s, 2, 3000.0)
        @test χ_mr > 0
    end

    @testset "extract_common_solver_params" begin
        solver = IzzoSolver(M = 2, prograde = false, maxiter = 50, atol = 1e-6, rtol = 1e-8)
        params = Lambert.extract_common_solver_params(solver)
        @test params.M == 2
        @test params.prograde == false
        @test params.maxiter == 50
        @test params.atol == 1e-6
        @test params.rtol == 1e-8

        solver_min = AvanziniSolver()
        params_min = Lambert.extract_common_solver_params(solver_min)
        @test params_min.M == 0
        @test params_min.prograde == true
    end

    @testset "transfer_classification" begin
        μ = 398600.4418
        r1 = [7000.0, 0.0, 0.0]

        r2_short = [6062.2, 3500.0, 0.0]
        @test Lambert.transfer_classification(r1, r2_short, true, 1800.0, μ) isa Symbol

        r2_pi = [-6900.0, 500.0, 0.0]
        @test Lambert.transfer_classification(r1, r2_pi, true, 5400.0, μ) ==
              :near_pi_transfer

        r2_med = [-3500.0, 6062.2, 0.0]
        @test Lambert.transfer_classification(r1, r2_med, true, 3600.0, μ) isa Symbol

        r2_large = [6062.2, -3500.0, 0.0]
        @test Lambert.transfer_classification(r1, r2_large, false, 7200.0, μ) isa Symbol
    end

    @testset "normalized_time_of_flight_function" begin
        params = (a = 1.0,)
        func(x, p) = x^2 + p.a
        residual = Lambert.normalized_time_of_flight_function(2.0, 3.0, params, func)
        @test residual ≈ (4.0 + 1.0) - 3.0
    end
end
