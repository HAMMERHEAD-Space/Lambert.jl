@testset "Solver Internals" begin
    # Reference problem for cross-solver velocity comparisons
    _REF_μ = 3.986004418e5
    _REF_r1 = [15945.34, 0.0, 0.0]
    _REF_r2 = [12214.83899, 10249.46731, 0.0]
    _REF_tof = 76.0 * 60

    ref_sol = solve(LambertProblem(_REF_μ, _REF_r1, _REF_r2, _REF_tof), IzzoSolver())

    # ═══════════════════════════════════════════════════════════════════════
    # Izzo
    # ═══════════════════════════════════════════════════════════════════════
    @testset "Izzo hyp2f1b" begin
        @test Lambert.hyp2f1b(1.0) == Inf
        @test Lambert.hyp2f1b(2.0) == Inf
        @test Lambert.hyp2f1b(0.0) ≈ 1.0
        @test Lambert.hyp2f1b(0.1) ≈ 1.1354 atol = 1e-3
        @test Lambert.hyp2f1b(0.5) > Lambert.hyp2f1b(0.1)
    end

    @testset "Izzo compute_psi branches" begin
        # Elliptic: ψ = acos(x*y + ll*(1-x²)), should be in (0, π)
        ψ_ell = Lambert.compute_psi(0.5, Lambert.compute_y(0.5, 0.3), 0.3)
        @test 0 < ψ_ell < π

        # Hyperbolic: ψ = asinh(...), should be positive
        ψ_hyp = Lambert.compute_psi(1.5, Lambert.compute_y(1.5, 0.3), 0.3)
        @test ψ_hyp > 0
    end

    @testset "Izzo compute_y" begin
        @test Lambert.compute_y(0.0, 0.5) ≈ sqrt(1 - 0.5^2)
        @test Lambert.compute_y(1.0, 0.0) ≈ 1.0
        # y decreases as |ll| increases for fixed x
        @test Lambert.compute_y(0.5, 0.8) < Lambert.compute_y(0.5, 0.3)
    end

    @testset "Izzo compute_T_min" begin
        # ll == 1: x_T_min = 0, T_min = π for M=1
        x_Tmin, T_min = Lambert.compute_T_min(1.0, 1, 35, 1e-5, 1e-7)
        @test x_Tmin ≈ 0.0
        @test T_min ≈ π atol = 1e-6

        # M == 0: no minimum needed
        x_Tmin0, T_min0 = Lambert.compute_T_min(0.5, 0, 35, 1e-5, 1e-7)
        @test x_Tmin0 == Inf
        @test T_min0 == 0.0

        # Multi-rev: should find a valid minimum with T_min > 0
        x_Tmin_mr, T_min_mr = Lambert.compute_T_min(0.5, 1, 35, 1e-5, 1e-7)
        @test -1 < x_Tmin_mr < 1
        @test T_min_mr > 0
    end

    @testset "Izzo initial_guess branches" begin
        T_0 = acos(0.3) + 0.3 * sqrt(1 - 0.3^2)
        T_1 = 2 * (1 - 0.3^3) / 3

        # T >= T_0: x_0 = (T_0/T)^(2/3) - 1, should be negative
        x_0 = Lambert.initial_guess(5.0, 0.3, 0, true)
        @test x_0 ≈ (T_0 / 5.0)^(2 / 3) - 1 atol = 1e-10

        # T < T_1: x_0 > 1 (hyperbolic initial guess)
        x_0_small = Lambert.initial_guess(0.01, 0.3, 0, true)
        @test x_0_small > 1

        # T_1 <= T < T_0: intermediate branch
        T_mid = (T_0 + T_1) / 2
        x_0_mid = Lambert.initial_guess(T_mid, 0.3, 0, true)
        @test -1 < x_0_mid < 1

        # M > 0 low_path: max(x_0l, x_0r)
        x_0_mr_low = Lambert.initial_guess(5.0, 0.3, 1, true)
        x_0_mr_high = Lambert.initial_guess(5.0, 0.3, 1, false)
        @test x_0_mr_low >= x_0_mr_high
    end

    @testset "Izzo tof_equation consistency" begin
        ll = 0.3
        T0_ref = 0.5

        # M=0 hyp2f1b path (x in sqrt(0.6)..sqrt(1.4))
        x_mid = 0.9
        y_mid = Lambert.compute_y(x_mid, ll)
        res = Lambert.tof_equation_y(x_mid, y_mid, T0_ref, ll, 0)
        @test res isa Float64

        # M > 0 compute_psi path — different from M=0 path
        x_out = 0.3
        y_out = Lambert.compute_y(x_out, ll)
        res0 = Lambert.tof_equation_y(x_out, y_out, T0_ref, ll, 0)
        res1 = Lambert.tof_equation_y(x_out, y_out, T0_ref, ll, 1)
        @test res1 > res0
    end

    @testset "Izzo multi-rev infeasible M" begin
        @test_throws ErrorException Lambert.izzo2015(
            _REF_μ,
            _REF_r1,
            _REF_r2,
            _REF_tof;
            M = 100,
        )
    end

    # ═══════════════════════════════════════════════════════════════════════
    # Gauss
    # ═══════════════════════════════════════════════════════════════════════
    @testset "Gauss solver agrees with reference" begin
        v1, v2, numiter, status = Lambert.gauss1809(_REF_μ, _REF_r1, _REF_r2, _REF_tof)
        @test status == :SUCCESS
        @test v1 ≈ ref_sol.v1 atol = 0.5
        @test v2 ≈ ref_sol.v2 atol = 0.5
    end

    @testset "Gauss X_at_x" begin
        # X_at_x(0) = (4/3) * 1.0 since all series terms vanish
        @test Lambert.X_at_x(0.0) ≈ 4 / 3
        @test Lambert.X_at_x(0.1) > 4 / 3
    end

    # ═══════════════════════════════════════════════════════════════════════
    # Avanzini
    # ═══════════════════════════════════════════════════════════════════════
    @testset "Avanzini M > 0 error" begin
        @test_throws ErrorException Lambert.avanzini2008(
            _REF_μ,
            _REF_r1,
            _REF_r2,
            _REF_tof;
            M = 1,
        )
    end

    @testset "Avanzini custom parameters" begin
        prob = LambertProblem(_REF_μ, _REF_r1, _REF_r2, _REF_tof)
        sol = solve(prob, AvanziniSolver(maxiter = 50, atol = 1e-6))
        @test sol.retcode == :SUCCESS
        @test sol.v1 ≈ ref_sol.v1 atol = 0.1
    end

    # ═══════════════════════════════════════════════════════════════════════
    # Battin
    # ═══════════════════════════════════════════════════════════════════════
    @testset "Battin helper functions" begin
        # _get_λ: positive when dtheta < π, negative when dtheta > π
        λ_less = Lambert._get_λ(5000.0, 10000.0, π / 3)
        @test λ_less > 0
        λ_more = Lambert._get_λ(5000.0, 10000.0, 4.0)
        @test λ_more < 0

        # _get_ll: known formula
        @test Lambert._get_ll(0.5) ≈ ((1 - 0.5) / (1 + 0.5))^2

        # _get_m: m = 8μt²/(s³(1+λ)⁶)
        m = Lambert._get_m(398600.0, 3600.0, 10000.0, 0.5)
        @test m ≈ 8 * 398600.0 * 3600.0^2 / (10000.0^3 * 1.5^6) rtol = 1e-10
        @test m > 0

        # _xi_at_x: at x=0, η=0, σ→1, ξ ≈ 8*2/(3+1/(5+0)) = 16/3.2 = 5
        @test Lambert._xi_at_x(0.0) ≈ 5.0 atol = 1e-6
        @test Lambert._xi_at_x(0.1) > 0

        # _K_at_u: continued fraction, should be positive for small positive u
        K = Lambert._K_at_u(0.1)
        @test K > 0
        @test K < 10

        # _B_at_h: B = 27h₂/(4(1+h₁)³)
        B = Lambert._B_at_h(0.1, 0.2)
        @test B ≈ 27 * 0.2 / (4 * 1.1^3) rtol = 1e-10

        # _u_at_B roundtrips through B
        u = Lambert._u_at_B(B)
        @test u < 0
        @test Lambert._u_at_h(0.1, 0.2) ≈ u
    end

    @testset "Battin dtheta > π" begin
        prob = LambertProblem(
            3.986004418e5,
            [7000.0, 0.0, 0.0],
            [6062.2, -3500.0, 0.0],
            5400.0,
        )
        sol = solve(prob, BattinSolver(prograde = false))
        @test sol isa LambertSolution
    end

    # ═══════════════════════════════════════════════════════════════════════
    # Arora
    # ═══════════════════════════════════════════════════════════════════════
    @testset "Arora get_W branches" begin
        # W(k=0, M=0) = π/(2√2) - 0/2 ≈ 1.1107
        # From the elliptic formula: W = (π+0)/(2^(3/2)) - 0/2
        W_0 = Lambert.get_W(0.5, 0)
        @test W_0 > 0

        # Parabolic band: series should match smoothly at boundaries
        W_para_lo = Lambert.get_W(√2 - 0.01, 0)
        W_para_hi = Lambert.get_W(√2 + 0.01, 0)
        @test W_para_lo > W_para_hi

        # Hyperbolic: W > 0 and decreasing
        W_hyp = Lambert.get_W(2.0, 0)
        @test W_hyp > 0
        @test W_hyp < W_0

        # Negative k: larger W (longer TOF)
        W_neg = Lambert.get_W(-1.0, 0)
        @test W_neg > W_0

        # Multi-rev adds 2πM to the numerator — W should be larger
        W_mr = Lambert.get_W(0.5, 1)
        @test W_mr > W_0
    end

    @testset "Arora get_Wprime sign" begin
        # W is decreasing in k for k < √2 → W' < 0
        for k in [0.5, √2 - 0.005]
            W = Lambert.get_W(k, 0)
            Wp = Lambert.get_Wprime(k, W)
            @test Wp < 0
        end
        # Hyperbolic
        W_hyp = Lambert.get_W(2.0, 0)
        Wp_hyp = Lambert.get_Wprime(2.0, W_hyp)
        @test Wp_hyp < 0
    end

    @testset "Arora get_W2prime" begin
        for k in [0.5, √2 - 0.005, 2.0]
            W = Lambert.get_W(k, 0)
            Wp = Lambert.get_Wprime(k, W)
            Wpp = Lambert.get_W2prime(k, W, Wp)
            # Second derivative should be positive (W is convex)
            @test Wpp > 0
        end
    end

    @testset "Arora get_gammas" begin
        γ1, γ2, γ3 = Lambert.get_gammas(0.5, 0.3, 0.7)
        # γ₁ = F_i(F_star - F_n) = 0.5*(0.7-0.3) = 0.2
        @test γ1 ≈ 0.5 * (0.7 - 0.3)
        # γ₂ = F_star(F_n - F_i) = 0.7*(0.3-0.5) = -0.14
        @test γ2 ≈ 0.7 * (0.3 - 0.5)
        # γ₃ = F_n(F_star - F_i) = 0.3*(0.7-0.5) = 0.06
        @test γ3 ≈ 0.3 * (0.7 - 0.5)
    end

    @testset "Arora get_x" begin
        x = Lambert.get_x(1.0, 0.5, 0.75, 0.8, 0.5, 1.0)
        @test 0 < x < 1
    end

    @testset "Arora get_TOF" begin
        # TOF = S*√(1-kτ)*(τ + (1-kτ)*W)
        W = Lambert.get_W(0.5, 0)
        T = Lambert.get_TOF(0.5, 0.3, 10.0, W)
        expected = 10.0 * √(1 - 0.5 * 0.3) * (0.3 + (1 - 0.5 * 0.3) * W)
        @test T ≈ expected
        @test T > 0
    end

    @testset "Arora multi-rev" begin
        prob = LambertProblem(_REF_μ, _REF_r1, _REF_r2, _REF_tof * 7)
        sol = solve(prob, AroraSolver(M = 1))
        @test sol isa LambertSolution
    end

    @testset "Arora elliptic TOF regions" begin
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 8000.0, 0.0]

        for tof in [2000.0, 4000.0, 20000.0, 100000.0]
            sol = solve(LambertProblem(_REF_μ, r1, r2, tof), AroraSolver())
            @test sol isa LambertSolution
            @test sol.retcode == :SUCCESS
        end
    end

    @testset "Arora hyperbolic d=-1 regions" begin
        sol = solve(
            LambertProblem(_REF_μ, [7000.0, 0.0, 0.0], [6062.2, -3500.0, 0.0], 2000.0),
            AroraSolver(prograde = false),
        )
        @test sol isa LambertSolution
    end

    # ═══════════════════════════════════════════════════════════════════════
    # Russell
    # ═══════════════════════════════════════════════════════════════════════
    @testset "Russell _russell_W regions" begin
        # Near-zero Taylor series
        W0_zero, dW1_0, _, _ = Lambert._russell_W(0.001, 0, 3)
        # At k=0, W ≈ π/(2√2) ≈ 1.1107 for M=0
        @test W0_zero ≈ 1.1107 atol = 0.01

        # Normal ellipse k > 0
        W_pos, _, _, _ = Lambert._russell_W(0.5, 0, 3)
        @test 0 < W_pos < W0_zero

        # Normal ellipse k < 0
        W_neg, _, _, _ = Lambert._russell_W(-0.5, 0, 3)
        @test W_neg > W0_zero

        # Parabola band
        W_para, _, _, _ = Lambert._russell_W(√2 - 0.01, 0, 3)
        @test W_para > 0
        @test W_para < W_pos

        # Hyperbola
        W_hyp, _, _, _ = Lambert._russell_W(2.0, 0, 3)
        @test W_hyp > 0
        @test W_hyp < W_para

        # Multi-rev
        W_mr, _, _, _ = Lambert._russell_W(0.5, 1, 3)
        @test W_mr > W_pos
    end

    @testset "Russell _russell_W_parabola_series" begin
        # At ν=0 (k=√2), W should equal √2/3 ≈ 0.4714
        W0, dW1_0, dW2_0, dW3_0 = Lambert._russell_W_parabola_series(0.0, 3)
        @test W0 ≈ 0.47140452079103168 atol = 1e-12
        @test dW1_0 ≈ -0.2 atol = 1e-12
        @test dW2_0 ≈ 0.16162440712835372 atol = 1e-10

        # Continuity: W(ν) should be close at small ν
        Wp, _, _, _ = Lambert._russell_W_parabola_series(0.01, 3)
        Wm, _, _, _ = Lambert._russell_W_parabola_series(-0.01, 3)
        @test abs(Wp - W0) < 0.01
        @test abs(Wm - W0) < 0.01

        W1, dW1_1, _, _ = Lambert._russell_W_parabola_series(0.0, 1)
        @test W1 ≈ W0
    end

    @testset "Russell _russell_W order variations" begin
        # All orders should give the same W value — only derivatives differ
        W_o1, _, _, _ = Lambert._russell_W(0.5, 0, 1)
        W_o2, _, _, _ = Lambert._russell_W(0.5, 0, 2)
        W_o3, _, _, _ = Lambert._russell_W(0.5, 0, 3)
        @test W_o1[1] ≈ W_o2[1]
        @test W_o2[1] ≈ W_o3[1]
    end

    @testset "Russell solver order sweep" begin
        prob = LambertProblem(_REF_μ, _REF_r1, _REF_r2, _REF_tof)
        for order in [1, 2, 3]
            sol = solve(prob, RussellSolver(order = order))
            @test sol.retcode == :SUCCESS
            @test sol.v1 ≈ ref_sol.v1 atol = 1e-4
        end
    end

    @testset "Russell multi-rev low/high path" begin
        prob = LambertProblem(_REF_μ, _REF_r1, _REF_r2, _REF_tof * 5)
        sol_low = solve(prob, RussellSolver(M = 1, low_path = true))
        sol_high = solve(prob, RussellSolver(M = 1, low_path = false))
        @test sol_low isa LambertSolution
        @test sol_high isa LambertSolution
        # Low and high paths should give different velocities
        @test !(sol_low.v1 ≈ sol_high.v1)
    end

    @testset "Russell _russell_initial_guess_multirev no-solution" begin
        # Very short tof → no multi-rev solution exists → NaN
        k0, _ = Lambert._russell_initial_guess_multirev(0.5, 1.0, 0.001, 5, true)
        @test isnan(k0)
    end

    @testset "Russell _russell_correction trust region" begin
        info, dv1, _, _ = Lambert._russell_correction(1.0, 0.01, 0.0, 0.0, 1, false)
        @test info == 3 || dv1 ≈ Lambert._RUSSELL_MAX_K_STEP_MREV

        info2, dv1_2, _, _ = Lambert._russell_correction(1.0, -0.01, 0.0, 0.0, 1, false)
        @test info2 == 3 || dv1_2 ≈ -Lambert._RUSSELL_MAX_K_STEP_MREV
    end

    @testset "Russell near-180 transfer" begin
        prob = LambertProblem(_REF_μ, [7000.0, 0.0, 0.0], [-6900.0, 500.0, 0.0], 10000.0)
        sol = solve(prob, RussellSolver())
        @test sol.retcode == :SUCCESS
    end

    @testset "Russell _russell_W near -√2" begin
        # Near -√2 the kps2 series expansion should give positive W
        W1, _, _, _ = Lambert._russell_W(-√2 + 0.005, 1, 3)
        W2, _, _, _ = Lambert._russell_W(-√2 + 0.001, 1, 3)
        @test W1[1] > 0
        @test W2[1] > W1[1]
    end

    @testset "Russell hyperbolic" begin
        prob = LambertProblem(_REF_μ, [7000.0, 0.0, 0.0], [50000.0, 50000.0, 0.0], 8000.0)
        sol_r = solve(prob, RussellSolver())
        sol_i = solve(prob, IzzoSolver())
        @test sol_r.retcode == :SUCCESS
        @test sol_r.v1 ≈ sol_i.v1 atol = 1e-4
    end

    # ═══════════════════════════════════════════════════════════════════════
    # McElreath
    # ═══════════════════════════════════════════════════════════════════════
    @testset "McElreath universal functions" begin
        # Elliptic: verify trig identities
        U0, U1, U2, U3, U0s, U1s = Lambert._mcelreath_universal_elliptic(1.0, 0.5)
        @test U0 ≈ cos(1.0)
        @test U1 ≈ sin(1.0) / √0.5
        @test U2 ≈ (1 - cos(1.0)) / 0.5
        # U0² + α*U1² = 1 (Pythagorean identity)
        @test U0^2 + 0.5 * U1^2 ≈ 1.0 atol = 1e-10

        # Hyperbolic: verify hyperbolic identities
        U0h, U1h, U2h, U3h, U0sh, U1sh = Lambert._mcelreath_universal_hyperbolic(1.0, -0.5)
        @test U0h ≈ cosh(1.0)
        @test U1h ≈ sinh(1.0) / √0.5
        @test U2h ≈ (cosh(1.0) - 1) / 0.5
        # U0² - |α|*U1² = 1 (hyperbolic identity)
        @test U0h^2 - 0.5 * U1h^2 ≈ 1.0 atol = 1e-10
    end

    @testset "McElreath alpha functions" begin
        α_ell = Lambert._mcelreath_alpha_elliptic(1.0, 7000.0, 8000.0, 5000.0)
        @test α_ell > 0

        α_hyp = Lambert._mcelreath_alpha_hyperbolic(1.0, 7000.0, 8000.0, 5000.0)
        @test α_hyp < 0
    end

    @testset "McElreath multi-rev high M" begin
        prob = LambertProblem(_REF_μ, _REF_r1, _REF_r2, _REF_tof * 50)
        sol = solve(prob, McElreathSolver(M = 10))
        @test sol isa LambertSolution
        @test sol.retcode == :SUCCESS
    end

    @testset "McElreath low_path and high_path" begin
        prob = LambertProblem(_REF_μ, _REF_r1, _REF_r2, _REF_tof * 5)
        sol_low = solve(prob, McElreathSolver(M = 1, low_path = true))
        sol_high = solve(prob, McElreathSolver(M = 1, low_path = false))
        @test sol_low isa LambertSolution
        @test sol_high isa LambertSolution
        @test !(sol_low.v1 ≈ sol_high.v1)
    end

    # ═══════════════════════════════════════════════════════════════════════
    # Vallado
    # ═══════════════════════════════════════════════════════════════════════
    @testset "Vallado multi-rev" begin
        prob = LambertProblem(_REF_μ, _REF_r1, _REF_r2, _REF_tof * 7)
        sol = solve(prob, ValladoSolver(M = 1))
        @test sol isa LambertSolution
    end

    @testset "Vallado helper functions" begin
        # get_A: positive for dtheta < π, negative for dtheta > π
        A1 = Lambert.get_A(7000.0, 8000.0, π / 2)
        @test A1 > 0
        @test A1 ≈ sin(π / 2) * √(7000.0 * 8000.0 / (1 - cos(π / 2)))

        A2 = Lambert.get_A(7000.0, 8000.0, 4.0)
        @test A2 < 0

        # X_at_psi(0, y, tol): C2(0)=0.5, X = √(y/C2) = √(1/0.5) = √2
        X = Lambert.X_at_psi(0.0, 1.0, 1e-6)
        @test X ≈ √2 atol = 1e-6

        # y_at_psi should be positive for physical inputs
        y = Lambert.y_at_psi(1.0, 7000.0, 8000.0, 5000.0, 1e-6)
        @test y > 0
    end

    # ═══════════════════════════════════════════════════════════════════════
    # Gooding
    # ═══════════════════════════════════════════════════════════════════════
    @testset "Gooding tlamb derivative mode (n=-1)" begin
        t, dt, d2t, d3t = Lambert.tlamb(0, 0.5, 0.75, 0.5, -1)
        # Derivative should be non-zero at this point
        @test dt != 0.0
        @test d2t != 0.0
    end

    @testset "Gooding tlamb n=0 through n=3" begin
        results = Float64[]
        for n in [0, 1, 2, 3]
            t, _, _, _ = Lambert.tlamb(0, 0.5, 0.75, 0.5, n)
            push!(results, t)
        end
        @test all(r -> r > 0, results)
        # Different evaluation modes should produce the same t
        @test all(≈(results[1]), results)
    end

    @testset "Gooding tlamb series branch" begin
        # x near 1 triggers the series expansion
        t, _, _, _ = Lambert.tlamb(0, 0.5, 0.75, 0.95, 3)
        @test t > 0
    end

    @testset "Gooding multi-rev low/high path" begin
        prob = LambertProblem(_REF_μ, _REF_r1, _REF_r2, _REF_tof * 5)
        sol_low = solve(prob, GoodingSolver(M = 1, low_path = true))
        sol_high = solve(prob, GoodingSolver(M = 1, low_path = false))
        @test sol_low isa LambertSolution
        @test sol_high isa LambertSolution
        @test !(sol_low.v1 ≈ sol_high.v1)
    end

    @testset "Gooding xlamb paths" begin
        # m=0: should find at least one solution
        n, x1, x2, numiter = Lambert.xlamb(0, 0.5, 0.75, 5.0, 35, 1e-5, 1e-7)
        @test n >= 1
        @test -1 < x1 < 2

        # m>0: may find two solutions
        n2, x1_mr, x2_mr, numiter2 = Lambert.xlamb(1, 0.5, 0.75, 15.0, 35, 1e-5, 1e-7)
        @test -1 < x1_mr < 2
    end

    # ═══════════════════════════════════════════════════════════════════════
    # Direct function-level calls (all solvers, cross-validated)
    # ═══════════════════════════════════════════════════════════════════════
    @testset "Direct solver function calls agree" begin
        v1_izzo, v2_izzo = Lambert.izzo2015(_REF_μ, _REF_r1, _REF_r2, _REF_tof)

        _, _, _, rc_g = Lambert.gooding1990(_REF_μ, _REF_r1, _REF_r2, _REF_tof)
        @test rc_g == :SUCCESS

        v1_v, v2_v, _, rcv = Lambert.vallado2013(_REF_μ, _REF_r1, _REF_r2, _REF_tof)
        @test rcv == :SUCCESS
        @test v1_v ≈ v1_izzo atol = 0.01

        v1_b, v2_b, _, rcb = Lambert.battin1984(_REF_μ, _REF_r1, _REF_r2, _REF_tof)
        @test rcb == :SUCCESS
        @test v1_b ≈ v1_izzo atol = 0.01

        v1_r, v2_r, _, rcr = Lambert.russell2021(_REF_μ, _REF_r1, _REF_r2, _REF_tof)
        @test rcr == :SUCCESS
        @test v1_r ≈ v1_izzo atol = 1e-4

        v1_m, v2_m, _, rcm = Lambert.mcelreath2025(_REF_μ, _REF_r1, _REF_r2, _REF_tof)
        @test rcm == :SUCCESS
        @test v1_m ≈ v1_izzo atol = 1e-4

        v1_a, v2_a, _, rca = Lambert.arora2013(_REF_μ, _REF_r1, _REF_r2, _REF_tof)
        @test rca == :SUCCESS
        @test v1_a ≈ v1_izzo atol = 0.01

        v1_ga, v2_ga, _, rcga = Lambert.gauss1809(_REF_μ, _REF_r1, _REF_r2, _REF_tof)
        @test rcga == :SUCCESS
        @test v1_ga ≈ v1_izzo atol = 0.5
    end

    # ═══════════════════════════════════════════════════════════════════════
    # Cross-solver edge cases
    # ═══════════════════════════════════════════════════════════════════════
    @testset "Retrograde transfer — all compatible solvers" begin
        prob = LambertProblem(
            _REF_μ,
            [7100.0, 200.0, 1300.0],
            [-47332.7499, -54840.2027, -37100.17067],
            12000.0,
        )

        solutions = LambertSolution[]
        for solver_type in [
            IzzoSolver,
            GoodingSolver,
            BattinSolver,
            ValladoSolver,
            RussellSolver,
            McElreathSolver,
        ]
            solver = solver_type(prograde = false)
            sol = solve(prob, solver)
            @test sol.retcode == :SUCCESS
            push!(solutions, sol)
        end

        # All solvers should agree on the departure velocity
        for sol in solutions[2:end]
            @test sol.v1 ≈ solutions[1].v1 atol = 0.1
        end
    end

    @testset "Hyperbolic transfer — multi-solver" begin
        prob = LambertProblem(
            _REF_μ,
            [7100.0, 200.0, 1300.0],
            [-38113.587, 67274.1946, 29309.5799],
            12000.0,
        )

        solutions = LambertSolution[]
        for solver in [IzzoSolver(), GoodingSolver(), RussellSolver(), McElreathSolver()]
            sol = solve(prob, solver)
            @test sol.retcode == :SUCCESS
            push!(solutions, sol)
        end

        for sol in solutions[2:end]
            @test sol.v1 ≈ solutions[1].v1 atol = 1e-3
        end
    end
end
