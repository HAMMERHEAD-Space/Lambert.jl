@testset "SciMLBase Interface" begin
    # ── LambertIterator init / solve! ─────────────────────────────────────
    @testset "LambertIterator init/solve!" begin
        μ = 3.986004418e5
        r1 = [15945.34, 0.0, 0.0]
        r2 = [12214.83899, 10249.46731, 0.0]
        tof = 76.0 * 60
        prob = LambertProblem(μ, r1, r2, tof)

        iter = SciMLBase.init(prob, IzzoSolver())
        @test iter isa Lambert.LambertIterator
        @test iter.prob === prob
        @test iter.alg isa IzzoSolver

        sol = SciMLBase.solve!(iter)
        @test sol isa LambertSolution
        @test sol.retcode == :SUCCESS
    end

    @testset "init with heuristic (no alg)" begin
        μ = 3.986004418e5
        r1 = [15945.34, 0.0, 0.0]
        r2 = [12214.83899, 10249.46731, 0.0]
        tof = 76.0 * 60
        prob = LambertProblem(μ, r1, r2, tof)

        iter = SciMLBase.init(prob)
        @test iter isa Lambert.LambertIterator

        sol = SciMLBase.solve!(iter)
        @test sol isa LambertSolution
        @test sol.retcode == :SUCCESS
    end

    @testset "init with M and prograde kwargs" begin
        μ = 3.986004418e5
        r1 = [15945.34, 0.0, 0.0]
        r2 = [12214.83899, 10249.46731, 0.0]
        tof = 76.0 * 60 * 5
        prob = LambertProblem(μ, r1, r2, tof)

        iter = SciMLBase.init(prob; M = 1, prograde = true)
        @test iter isa Lambert.LambertIterator

        sol = SciMLBase.solve!(iter)
        @test sol isa LambertSolution
        @test sol.retcode == :SUCCESS
    end

    # ── _apply_kwargs_to_solver ───────────────────────────────────────────
    @testset "_apply_kwargs_to_solver" begin
        solver = IzzoSolver()
        new_solver = Lambert._apply_kwargs_to_solver(solver, (M = 2, prograde = false))
        @test new_solver.M == 2
        @test new_solver.prograde == false

        same_solver = Lambert._apply_kwargs_to_solver(solver, NamedTuple())
        @test same_solver.M == solver.M
        @test same_solver.prograde == solver.prograde

        extra_solver = Lambert._apply_kwargs_to_solver(solver, (nonexistent_field = 42,))
        @test extra_solver.M == solver.M
    end

    # ── solve with kwargs passthrough ─────────────────────────────────────
    @testset "solve with kwargs passthrough" begin
        μ = 3.986004418e5
        r1 = [15945.34, 0.0, 0.0]
        r2 = [12214.83899, 10249.46731, 0.0]
        tof = 76.0 * 60
        prob = LambertProblem(μ, r1, r2, tof)

        sol = solve(prob, IzzoSolver(); maxiter = 50)
        @test sol.retcode == :SUCCESS

        sol_retro = solve(prob, GoodingSolver(); prograde = false)
        @test sol_retro isa LambertSolution
    end

    # ── Heuristic algorithm selection ─────────────────────────────────────
    @testset "select_lambert_algorithm ultra-high M" begin
        μ = 398600.4418
        r1 = [7000.0, 0.0, 0.0]
        r2 = [-3500.0, 6062.2, 0.0]
        prob = LambertProblem(μ, r1, r2, 3600.0)

        alg = select_lambert_algorithm(prob, 15)
        @test alg isa McElreathSolver
        @test alg.M == 15

        alg100 = select_lambert_algorithm(prob, 100)
        @test alg100 isa McElreathSolver
        @test alg100.M == 100
    end

    @testset "select_lambert_algorithm retrograde" begin
        μ = 398600.4418
        r1 = [7000.0, 0.0, 0.0]
        r2 = [-3500.0, 6062.2, 0.0]
        prob = LambertProblem(μ, r1, r2, 3600.0)

        alg = select_lambert_algorithm(prob, 0, false)
        @test alg isa Lambert.AbstractLambertSolver
        @test alg.prograde == false
    end

    # ── Auto-solve (no solver specified) ──────────────────────────────────
    @testset "Solve with automatic algorithm" begin
        μ = 3.986004418e5
        r1 = [15945.34, 0.0, 0.0]
        r2 = [12214.83899, 10249.46731, 0.0]
        tof = 76.0 * 60
        prob = LambertProblem(μ, r1, r2, tof)

        sol = solve(prob)
        @test sol isa LambertSolution
        @test sol.retcode == :SUCCESS
    end

    @testset "Solve with M kwarg auto-select" begin
        μ = 3.986004418e5
        r1 = [15945.34, 0.0, 0.0]
        r2 = [12214.83899, 10249.46731, 0.0]
        tof = 76.0 * 60 * 5
        prob = LambertProblem(μ, r1, r2, tof)

        sol = solve(prob; M = 1)
        @test sol isa LambertSolution
        @test sol.retcode == :SUCCESS
    end

    @testset "Solve with prograde=false kwarg auto-select" begin
        μ = 3.986004418e5
        r1 = [7100.0, 200.0, 1300.0]
        r2 = [-47332.7499, -54840.2027, -37100.17067]
        tof = 12000.0
        prob = LambertProblem(μ, r1, r2, tof)

        sol = solve(prob; prograde = false)
        @test sol isa LambertSolution
        @test sol.retcode == :SUCCESS
    end

    # ── remake ────────────────────────────────────────────────────────────
    @testset "remake with all fields" begin
        μ = 3.986004418e5
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 7000.0, 0.0]
        tof = 3600.0
        prob = LambertProblem(μ, r1, r2, tof)

        new_prob =
            remake(prob; μ = 1.0, r1 = [1.0, 0.0, 0.0], r2 = [0.0, 1.0, 0.0], tof = 1.0)
        @test new_prob.μ == 1.0
        @test new_prob.r1 == [1.0, 0.0, 0.0]
        @test new_prob.r2 == [0.0, 1.0, 0.0]
        @test new_prob.tof == 1.0
    end
end
