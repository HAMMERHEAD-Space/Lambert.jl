using AllocCheck

const μ_perf = 3.986004418e5
const r1_perf = [15945.34, 0.0, 0.0]
const r2_perf = [12214.83899, 10249.46731, 0.0]
const tof_perf = 76.0 * 60

@testset "Non-Allocating Performance Tests" begin

    @testset "Problem Construction — Zero Allocation" begin
        allocs_vec = check_allocs(
            LambertProblem,
            (typeof(μ_perf), typeof(r1_perf), typeof(r2_perf), typeof(tof_perf)),
        )
        @test length(allocs_vec) == 0
    end

    @testset "Heuristic Algorithm Selection — Zero Allocation" begin
        problem = LambertProblem(μ_perf, r1_perf, r2_perf, tof_perf)
        allocs_vec = check_allocs(select_lambert_algorithm, (typeof(problem),))
        @test length(allocs_vec) == 0
    end

    @testset "Utility Functions — Zero Allocation" begin
        @testset "Stumpff c2" begin
            allocs_vec = check_allocs(Lambert.c2, (Float64,))
            @test length(allocs_vec) == 0
        end

        @testset "Stumpff c3" begin
            allocs_vec = check_allocs(Lambert.c3, (Float64,))
            @test length(allocs_vec) == 0
        end

        @testset "Scalar check_convergence" begin
            allocs_vec = check_allocs(
                Lambert.check_convergence,
                (Float64, Float64, Float64, Float64),
            )
            @test length(allocs_vec) == 0
        end

        @testset "compute_semiperimeter" begin
            allocs_vec =
                check_allocs(Lambert.compute_semiperimeter, (Float64, Float64, Float64))
            @test length(allocs_vec) == 0
        end

        @testset "handle_max_iterations" begin
            allocs_vec = check_allocs(Lambert.handle_max_iterations, (Int, Int))
            @test length(allocs_vec) == 0
        end

        @testset "assert_transfer_angle_not_zero" begin
            allocs_vec = check_allocs(Lambert.assert_transfer_angle_not_zero, (Float64,))
            @test length(allocs_vec) == 0
        end

        @testset "assert_transfer_angle_not_pi" begin
            allocs_vec = check_allocs(Lambert.assert_transfer_angle_not_pi, (Float64,))
            @test length(allocs_vec) == 0
        end
    end

    @testset "Izzo Solver Utilities — Zero Allocation" begin
        allocs_vec = check_allocs(Lambert.compute_y, (Float64, Float64))
        @test length(allocs_vec) == 0

        allocs_vec = check_allocs(
            Lambert.reconstruct,
            (Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64),
        )
        @test length(allocs_vec) == 0
    end

    @testset "Vector Reconstruction — Zero Allocation" begin
        @testset "reconstruct_velocities_from_components" begin
            allocs_vec = check_allocs(
                Lambert.reconstruct_velocities_from_components,
                (
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    SVector{3,Float64},
                    SVector{3,Float64},
                    SVector{3,Float64},
                    SVector{3,Float64},
                ),
            )
            @test length(allocs_vec) == 0
        end

        @testset "reconstruct_velocities_fg" begin
            allocs_vec = check_allocs(
                Lambert.reconstruct_velocities_fg,
                (Float64, Float64, Float64, SVector{3,Float64}, SVector{3,Float64}),
            )
            @test length(allocs_vec) == 0
        end
    end

    @testset "Solver Zero-Allocation Tests" begin
        problem = LambertProblem(μ_perf, r1_perf, r2_perf, tof_perf)

        zero_alloc_solvers = [
            IzzoSolver(),
            GaussSolver(),
            AroraSolver(),
            ValladoSolver(),
            RussellSolver(),
            GoodingSolver(),
            BattinSolver(),
        ]

        for solver in zero_alloc_solvers
            name = string(typeof(solver))
            @testset "$name" begin
                allocs = check_allocs(solve, (typeof(problem), typeof(solver)))
                @test length(allocs) == 0
            end
        end
    end

    @testset "Solver Parameter Variation — Zero Allocation" begin
        problem = LambertProblem(μ_perf, r1_perf, r2_perf, tof_perf)

        for (atol, rtol) in [(1e-4, 1e-6), (1e-6, 1e-8), (1e-8, 1e-10)]
            @testset "IzzoSolver atol=$atol rtol=$rtol" begin
                solver = IzzoSolver(atol = atol, rtol = rtol)
                allocs = check_allocs(solve, (typeof(problem), typeof(solver)))
                @test length(allocs) == 0
            end
        end
    end

    @testset "Prograde vs Retrograde — Zero Allocation" begin
        problem = LambertProblem(μ_perf, r1_perf, r2_perf, tof_perf)

        for solver in ROBUST_SOLVERS
            name = string(typeof(solver))
            for prograde in [true, false]
                @testset "$name prograde=$prograde" begin
                    solver_dir = typeof(solver)(prograde = prograde)
                    allocs = check_allocs(solve, (typeof(problem), typeof(solver_dir)))
                    @test length(allocs) == 0
                end
            end
        end
    end

    @testset "Analytical Jacobian — Zero Allocation" begin
        problem = LambertProblem(μ_perf, r1_perf, r2_perf, tof_perf)
        sol = solve(problem, IzzoSolver())

        @testset "lambert_jacobian" begin
            allocs = check_allocs(lambert_jacobian, (typeof(problem), typeof(sol)))
            @test length(allocs) == 0
        end

        @testset "_twobody_stm" begin
            allocs = check_allocs(
                Lambert._twobody_stm,
                (SVector{3,Float64}, SVector{3,Float64}, Float64, Float64),
            )
            @test length(allocs) == 0
        end

        @testset "_solve_kepler_universal" begin
            allocs = check_allocs(
                Lambert._solve_kepler_universal,
                (Float64, Float64, Float64, Float64, Float64),
            )
            @test length(allocs) == 0
        end

        @testset "_mu_orbit_sensitivity" begin
            allocs = check_allocs(
                Lambert._mu_orbit_sensitivity,
                (
                    SVector{3,Float64},
                    SVector{3,Float64},
                    Float64,
                    Float64,
                    SMatrix{3,3,Float64,9},
                    SMatrix{3,3,Float64,9},
                    SMatrix{3,3,Float64,9},
                    SMatrix{3,3,Float64,9},
                ),
            )
            @test length(allocs) == 0
        end
    end

    @testset "Multi-Revolution — Zero Allocation" begin
        for solver in MULTIREV_SOLVERS
            name = string(typeof(solver))
            for M in [1, 2]
                problem = LambertProblem(μ_perf, r1_perf, r2_perf, tof_perf * (3 + 2*M))
                @testset "$name M=$M" begin
                    solver_mr = typeof(solver)(M = M)
                    allocs = check_allocs(solve, (typeof(problem), typeof(solver_mr)))
                    @test length(allocs) == 0
                end
            end
        end
    end
end
