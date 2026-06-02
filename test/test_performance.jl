using AllocCheck

const μ_perf = 3.986004418e5
const r1_perf = [15945.34, 0.0, 0.0]
const r2_perf = [12214.83899, 10249.46731, 0.0]
const tof_perf = 76.0 * 60

# AllocCheck reports spurious `jl_get_pgcstack_static` "allocating runtime
# call"s on macOS aarch64 with Julia 1.12+. These are not real heap allocations:
# the analyzed code is allocation-free on every other platform/version (Linux,
# Windows, and macOS on Julia 1.10/1.11). This is a known AllocCheck/Julia
# limitation, so the checks are skipped on the affected platform.
# Ref: https://github.com/SciML/SciMLStructures.jl/issues/59
const _SKIP_ALLOCCHECK = Sys.isapple() && Sys.ARCH === :aarch64 && VERSION >= v"1.12"

if _SKIP_ALLOCCHECK
    @info "Skipping AllocCheck allocation tests (spurious jl_get_pgcstack_static reports on macOS aarch64 + Julia 1.12+; see SciML/SciMLStructures.jl#59)."
end

function checked_allocs(f, types)
    _SKIP_ALLOCCHECK && return ()
    allocs = check_allocs(f, types)
    if !isempty(allocs)
        printstyled(stdout, "\n[ALLOC] "; color = :red, bold = true)
        println(stdout, f, " with ", types, " => ", length(allocs), " allocation(s)")
        for (i, a) in enumerate(allocs)
            println(stdout, "  ──────── allocation ", i, " ────────")
            show(stdout, MIME"text/plain"(), a)
            println(stdout)
        end
        flush(stdout)
    end
    return allocs
end

@testset "Non-Allocating Performance Tests" begin

    @testset "Problem Construction — Zero Allocation" begin
        @test length(
            checked_allocs(
                LambertProblem,
                (typeof(μ_perf), typeof(r1_perf), typeof(r2_perf), typeof(tof_perf)),
            ),
        ) == 0
    end

    @testset "Heuristic Algorithm Selection — Zero Allocation" begin
        problem = LambertProblem(μ_perf, r1_perf, r2_perf, tof_perf)
        @test length(checked_allocs(select_lambert_algorithm, (typeof(problem),))) == 0
    end

    @testset "Utility Functions — Zero Allocation" begin
        @testset "Stumpff c2" begin
            @test length(checked_allocs(Lambert.c2, (Float64,))) == 0
        end

        @testset "Stumpff c3" begin
            @test length(checked_allocs(Lambert.c3, (Float64,))) == 0
        end

        @testset "Scalar check_convergence" begin
            @test length(
                checked_allocs(
                    Lambert.check_convergence,
                    (Float64, Float64, Float64, Float64),
                ),
            ) == 0
        end

        @testset "compute_semiperimeter" begin
            @test length(
                checked_allocs(Lambert.compute_semiperimeter, (Float64, Float64, Float64)),
            ) == 0
        end

        @testset "handle_max_iterations" begin
            @test length(checked_allocs(Lambert.handle_max_iterations, (Int, Int))) == 0
        end

        @testset "assert_transfer_angle_not_zero" begin
            @test length(
                checked_allocs(Lambert.assert_transfer_angle_not_zero, (Float64,)),
            ) == 0
        end

        @testset "assert_transfer_angle_not_pi" begin
            @test length(
                checked_allocs(Lambert.assert_transfer_angle_not_pi, (Float64,)),
            ) == 0
        end
    end

    @testset "Izzo Solver Utilities — Zero Allocation" begin
        @test length(checked_allocs(Lambert.compute_y, (Float64, Float64))) == 0

        @test length(
            checked_allocs(
                Lambert.reconstruct,
                (Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64),
            ),
        ) == 0
    end

    @testset "Vector Reconstruction — Zero Allocation" begin
        @testset "reconstruct_velocities_from_components" begin
            @test length(
                checked_allocs(
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
                ),
            ) == 0
        end

        @testset "reconstruct_velocities_fg" begin
            @test length(
                checked_allocs(
                    Lambert.reconstruct_velocities_fg,
                    (Float64, Float64, Float64, SVector{3,Float64}, SVector{3,Float64}),
                ),
            ) == 0
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
                @test length(checked_allocs(solve, (typeof(problem), typeof(solver)))) == 0
            end
        end
    end

    @testset "Solver Parameter Variation — Zero Allocation" begin
        problem = LambertProblem(μ_perf, r1_perf, r2_perf, tof_perf)

        for (atol, rtol) in [(1e-4, 1e-6), (1e-6, 1e-8), (1e-8, 1e-10)]
            @testset "IzzoSolver atol=$atol rtol=$rtol" begin
                solver = IzzoSolver(atol = atol, rtol = rtol)
                @test length(checked_allocs(solve, (typeof(problem), typeof(solver)))) == 0
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
                    @test length(
                        checked_allocs(solve, (typeof(problem), typeof(solver_dir))),
                    ) == 0
                end
            end
        end
    end

    @testset "Analytical Jacobian — Zero Allocation" begin
        problem = LambertProblem(μ_perf, r1_perf, r2_perf, tof_perf)
        sol = solve(problem, IzzoSolver())

        @testset "lambert_jacobian" begin
            @test length(
                checked_allocs(lambert_jacobian, (typeof(problem), typeof(sol))),
            ) == 0
        end

        @testset "_twobody_stm" begin
            @test length(
                checked_allocs(
                    Lambert._twobody_stm,
                    (SVector{3,Float64}, SVector{3,Float64}, Float64, Float64),
                ),
            ) == 0
        end

        @testset "_solve_kepler_universal" begin
            @test length(
                checked_allocs(
                    Lambert._solve_kepler_universal,
                    (Float64, Float64, Float64, Float64, Float64),
                ),
            ) == 0
        end

        @testset "_mu_orbit_sensitivity" begin
            @test length(
                checked_allocs(
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
                ),
            ) == 0
        end
    end

    @testset "Multi-Revolution — Zero Allocation" begin
        for solver in MULTIREV_SOLVERS
            name = string(typeof(solver))
            for M in [1, 2]
                problem = LambertProblem(μ_perf, r1_perf, r2_perf, tof_perf * (3 + 2*M))
                @testset "$name M=$M" begin
                    solver_mr = typeof(solver)(M = M)
                    @test length(
                        checked_allocs(solve, (typeof(problem), typeof(solver_mr))),
                    ) == 0
                end
            end
        end
    end
end
