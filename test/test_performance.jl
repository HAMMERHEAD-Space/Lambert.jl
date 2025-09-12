# Standard test problem for non-allocation testing
const μ_earth = 3.986004418e5  # [km³/s²]
const r1_std = [15945.34, 0.0, 0.0]  # [km]
const r2_std = [12214.83899, 10249.46731, 0.0]  # [km]
const tof_std = 76.0 * 60  # [s]

using AllocCheck

@testset "Non-Allocating Performance Tests" begin
    
    @testset "Solver Non-Allocation Tests" begin
        problem = LambertProblem(μ_earth, r1_std, r2_std, tof_std)
        
        for solver in ALL_SOLVERS
            @testset "$(typeof(solver)) Non-Allocating" begin
                solution = solve(problem, solver)
                @test solution.retcode == :SUCCESS

                # Use AllocCheck to verify zero allocations
                # These tests are expected to fail initially until optimizations are complete
                allocs_vec = check_allocs(solve, (typeof(problem), typeof(solver)))
                @test length(allocs_vec) == 0

            end
        end
    end
    
    @testset "Problem Construction Non-Allocation" begin
        # Test non-allocation for problem construction
        allocs_vec = check_allocs(LambertProblem, (typeof(μ_earth), typeof(r1_std), typeof(r2_std), typeof(tof_std)))
        @test length(allocs_vec) == 0
    end
    
    @testset "Heuristic Algorithm Selection Non-Allocation" begin
        problem = LambertProblem(μ_earth, r1_std, r2_std, tof_std)
        
        allocs_vec = check_allocs(select_lambert_algorithm, (typeof(problem),))
        @test length(allocs_vec) == 0
    end
    
    @testset "Multi-Revolution Non-Allocation" begin
        
        for solver in MULTIREV_SOLVERS
            for M in [1, 2]
                # Use longer time of flight for multi-revolution cases to ensure feasible solutions
                problem = LambertProblem(μ_earth, r1_std, r2_std, tof_std * (3 + 2*M))
                @testset "$(typeof(solver)) M=$M Non-Allocating" begin
                    solver_multirev = typeof(solver)(M=M)

                    solution = solve(problem, solver_multirev)
                    # Some solvers may have convergence issues for multi-rev cases
                    @test solution.retcode in [:SUCCESS, :MAXIMUM_ITERATIONS]
                    
                    allocs_vec = check_allocs(solve, (typeof(problem), typeof(solver_multirev)))
                    @test length(allocs_vec) == 0
                end
            end
        end
    end
    
    @testset "Prograde vs Retrograde Non-Allocation" begin
        problem = LambertProblem(μ_earth, r1_std, r2_std, tof_std)
        
        # Test both prograde and retrograde for non-allocation
        for solver in ROBUST_SOLVERS
            for prograde in [true, false]
                @testset "$(typeof(solver)) Prograde=$prograde Non-Allocating" begin
                    solver_direction = typeof(solver)(prograde=prograde)

                    solution = solve(problem, solver_direction)
                    @test solution.retcode == :SUCCESS
                    
                    allocs_vec = check_allocs(solve, (typeof(problem), typeof(solver_direction)))
                    @test length(allocs_vec) == 0
                end
            end
        end
    end
    
    @testset "Utility Function Non-Allocation" begin
        # Test key utility functions that should be non-allocating
        
        @testset "Stumpff Functions Non-Allocation" begin
            allocs_vec = check_allocs(AstroProblemsLambert.c2, (Float64,))
            @test length(allocs_vec) == 0

            allocs_vec = check_allocs(AstroProblemsLambert.c3, (Float64,))
            @test length(allocs_vec) == 0
        end
        
        @testset "Convergence Check Non-Allocation" begin
            allocs_vec = check_allocs(AstroProblemsLambert.check_convergence, (Float64, Float64, Float64, Float64))
            @test length(allocs_vec) == 0
        end
        
        @testset "Geometry Functions Non-Allocation" begin
            allocs_vec = check_allocs(AstroProblemsLambert.compute_semiperimeter, (Float64, Float64, Float64))
            @test length(allocs_vec) == 0
        end

        @testset "Other Inline Utility Functions" begin
            # handle_max_iterations - should be non-allocating
            allocs_vec = check_allocs(AstroProblemsLambert.handle_max_iterations, (Int, Int))
            @test length(allocs_vec) == 0
        end
    end
    
    @testset "Solver Parameter Variation Non-Allocation" begin
        problem = LambertProblem(μ_earth, r1_std, r2_std, tof_std)
        
        # Test with different tolerance settings  
        tolerance_settings = [
            (1e-4, 1e-6),   # Relaxed
            (1e-6, 1e-8),   # Standard  
            (1e-8, 1e-10),  # Tight
        ]
        
        for (atol, rtol) in tolerance_settings
            @testset "Tolerance atol=$atol rtol=$rtol Non-Allocating" begin
                solver = IzzoSolver(atol=atol, rtol=rtol)
                
                solution = solve(problem, solver)
                @test solution.retcode == :SUCCESS

                allocs_vec = check_allocs(solve, (typeof(problem), typeof(solver)))
                @test length(allocs_vec) == 0
            end
        end
    end
    
    @testset "Specialized Solver Functions Non-Allocation" begin
        # Test solver-specific utility functions that should be non-allocating
        
        @testset "Izzo Solver Utilities" begin
            # Test compute_y function from Izzo solver
            allocs_vec = check_allocs(AstroProblemsLambert.compute_y, (Float64, Float64))
            @test length(allocs_vec) == 0
        
            allocs_vec = check_allocs(AstroProblemsLambert.reconstruct, (Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64))
            @test length(allocs_vec) == 0
        end
    end
        
    @testset "Vector Reconstruction Functions" begin

        allocs_vec = check_allocs(AstroProblemsLambert.reconstruct_velocities_from_components, (Float64, Float64, Float64, Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}))
        @test length(allocs_vec) == 0
        
        allocs_vec = check_allocs(AstroProblemsLambert.reconstruct_velocities_fg, (Float64, Float64, Float64, Vector{Float64}, Vector{Float64}))
        @test length(allocs_vec) == 0
    end
end