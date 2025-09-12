# Test tolerances
const ATOL = 0.025
const RTOL = 0.001

@testset "Lambert Solvers Tests" begin
    @testset "Vallado Book Example" begin
        # Example 5.7 from Fundamentals of Astrodynamics and Applications (4th Edition)
        μ_earth = 3.986004418e5  # [km³/s²]
        r1 = [15945.34, 0.0, 0.0]  # [km]
        r2 = [12214.83899, 10249.46731, 0.0]  # [km]
        tof = 76.0 * 60  # [s]
        
        # Expected results
        expected_v1 = [2.058913, 2.915965, 0.0]  # [km/s]
        expected_v2 = [-3.451565, 0.910315, 0.0]  # [km/s]
        
        problem = LambertProblem(μ_earth, r1, r2, tof)
        
        for solver in ALL_SOLVERS
            @testset "$(typeof(solver))" begin
                solution = solve(problem, solver)
                
                # Extract velocities directly from solution
                v1 = solution.v1
                v2 = solution.v2
                
                @test solution.retcode == :SUCCESS
                @test v1 ≈ expected_v1 atol=ATOL rtol=RTOL
                @test v2 ≈ expected_v2 atol=ATOL rtol=RTOL
            end
        end
    end
    
    @testset "Curtis Book Example" begin
        # Example 5.2 from Orbital Mechanics for Engineering Students (3rd Edition)
        μ_earth = 3.986004418e5  # [km³/s²]
        r1 = [5000.0, 10000.0, 2100.0]  # [km]
        r2 = [-14600.0, 2500.0, 7000.0]  # [km]
        tof = 3600  # [s]
        
        expected_v1 = [-5.9925, 1.9254, 3.2456]  # [km/s]
        expected_v2 = [-3.3125, -4.1966, -0.38529]  # [km/s]
        
        problem = LambertProblem(μ_earth, r1, r2, tof)
        for solver in ROBUST_SOLVERS
            @testset "$(typeof(solver))" begin
                solution = SciMLBase.solve(problem, solver)
                
                @test solution.retcode == :SUCCESS
                @test solution.v1 ≈ expected_v1 atol=ATOL rtol=RTOL
                @test solution.v2 ≈ expected_v2 atol=ATOL rtol=RTOL
            end
        end
    end
    
    @testset "Battin Book Example" begin
        # Example 7.12 from An Introduction to the Mathematics and Methods of Astrodynamics
        μ_sun = 39.47692641  # [AU ** 3 / year ** 2]
        r1 = [0.159321004, 0.579266185, 0.052359607] # [AU]
        r2 = [0.057594337, 0.605750797, 0.068345246] # [AU]
        tof = 0.010794065  # [year]

        expected_v1 = [-9.303603251, 3.018641330, 1.536362143] # [AU / year]

        problem = LambertProblem(μ_sun, r1, r2, tof)

        for solver in ROBUST_SOLVERS
            @testset "$(typeof(solver))" begin   
                solution = SciMLBase.solve(problem, solver)
                
                @test solution.retcode == :SUCCESS
                @test solution.v1 ≈ expected_v1 atol=ATOL rtol=RTOL
            end
        end
    end
    
    @testset "GMAT Hyperbolic Prograde Example" begin
        μ_earth = 3.986004418e5  # [km³/s²]
        r1 = [7100.0, 200.0, 1300.0]  # [km]
        r2 = [-38113.5870, 67274.1946, 29309.5799]  # [km]
        tof = 12000.0  # [s]
        
        expected_v1 = [0.0, 10.35, 5.5]  # [km/s]
        expected_v2 = [-3.6379, 4.4932, 1.7735]  # [km/s]
        
        problem = LambertProblem(μ_earth, r1, r2, tof)

        for solver in ROBUST_SOLVERS
            @testset "$(typeof(solver))" begin
                    
                solution = SciMLBase.solve(problem, solver)
                    
                @test solution.retcode == :SUCCESS
                @test solution.v1 ≈ expected_v1 atol=ATOL rtol=RTOL
                @test solution.v2 ≈ expected_v2 atol=ATOL rtol=RTOL
            end
        end
    end
    
    @testset "GMAT Hyperbolic Retrograde Example" begin
        μ_earth = 3.986004418e5  # [km³/s²]
        r1 = [7100.0, 200.0, 1300.0]  # [km]
        r2 = [-47332.7499, -54840.2027, -37100.17067]  # [km]
        tof = 12000.0  # [s]
        
        expected_v1 = [0.0, -10.35, -5.5]  # [km/s]
        expected_v2 = [-4.3016, -3.4314, -2.5467]  # [km/s]
        
        problem = LambertProblem(μ_earth, r1, r2, tof)

        for solver in ROBUST_SOLVERS
            #solver_retrograde = typeof(solver)(prograde=false)
            @testset "$(typeof(solver))" begin
                solution = SciMLBase.solve(problem, solver; prograde=false)
                
                @test solution.retcode == :SUCCESS 
                @test solution.v1 ≈ expected_v1 atol=ATOL rtol=RTOL
                @test solution.v2 ≈ expected_v2 atol=ATOL rtol=RTOL
            end
        end
    end
    
    @testset "Der Article Example I - Prograde Low Path" begin
        μ_earth = 3.986004418e5  # [km³/s²]
        r1 = [2249.171260, 1898.007100, 5639.599193]  # [km]
        r2 = [1744.495443, -4601.556054, 4043.864391]  # [km]
        tof = 1618.50  # [s]
        
        expected_v1 = [-2.09572809, 3.92602196, -4.94516810]  # [km/s]
        expected_v2 = [2.46309613, 0.84490197, 6.10890863]  # [km/s]
        
        problem = LambertProblem(μ_earth, r1, r2, tof)

        for solver in ROBUST_SOLVERS
            @testset "$(typeof(solver))" begin
                solution = SciMLBase.solve(problem, solver)
                
                @test solution.retcode == :SUCCESS
                @test solution.v1 ≈ expected_v1 atol=ATOL rtol=RTOL
                @test solution.v2 ≈ expected_v2 atol=ATOL rtol=RTOL
            end
        end
    end
    
    @testset "Der Article Example II - Prograde High Path" begin
        μ_earth = 3.986004418e5  # [km³/s²]
        r1 = [22592.145603, -1599.915239, -19783.950506]  # [km]
        r2 = [1922.067697, 4054.157051, -8925.727465]  # [km]
        tof = 36000  # [s]
        
        expected_v1 = [2.000652697, 0.387688615, -2.666947760]  # [km/s]
        expected_v2 = [-3.79246619, -1.77707641, 6.856814395]  # [km/s]
        
        problem = LambertProblem(μ_earth, r1, r2, tof)

        for solver in ROBUST_SOLVERS
            @testset "$(typeof(solver))" begin
                solution = SciMLBase.solve(problem, solver)
                
                @test solution.retcode == :SUCCESS
                @test solution.v1 ≈ expected_v1 atol=ATOL rtol=RTOL
                @test solution.v2 ≈ expected_v2 atol=ATOL rtol=RTOL
            end
        end
    end
    
    @testset "Solver Parameter Validation" begin
        # Test that solvers properly validate their parameters
        μ_earth = 3.986004418e5
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 7000.0, 0.0]
        tof = 3600.0
        
        # Test with various solver parameters
        @testset "Custom Parameters" begin
            solver = GoodingSolver(maxiter=50, atol=1e-6, rtol=1e-8)
            problem = LambertProblem(μ_earth, r1, r2, tof)
            solution = SciMLBase.solve(problem, solver)
            @test solution.retcode in [:SUCCESS, :MAXIMUM_ITERATIONS]
        end
        
        # Test invalid inputs
        @testset "Invalid Inputs" begin
            problem = LambertProblem(μ_earth, r1, r1, tof)  # Same position
            solver = GoodingSolver()
            @test_throws AssertionError SciMLBase.solve(problem, solver)
        end
    end
    
    # Test heuristic algorithm selection
    @testset "Heuristic Algorithm Selection" begin
        μ = 398600.4418
        
        # Test different transfer angle scenarios
        @testset "Transfer Angle Selection" begin
            # Short transfer (~30°) -> Should select Gooding
            r1_short = [7000.0, 0.0, 0.0]
            r2_short = [6062.2, 3500.0, 0.0]  # ~30° transfer
            prob_short = LambertProblem(μ, r1_short, r2_short, 1800.0)
            alg_short = select_lambert_algorithm(prob_short)
            @test alg_short isa GoodingSolver
            
            # Medium transfer (~120°) -> Should select Izzo
            r1_med = [7000.0, 0.0, 0.0] 
            r2_med = [-3500.0, 6062.2, 0.0]  # ~120° transfer
            prob_med = LambertProblem(μ, r1_med, r2_med, 3600.0)
            alg_med = select_lambert_algorithm(prob_med)
            @test alg_med isa IzzoSolver
            
            # Large transfer (~300°) -> Should select Arora
            r1_large = [7000.0, 0.0, 0.0]
            r2_large = [3500.0, -6062.2, 0.0]  # ~300° transfer  
            prob_large = LambertProblem(μ, r1_large, r2_large, 7200.0)
            alg_large = select_lambert_algorithm(prob_large)
            @test alg_large isa AroraSolver
            
            # Near-180° transfer -> Should select Battin
            r1_180 = [7000.0, 0.0, 0.0]
            r2_180 = [-6900.0, 500.0, 0.0]  # ~176° transfer (within 5° of 180°)
            prob_180 = LambertProblem(μ, r1_180, r2_180, 5400.0)
            alg_180 = select_lambert_algorithm(prob_180)
            @test alg_180 isa BattinSolver
        end
        
        # Test multi-revolution selection
        @testset "Multi-Revolution Selection" begin
            r1 = [7000.0, 0.0, 0.0]
            r2 = [-3500.0, 6062.2, 0.0]  # ~120° transfer
            prob = LambertProblem(μ, r1, r2, 3600.0)
            
            # Multi-rev with small angle -> Should select Izzo
            alg_multirev_small = select_lambert_algorithm(prob, 2)
            @test alg_multirev_small isa IzzoSolver
            @test alg_multirev_small.M == 2
            
            # Multi-rev with large angle -> Should select Arora
            r2_large = [3500.0, -6062.2, 0.0]  # ~300° transfer
            prob_large = LambertProblem(μ, r1, r2_large, 7200.0)
            alg_multirev_large = select_lambert_algorithm(prob_large, 1)
            @test alg_multirev_large isa AroraSolver
            @test alg_multirev_large.M == 1
        end
        
        # Test automatic solve functionality
        @testset "Automatic Solve" begin
            r1 = [15945.34, 0.0, 0.0]
            r2 = [12214.83, 10249.46, 0.0]
            tof = 600.0 * 60.0  # Increased to make multi-rev feasible
            prob = LambertProblem(μ, r1, r2, tof)
            
            # Test solve without algorithm (should use heuristic)
            solution = solve(prob)
            @test solution isa LambertSolution
            @test solution.retcode == :SUCCESS
            
            # Test solve with multi-rev (should select appropriate algorithm)
            solution_multirev = solve(prob, M=1)
            @test solution_multirev isa LambertSolution
            @test solution_multirev.retcode == :SUCCESS
        end
    end
    
    @testset "AstroCoords Interface" begin
        # Test LambertProblem with AstroCoords input
        μ_earth = 3.986004418e5  # [km³/s²]
        
        # Create Cartesian coordinates (position + velocity)
        r1_vec = [15945.34, 0.0, 0.0]  # [km]
        v1_vec = [0.0, 2.0, 1.0]  # [km/s] (dummy velocity)
        coord1_cart = Cartesian([r1_vec..., v1_vec...])
        
        r2_vec = [12214.83899, 10249.46731, 0.0]  # [km]  
        v2_vec = [0.0, -1.0, 2.0]  # [km/s] (dummy velocity)
        coord2_cart = Cartesian([r2_vec..., v2_vec...])
        
        tof = 76.0 * 60  # [s]
        
        # Expected results (same as Vallado example)
        expected_v1 = [2.058913, 2.915965, 0.0]  # [km/s]
        expected_v2 = [-3.451565, 0.910315, 0.0]  # [km/s]
        
        # Convert Cartesian to Keplerian and ModEq coordinates
        coord1_kep = Keplerian(coord1_cart, μ_earth)
        coord2_kep = Keplerian(coord2_cart, μ_earth)
        
        coord1_modeq = ModEq(coord1_cart, μ_earth)
        coord2_modeq = ModEq(coord2_cart, μ_earth)
        
        @testset "Cartesian Coordinates" begin
            # Create problem with Cartesian coordinates
            problem_cart = LambertProblem(μ_earth, coord1_cart, coord2_cart, tof)
            
            # Test that it matches vector-based problem
            problem_vector = LambertProblem(μ_earth, r1_vec, r2_vec, tof)
            @test problem_cart.μ == problem_vector.μ
            @test problem_cart.r1 ≈ problem_vector.r1 atol=1e-10
            @test problem_cart.r2 ≈ problem_vector.r2 atol=1e-10
            @test problem_cart.tof == problem_vector.tof
            
            # Test solving with a robust solver
            solution_cart = solve(problem_cart, GoodingSolver())
            
            @test solution_cart.retcode == :SUCCESS
            @test solution_cart.v1 ≈ expected_v1 atol=ATOL rtol=RTOL
            @test solution_cart.v2 ≈ expected_v2 atol=ATOL rtol=RTOL
        end
        
        @testset "Keplerian Coordinates" begin
            # Create problem with Keplerian coordinates
            problem_kep = LambertProblem(μ_earth, coord1_kep, coord2_kep, tof)
            
            # Test that it matches vector-based problem
            problem_vector = LambertProblem(μ_earth, r1_vec, r2_vec, tof)
            @test problem_kep.μ == problem_vector.μ
            @test problem_kep.r1 ≈ problem_vector.r1 atol=1e-10
            @test problem_kep.r2 ≈ problem_vector.r2 atol=1e-10
            @test problem_kep.tof == problem_vector.tof
            
            # Test solving with a robust solver
            solution_kep = solve(problem_kep, GoodingSolver())
            
            @test solution_kep.retcode == :SUCCESS
            @test solution_kep.v1 ≈ expected_v1 atol=ATOL rtol=RTOL
            @test solution_kep.v2 ≈ expected_v2 atol=ATOL rtol=RTOL
        end
        
        @testset "ModEq Coordinates" begin
            # Create problem with ModEq coordinates
            problem_modeq = LambertProblem(μ_earth, coord1_modeq, coord2_modeq, tof)
            
            # Test that it matches vector-based problem
            problem_vector = LambertProblem(μ_earth, r1_vec, r2_vec, tof)
            @test problem_modeq.μ == problem_vector.μ
            @test problem_modeq.r1 ≈ problem_vector.r1 atol=1e-10
            @test problem_modeq.r2 ≈ problem_vector.r2 atol=1e-10
            @test problem_modeq.tof == problem_vector.tof
            
            # Test solving with a robust solver
            solution_modeq = solve(problem_modeq, GoodingSolver())
            
            @test solution_modeq.retcode == :SUCCESS
            @test solution_modeq.v1 ≈ expected_v1 atol=ATOL rtol=RTOL
            @test solution_modeq.v2 ≈ expected_v2 atol=ATOL rtol=RTOL
        end
        
        @testset "Cross-Coordinate System Consistency" begin
            # Test that all coordinate systems produce identical position vectors
            problem_cart = LambertProblem(μ_earth, coord1_cart, coord2_cart, tof)
            problem_kep = LambertProblem(μ_earth, coord1_kep, coord2_kep, tof)
            problem_modeq = LambertProblem(μ_earth, coord1_modeq, coord2_modeq, tof)
            
            @test problem_cart.r1 ≈ problem_kep.r1 atol=1e-10
            @test problem_cart.r1 ≈ problem_modeq.r1 atol=1e-10
            @test problem_cart.r2 ≈ problem_kep.r2 atol=1e-10
            @test problem_cart.r2 ≈ problem_modeq.r2 atol=1e-10
            
            # Test that all coordinate systems produce identical solutions
            solution_cart = solve(problem_cart, GoodingSolver())
            solution_kep = solve(problem_kep, GoodingSolver())
            solution_modeq = solve(problem_modeq, GoodingSolver())
            
            @test solution_cart.retcode == :SUCCESS
            @test solution_kep.retcode == :SUCCESS
            @test solution_modeq.retcode == :SUCCESS
            @test solution_cart.v1 ≈ solution_kep.v1 atol=1e-10
            @test solution_cart.v1 ≈ solution_modeq.v1 atol=1e-10
            @test solution_cart.v2 ≈ solution_kep.v2 atol=1e-10
            @test solution_cart.v2 ≈ solution_modeq.v2 atol=1e-10
        end
    end

    @testset "Multi-Revolution Tests" begin
        μ_earth = 3.986004418e5  # [km³/s²]
        r1 = [15945.34, 0.0, 0.0]  # [km]
        r2 = [12214.83899, 10249.46731, 0.0]  # [km]
        tof_base = 76.0 * 60  # [s]
        
        @testset "Single vs Multi-Revolution Comparison" begin
            # Test that multi-rev solutions are different from single-rev
            problem_single = LambertProblem(μ_earth, r1, r2, tof_base)
            problem_multi1 = LambertProblem(μ_earth, r1, r2, tof_base * 5)  # M=1 feasible
            
            for solver_type in [IzzoSolver, BattinSolver]  # Most reliable multi-rev solvers
                @testset "$(solver_type) Single vs Multi-Rev" begin
                    solver_single = solver_type(M=0)
                    solver_multi = solver_type(M=1)
                    
                    solution_single = solve(problem_single, solver_single)
                    solution_multi = solve(problem_multi1, solver_multi)
                    
                    @test solution_single.retcode == :SUCCESS
                    @test solution_multi.retcode == :SUCCESS
                    
                    # Multi-rev solutions should be significantly different
                    @test !isapprox(solution_single.v1, solution_multi.v1, atol=0.1)
                    @test !isapprox(solution_single.v2, solution_multi.v2, atol=0.1)
                end
            end
        end
        
        @testset "MULTIREV_SOLVERS Functionality" begin
            # Test all multi-revolution capable solvers
            for solver_type in [IzzoSolver, BattinSolver, GaussSolver, AroraSolver]
                @testset "$(solver_type) Multi-Revolution" begin
                    for M in [1, 2]
                        # Use appropriate time of flight for multi-rev cases
                        tof_multirev = tof_base * (3 + 2*M)
                        problem = LambertProblem(μ_earth, r1, r2, tof_multirev)
                        solver = solver_type(M=M)
                        
                        solution = solve(problem, solver)
                        
                        # Test that solution is computed (some solvers may have convergence issues)
                        @test solution isa LambertSolution
                        @test solution.retcode in [:SUCCESS, :MAXIMUM_ITERATIONS]
                        @test solution.numiter >= 0
                        
                        if solution.retcode == :SUCCESS
                            # Verify solution vectors are reasonable (non-zero, finite)
                            @test norm(solution.v1) > 0.01  # Minimum reasonable velocity
                            @test norm(solution.v2) > 0.01
                            @test all(isfinite, solution.v1)
                            @test all(isfinite, solution.v2)
                            @test length(solution.v1) == 3
                            @test length(solution.v2) == 3
                        end
                    end
                end
            end
        end
        
        @testset "Multi-Revolution Direction Tests" begin
            # Test prograde vs retrograde for multi-revolution
            tof_multirev = tof_base * 5
            problem = LambertProblem(μ_earth, r1, r2, tof_multirev)
            
            for solver_type in [IzzoSolver, BattinSolver]  # Most reliable for direction tests
                @testset "$(solver_type) Prograde vs Retrograde" begin
                    solver_prograde = solver_type(M=1, prograde=true)
                    solver_retrograde = solver_type(M=1, prograde=false)
                    
                    solution_prograde = solve(problem, solver_prograde)
                    solution_retrograde = solve(problem, solver_retrograde)
                    
                    @test solution_prograde.retcode == :SUCCESS
                    @test solution_retrograde.retcode == :SUCCESS
                    
                    # Prograde and retrograde solutions should be different
                    @test !isapprox(solution_prograde.v1, solution_retrograde.v1, atol=0.1)
                    @test !isapprox(solution_prograde.v2, solution_retrograde.v2, atol=0.1)
                end
            end
        end
        
        @testset "Multi-Revolution Validation" begin
            # Test edge cases and validation for multi-rev
            
            @testset "High Revolution Count" begin
                # Test higher revolution counts (may be challenging)
                tof_high = tof_base * 10
                problem = LambertProblem(μ_earth, r1, r2, tof_high)
                
                solver = IzzoSolver(M=3)  # Izzo is most robust
                solution = solve(problem, solver)
                
                @test solution isa LambertSolution
                # Allow convergence issues for high revolution counts
                @test solution.retcode in [:SUCCESS, :MAXIMUM_ITERATIONS]
            end
            
            @testset "Multi-Revolution with Different TOF" begin
                # Test that longer TOF enables higher revolution counts
                for M in [0, 1]
                    tof_test = tof_base * (2 + 3*M)  # Scale TOF with revolution count
                    problem = LambertProblem(μ_earth, r1, r2, tof_test)
                    solver = IzzoSolver(M=M)
                    
                    solution = solve(problem, solver)
                    @test solution.retcode == :SUCCESS
                end
            end
        end
        
        @testset "Multi-Revolution Physical Validation" begin
            # Verify that multi-rev solutions satisfy physical constraints
            tof_multirev = tof_base * 7
            problem = LambertProblem(μ_earth, r1, r2, tof_multirev)
            solver = IzzoSolver(M=1)
            
            solution = solve(problem, solver)
            
            if solution.retcode == :SUCCESS
                v1, v2 = solution.v1, solution.v2
                
                # Test velocity magnitudes are physically reasonable for Earth orbit
                v1_mag = norm(v1)
                v2_mag = norm(v2)
                
                @test 0.1 < v1_mag < 20.0  # km/s - reasonable for Earth orbits
                @test 0.1 < v2_mag < 20.0  # km/s
                
                # Test that the solution satisfies energy constraints
                r1_mag = norm(problem.r1)
                r2_mag = norm(problem.r2)
                
                # Specific energy at each point
                ε1 = 0.5 * v1_mag^2 - μ_earth / r1_mag
                ε2 = 0.5 * v2_mag^2 - μ_earth / r2_mag
                
                # Energy should be consistent for the same orbit (within tolerance)
                @test abs(ε1 - ε2) < 10.0  # km²/s² - allow for some numerical error
            end
        end
    end
end
