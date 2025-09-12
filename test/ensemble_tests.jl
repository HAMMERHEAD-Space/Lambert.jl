@testset "EnsembleProblem Integration Tests" begin
    
    @testset "LambertProblem SciMLBase Interface" begin
        # Test basic SciMLBase inheritance
        μ_earth = 3.986004418e5  # [km³/s²]
        r1 = [15945.34, 0.0, 0.0]  # [km]
        r2 = [12214.83899, 10249.46731, 0.0]  # [km]
        tof = 76.0 * 60  # [s]
        
        problem = LambertProblem(μ_earth, r1, r2, tof)
        
        @test problem isa SciMLBase.AbstractSciMLProblem
        @test problem isa LambertProblem
        @test problem.μ == μ_earth
        @test problem.r1 == r1
        @test problem.r2 == r2
        @test problem.tof == tof
    end
    
    @testset "remake Function Tests" begin
        # Test remake functionality
        μ_earth = 3.986004418e5
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 7000.0, 0.0]
        tof = 3600.0
        
        original_prob = LambertProblem(μ_earth, r1, r2, tof)
        
        # Test single parameter remake
        new_tof = 7200.0
        remade_prob = remake(original_prob, tof=new_tof)
        @test remade_prob.tof == new_tof
        @test remade_prob.μ == original_prob.μ
        @test remade_prob.r1 == original_prob.r1
        @test remade_prob.r2 == original_prob.r2
        
        # Test multiple parameter remake
        new_μ = 3.986004418e5 * 1.1
        new_r1 = [7500.0, 0.0, 0.0]
        remade_prob2 = remake(original_prob, μ=new_μ, r1=new_r1, tof=new_tof)
        @test remade_prob2.μ == new_μ
        @test remade_prob2.r1 == new_r1
        @test remade_prob2.r2 == original_prob.r2
        @test remade_prob2.tof == new_tof
        
        # Test remake with no changes (should return equivalent problem)
        unchanged_prob = remake(original_prob)
        @test unchanged_prob.μ == original_prob.μ
        @test unchanged_prob.r1 == original_prob.r1
        @test unchanged_prob.r2 == original_prob.r2
        @test unchanged_prob.tof == original_prob.tof
    end
    
    @testset "EnsembleProblem Creation and Basic Properties" begin
        # Create base Lambert problem
        μ_earth = 3.986004418e5
        r1 = [15945.34, 0.0, 0.0]
        r2 = [12214.83899, 10249.46731, 0.0]
        tof = 76.0 * 60
        
        base_prob = LambertProblem(μ_earth, r1, r2, tof)
        
        # Define problem function for ensemble
        function prob_func(prob, i, repeat)
            # Vary time of flight by ±10% with some systematic variation
            tof_variation = 1.0 + 0.1 * sin(2π * i / 10)  # Deterministic variation for testing
            new_tof = prob.tof * tof_variation
            return remake(prob, tof=new_tof)
        end
        
        # Create ensemble problem
        ensemble_prob = EnsembleProblem(base_prob, prob_func=prob_func)
        
        @test ensemble_prob isa EnsembleProblem
        @test ensemble_prob.prob == base_prob
        @test ensemble_prob.prob_func == prob_func
        
        # Test that prob_func works correctly
        test_prob = prob_func(base_prob, 1, 1)
        @test test_prob isa LambertProblem
        @test test_prob.μ == base_prob.μ
        @test test_prob.r1 == base_prob.r1
        @test test_prob.r2 == base_prob.r2
        @test test_prob.tof ≠ base_prob.tof  # Should be different due to variation
    end
    
    @testset "EnsembleProblem Parameter Variations" begin
        μ_earth = 3.986004418e5
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 7000.0, 0.0]
        tof = 3600.0
        
        base_prob = LambertProblem(μ_earth, r1, r2, tof)
        
        @testset "Time of Flight Variations" begin
            # Test systematic time of flight variations
            function tof_variation_func(prob, i, repeat)
                # Create variations from 90% to 110% of original tof
                factor = 0.9 + 0.2 * (i - 1) / 9  # For i = 1 to 10
                return remake(prob, tof=prob.tof * factor)
            end
            
            ensemble_prob = EnsembleProblem(base_prob, prob_func=tof_variation_func)
            
            # Test that different indices produce different problems
            prob1 = tof_variation_func(base_prob, 1, 1)
            prob5 = tof_variation_func(base_prob, 5, 1)
            prob10 = tof_variation_func(base_prob, 10, 1)
            
            @test prob1.tof ≈ base_prob.tof * 0.9
            @test prob5.tof ≈ base_prob.tof * (0.9 + 0.2 * 4 / 9)  # Middle value
            @test prob10.tof ≈ base_prob.tof * 1.1
            
            # Ensure other parameters remain unchanged
            for prob in [prob1, prob5, prob10]
                @test prob.μ == base_prob.μ
                @test prob.r1 == base_prob.r1
                @test prob.r2 == base_prob.r2
            end
        end
        
        @testset "Gravitational Parameter Variations" begin
            # Test μ variations (e.g., for different celestial bodies)
            μ_values = [3.986004418e5, 1.327124400e11, 4.902800e4]  # Earth, Sun, Mars
            
            function μ_variation_func(prob, i, repeat)
                μ_index = ((i - 1) % length(μ_values)) + 1
                return remake(prob, μ=μ_values[μ_index])
            end
            
            ensemble_prob = EnsembleProblem(base_prob, prob_func=μ_variation_func)
            
            # Test different μ values
            prob1 = μ_variation_func(base_prob, 1, 1)
            prob2 = μ_variation_func(base_prob, 2, 1)
            prob3 = μ_variation_func(base_prob, 3, 1)
            
            @test prob1.μ == μ_values[1]
            @test prob2.μ == μ_values[2]
            @test prob3.μ == μ_values[3]
        end
        
        @testset "Position Vector Variations" begin
            # Test small perturbations in position vectors
            Random.seed!(42)  # For reproducible results
            
            function position_variation_func(prob, i, repeat)
                # Add small random perturbations (1% of magnitude)
                r1_pert = prob.r1 .+ 0.01 * norm(prob.r1) * randn(3)
                r2_pert = prob.r2 .+ 0.01 * norm(prob.r2) * randn(3)
                return remake(prob, r1=r1_pert, r2=r2_pert)
            end
            
            ensemble_prob = EnsembleProblem(base_prob, prob_func=position_variation_func)
            
            # Test that positions are perturbed but close to original
            Random.seed!(42)  # Reset for consistent test
            prob_pert = position_variation_func(base_prob, 1, 1)
            
            @test !isapprox(prob_pert.r1, base_prob.r1, rtol=1e-6)
            @test !isapprox(prob_pert.r2, base_prob.r2, rtol=1e-6)
            @test isapprox(prob_pert.r1, base_prob.r1, rtol=0.02)  # Within 2%
            @test isapprox(prob_pert.r2, base_prob.r2, rtol=0.02)  # Within 2%
        end
        
        @testset "Multi-Parameter Variations" begin
            # Test simultaneous variation of multiple parameters
            function multi_param_func(prob, i, repeat)
                # Vary both tof and μ simultaneously
                tof_factor = 0.8 + 0.4 * rand()  # 80% to 120%
                μ_factor = 0.95 + 0.1 * rand()   # 95% to 105%
                
                return remake(prob, 
                            tof=prob.tof * tof_factor,
                            μ=prob.μ * μ_factor)
            end
            
            Random.seed!(123)
            ensemble_prob = EnsembleProblem(base_prob, prob_func=multi_param_func)
            
            Random.seed!(123)  # Reset for consistent test
            prob_varied = multi_param_func(base_prob, 1, 1)
            
            @test prob_varied.tof ≠ base_prob.tof
            @test prob_varied.μ ≠ base_prob.μ
            @test prob_varied.r1 == base_prob.r1  # Should remain unchanged
            @test prob_varied.r2 == base_prob.r2  # Should remain unchanged
        end
    end
    
    @testset "EnsembleProblem with Output Functions" begin
        μ_earth = 3.986004418e5
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 7000.0, 0.0]
        tof = 3600.0
        
        base_prob = LambertProblem(μ_earth, r1, r2, tof)
        
        function simple_prob_func(prob, i, repeat)
            return remake(prob, tof=prob.tof * (1.0 + 0.1 * (i - 1) / 9))
        end
        
        # Test with custom output function
        function custom_output_func(sol, i)
            # Extract just the velocity magnitudes for analysis
            if hasfield(typeof(sol), :v1) && hasfield(typeof(sol), :v2)
                return (v1_mag=norm(sol.v1), v2_mag=norm(sol.v2), index=i)
            else
                return (error="No velocity data", index=i)
            end
        end
        
        ensemble_prob = EnsembleProblem(base_prob, 
                                      prob_func=simple_prob_func,
                                      output_func=custom_output_func)
        
        @test ensemble_prob isa EnsembleProblem
        @test ensemble_prob.prob_func == simple_prob_func
        @test ensemble_prob.output_func == custom_output_func
    end
    
    @testset "Statistical Analysis Setup" begin
        # Demonstrate how to set up problems for statistical analysis
        μ_earth = 3.986004418e5
        r1 = [15945.34, 0.0, 0.0]
        r2 = [12214.83899, 10249.46731, 0.0]
        tof = 76.0 * 60
        
        base_prob = LambertProblem(μ_earth, r1, r2, tof)
        
        # Monte Carlo setup with random variations
        Random.seed!(456)
        function monte_carlo_func(prob, i, repeat)
            # Add uncertainties to all parameters
            μ_noise = prob.μ * (1.0 + 0.001 * randn())  # 0.1% uncertainty
            r1_noise = prob.r1 .+ 10.0 * randn(3)       # 10 km position uncertainty
            r2_noise = prob.r2 .+ 10.0 * randn(3)       # 10 km position uncertainty  
            tof_noise = prob.tof * (1.0 + 0.01 * randn()) # 1% time uncertainty
            
            return remake(prob, μ=μ_noise, r1=r1_noise, r2=r2_noise, tof=tof_noise)
        end
        
        ensemble_prob = EnsembleProblem(base_prob, prob_func=monte_carlo_func)
        
        # Generate a few sample problems to verify statistical properties
        Random.seed!(456)  # Reset for consistent results
        n_samples = 100
        tof_samples = Float64[]
        μ_samples = Float64[]
        
        for i in 1:n_samples
            sample_prob = monte_carlo_func(base_prob, i, 1)
            push!(tof_samples, sample_prob.tof)
            push!(μ_samples, sample_prob.μ)
        end
        
        # Verify statistical properties
        @test abs(mean(tof_samples) - base_prob.tof) < 0.01 * base_prob.tof  # Mean close to original
        @test std(tof_samples) > 0  # Non-zero standard deviation
        @test abs(mean(μ_samples) - base_prob.μ) < 0.001 * base_prob.μ  # Mean close to original
        @test std(μ_samples) > 0  # Non-zero standard deviation
    end
end
