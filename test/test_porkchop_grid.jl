@testset "Porkchop Grid Generation" begin
    # Earth's gravitational parameter
    μ_earth = 3.986004418e5  # [km³/s²]

    # Define simple circular orbit state functions returning [r; v]
    function state1_func(t)
        # Circular orbit at 7000 km
        ω = sqrt(μ_earth / 7000^3)  # Angular velocity
        θ = ω * t
        r = [7000 * cos(θ), 7000 * sin(θ), 0.0]
        v = [-7000 * ω * sin(θ), 7000 * ω * cos(θ), 0.0]
        return [r; v]  # Return 6-element vector [r; v]
    end

    function state2_func(t)
        # Circular orbit at 10000 km
        ω = sqrt(μ_earth / 10000^3)
        θ = ω * t + π / 4  # Phase offset
        r = [10000 * cos(θ), 10000 * sin(θ), 0.0]
        v = [-10000 * ω * sin(θ), 10000 * ω * cos(θ), 0.0]
        return [r; v]
    end

    # Time ranges (in seconds) - using small ranges for fast testing
    departure_times = range(0, 3600, length = 5)
    arrival_times = range(3600, 7200, length = 5)

    @testset "PorkchopGrid struct creation and basic functionality" begin
        # Test basic grid generation with default parameters
        grid =
            porkchop_grid(μ_earth, state1_func, state2_func, departure_times, arrival_times)

        @test grid isa PorkchopGrid
        @test grid.departure_times == departure_times
        @test grid.arrival_times == arrival_times
        @test grid.quantities == [:total_dv]
        @test haskey(grid.data, :total_dv)
        @test size(grid[:total_dv]) == (length(arrival_times), length(departure_times))

        # Test convenience methods
        @test haskey(grid, :total_dv)
        @test !haskey(grid, :nonexistent_quantity)
        @test :total_dv in keys(grid)
        @test grid[:total_dv] isa Matrix{Float64}
    end

    @testset "Multiple quantities" begin
        # Test with multiple quantities
        quantities = [:total_dv, :dv_departure, :dv_arrival]
        grid = porkchop_grid(
            μ_earth,
            state1_func,
            state2_func,
            departure_times,
            arrival_times;
            quantities = quantities,
        )

        @test grid.quantities == quantities
        @test length(grid.data) == 3

        for qty in quantities
            @test haskey(grid, qty)
            @test size(grid[qty]) == (length(arrival_times), length(departure_times))
        end

        # Test that all matrices have the same dimensions
        total_dv = grid[:total_dv]
        dv_dep = grid[:dv_departure]
        dv_arr = grid[:dv_arrival]

        @test size(total_dv) == size(dv_dep) == size(dv_arr)

        # Test that total_dv ≈ dv_departure + dv_arrival (where both are finite)
        for i in eachindex(total_dv)
            if isfinite(total_dv[i]) && isfinite(dv_dep[i]) && isfinite(dv_arr[i])
                @test total_dv[i] ≈ dv_dep[i] + dv_arr[i] rtol=1e-10
            end
        end
    end

    @testset "Excess velocity quantities" begin
        # Test excess velocity aliases
        quantities =
            [:total_excess_velocity, :excess_velocity_departure, :excess_velocity_arrival]
        grid = porkchop_grid(
            μ_earth,
            state1_func,
            state2_func,
            departure_times,
            arrival_times;
            quantities = quantities,
        )

        @test grid.quantities == quantities
        @test haskey(grid, :total_excess_velocity)
        @test haskey(grid, :excess_velocity_departure)
        @test haskey(grid, :excess_velocity_arrival)

        # Test that excess velocity aliases work the same as dv quantities
        grid_dv = porkchop_grid(
            μ_earth,
            state1_func,
            state2_func,
            departure_times,
            arrival_times;
            quantities = [:total_dv, :dv_departure, :dv_arrival],
        )

        # Compare element-wise, handling NaN values properly
        total_excess = grid[:total_excess_velocity]
        total_dv = grid_dv[:total_dv]
        for i in eachindex(total_excess)
            if isnan(total_excess[i]) && isnan(total_dv[i])
                continue  # Both NaN, that's fine
            elseif isnan(total_excess[i]) || isnan(total_dv[i])
                @test false  # One is NaN, the other isn't
            else
                @test total_excess[i] ≈ total_dv[i]
            end
        end

        # Same for departure
        dep_excess = grid[:excess_velocity_departure]
        dep_dv = grid_dv[:dv_departure]
        for i in eachindex(dep_excess)
            if isnan(dep_excess[i]) && isnan(dep_dv[i])
                continue
            elseif isnan(dep_excess[i]) || isnan(dep_dv[i])
                @test false
            else
                @test dep_excess[i] ≈ dep_dv[i]
            end
        end

        # Same for arrival
        arr_excess = grid[:excess_velocity_arrival]
        arr_dv = grid_dv[:dv_arrival]
        for i in eachindex(arr_excess)
            if isnan(arr_excess[i]) && isnan(arr_dv[i])
                continue
            elseif isnan(arr_excess[i]) || isnan(arr_dv[i])
                @test false
            else
                @test arr_excess[i] ≈ arr_dv[i]
            end
        end
    end

    @testset "Different solvers" begin
        # Test with different Lambert solvers
        for solver in ROBUST_SOLVERS
            grid = porkchop_grid(
                μ_earth,
                state1_func,
                state2_func,
                departure_times,
                arrival_times;
                solver = solver,
            )

            @test grid isa PorkchopGrid
            @test size(grid[:total_dv]) == (length(arrival_times), length(departure_times))

            # Should have some valid (finite) solutions
            @test any(isfinite, grid[:total_dv])
        end
    end

    @testset "Ensemble algorithms" begin
        # Test with different ensemble algorithms
        ensemble_algs = [SciMLBase.EnsembleSerial()]

        # Only test threaded if we have multiple threads
        if Threads.nthreads() > 1
            push!(ensemble_algs, SciMLBase.EnsembleThreads())
        end

        for ensemble_alg in ensemble_algs
            grid = porkchop_grid(
                μ_earth,
                state1_func,
                state2_func,
                departure_times,
                arrival_times;
                ensemble_alg = ensemble_alg,
            )

            @test grid isa PorkchopGrid
            @test size(grid[:total_dv]) == (length(arrival_times), length(departure_times))
        end
    end

    @testset "Max ΔV filtering" begin
        # Test max_deltav parameter
        max_dv = 5.0  # km/s
        grid = porkchop_grid(
            μ_earth,
            state1_func,
            state2_func,
            departure_times,
            arrival_times;
            max_deltav = max_dv,
        )

        # All finite values should be <= max_dv
        finite_values = filter(isfinite, grid[:total_dv])
        if !isempty(finite_values)
            @test all(v -> v <= max_dv, finite_values)
        end
    end

    @testset "Edge cases" begin
        # Test with no valid trajectories (all negative TOF)
        bad_departure_times = range(7200, 10800, length = 3)  # All after arrival times
        bad_arrival_times = range(0, 3600, length = 3)

        grid = porkchop_grid(
            μ_earth,
            state1_func,
            state2_func,
            bad_departure_times,
            bad_arrival_times,
        )

        @test grid isa PorkchopGrid
        @test all(isnan, grid[:total_dv])

        # Test with single time point
        single_dep = [1800.0]
        single_arr = [5400.0]

        grid_single =
            porkchop_grid(μ_earth, state1_func, state2_func, single_dep, single_arr)

        @test grid_single isa PorkchopGrid
        @test size(grid_single[:total_dv]) == (1, 1)
    end

    @testset "Invalid quantity error handling" begin
        # Test error for invalid quantity
        @test_throws ErrorException porkchop_grid(
            μ_earth,
            state1_func,
            state2_func,
            departure_times,
            arrival_times;
            quantities = [:invalid_quantity],
        )

        # Test mixed valid and invalid quantities
        @test_throws ErrorException porkchop_grid(
            μ_earth,
            state1_func,
            state2_func,
            departure_times,
            arrival_times;
            quantities = [:total_dv, :invalid_quantity],
        )
    end

    @testset "Data consistency" begin
        # Generate grid and verify data consistency
        grid = porkchop_grid(
            μ_earth,
            state1_func,
            state2_func,
            departure_times,
            arrival_times;
            quantities = [:total_dv, :dv_departure, :dv_arrival],
        )

        # Check that matrices have correct dimensions
        n_dep = length(departure_times)
        n_arr = length(arrival_times)

        for qty in grid.quantities
            matrix = grid[qty]
            @test size(matrix) == (n_arr, n_dep)
            @test eltype(matrix) == Float64
        end

        # Check that departure and arrival times are preserved
        @test grid.departure_times === departure_times
        @test grid.arrival_times === arrival_times
    end
end
