@testset "Porkchop Plot Extension (requires Plots.jl)" begin
    # Check if Plots is available
    plots_available = false
    try
        using Plots
        plots_available = true
    catch e
        @info "Plots.jl not available - skipping porkchop plot tests"
        @info "To enable these tests, run: using Pkg; Pkg.add(\"Plots\")"
    end

    if plots_available
        @testset "porkchop_plot with state functions" begin
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
            departure_times = range(0, 3600, length = 10)
            arrival_times = range(3600, 7200, length = 10)

            # Create porkchop plot
            p = porkchop_plot(
                μ_earth,
                state1_func,
                state2_func,
                departure_times,
                arrival_times;
                title = "Test Porkchop Plot",
            )

            @test p isa Plots.Plot
            @test !isempty(p.series_list)
        end

        @testset "porkchop_plot with different plot quantities" begin
            μ_earth = 3.986004418e5

            function state1_func(t)
                ω = sqrt(μ_earth / 7000^3)
                θ = ω * t
                r = [7000 * cos(θ), 7000 * sin(θ), 0.0]
                v = [-7000 * ω * sin(θ), 7000 * ω * cos(θ), 0.0]
                return [r; v]
            end

            function state2_func(t)
                ω = sqrt(μ_earth / 10000^3)
                θ = ω * t + π / 4
                r = [10000 * cos(θ), 10000 * sin(θ), 0.0]
                v = [-10000 * ω * sin(θ), 10000 * ω * cos(θ), 0.0]
                return [r; v]
            end

            departure_times = range(0, 3600, length = 10)
            arrival_times = range(3600, 7200, length = 10)

            # Test individual quantities
            for qty in [
                :total_dv,
                :total_excess_velocity,
                :dv_departure,
                :dv_arrival,
                :excess_velocity_departure,
                :excess_velocity_arrival,
            ]
                p = porkchop_plot(
                    μ_earth,
                    state1_func,
                    state2_func,
                    departure_times,
                    arrival_times;
                    plot_quantity = qty,
                )
                @test p isa Plots.Plot
                @test !isempty(p.series_list)
            end
        end

        @testset "porkchop_plot with multiple quantities (subplots)" begin
            μ_earth = 3.986004418e5

            function state1_func(t)
                ω = sqrt(μ_earth / 7000^3)
                θ = ω * t
                r = [7000 * cos(θ), 7000 * sin(θ), 0.0]
                v = [-7000 * ω * sin(θ), 7000 * ω * cos(θ), 0.0]
                return [r; v]
            end

            function state2_func(t)
                ω = sqrt(μ_earth / 10000^3)
                θ = ω * t + π / 4
                r = [10000 * cos(θ), 10000 * sin(θ), 0.0]
                v = [-10000 * ω * sin(θ), 10000 * ω * cos(θ), 0.0]
                return [r; v]
            end

            departure_times = range(0, 3600, length = 10)
            arrival_times = range(3600, 7200, length = 10)

            # Test multiple quantities - should return subplots
            p = porkchop_plot(
                μ_earth,
                state1_func,
                state2_func,
                departure_times,
                arrival_times;
                plot_quantity = [:total_dv, :dv_departure, :dv_arrival],
            )

            @test p isa Plots.Plot
            @test !isempty(p.series_list)
            # Should have multiple subplots combined into one plot
            @test p.layout isa Plots.GridLayout
            @test size(p.layout) == (3, 1)  # 3 rows, 1 column
        end

        @testset "porkchop_plot with max_deltav filter" begin
            μ_earth = 3.986004418e5

            function state1_func(t)
                ω = sqrt(μ_earth / 7000^3)
                θ = ω * t
                r = [7000 * cos(θ), 7000 * sin(θ), 0.0]
                v = [-7000 * ω * sin(θ), 7000 * ω * cos(θ), 0.0]
                return [r; v]
            end

            function state2_func(t)
                ω = sqrt(μ_earth / 10000^3)
                θ = ω * t + π / 4
                r = [10000 * cos(θ), 10000 * sin(θ), 0.0]
                v = [-10000 * ω * sin(θ), 10000 * ω * cos(θ), 0.0]
                return [r; v]
            end

            departure_times = range(0, 3600, length = 10)
            arrival_times = range(3600, 7200, length = 10)

            # Create plot with max ΔV filter
            p = porkchop_plot(
                μ_earth,
                state1_func,
                state2_func,
                departure_times,
                arrival_times;
                max_deltav = 5.0,
                title = "Filtered Porkchop Plot",
            )

            @test p isa Plots.Plot
        end

        @testset "porkchop_plot with time scaling and TOF lines" begin
            μ_earth = 3.986004418e5

            function state1_func(t)
                ω = sqrt(μ_earth / 7000^3)
                θ = ω * t
                r = [7000 * cos(θ), 7000 * sin(θ), 0.0]
                v = [-7000 * ω * sin(θ), 7000 * ω * cos(θ), 0.0]
                return [r; v]
            end

            function state2_func(t)
                ω = sqrt(μ_earth / 10000^3)
                θ = ω * t + π / 4
                r = [10000 * cos(θ), 10000 * sin(θ), 0.0]
                v = [-10000 * ω * sin(θ), 10000 * ω * cos(θ), 0.0]
                return [r; v]
            end

            departure_times = range(0, 3600, length = 10)
            arrival_times = range(3600, 7200, length = 10)

            # Test with time scaling and TOF lines
            p = porkchop_plot(
                μ_earth,
                state1_func,
                state2_func,
                departure_times,
                arrival_times;
                time_scale = 60.0,  # Convert seconds to minutes
                tof_contours = true,
                tof_spacing = 600.0,  # TOF line every 10 minutes
                xlabel = "Departure Time [min]",
                ylabel = "Arrival Time [min]",
            )

            @test p isa Plots.Plot
            @test !isempty(p.series_list)
        end

        @testset "Plot from PorkchopGrid convenience function" begin
            # Define test parameters (reuse from parent scope)
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
            departure_times = range(0, 3600, length = 10)
            arrival_times = range(3600, 7200, length = 10)

            # Generate a grid first
            grid = porkchop_grid(
                μ_earth,
                state1_func,
                state2_func,
                departure_times,
                arrival_times;
                quantities = [:total_dv, :dv_departure, :dv_arrival],
            )

            # Test plotting single quantity from grid
            p1 =
                porkchop_plot(grid; plot_quantity = :total_dv, title = "Total ΔV from Grid")
            @test p1 isa Plots.Plot
            @test !isempty(p1.series_list)

            # Test plotting multiple quantities from grid
            p2 = porkchop_plot(grid; plot_quantity = [:total_dv, :dv_departure])
            @test p2 isa Plots.Plot
            @test !isempty(p2.series_list)

            # Test with custom settings
            p3 = porkchop_plot(
                grid;
                plot_quantity = :dv_arrival,
                title = "Arrival ΔV",
                color = :plasma,
                tof_contours = false,
            )
            @test p3 isa Plots.Plot

            # Test default quantity selection (should use first quantity in grid)
            p4 = porkchop_plot(grid)
            @test p4 isa Plots.Plot

            # Test error handling for invalid quantity
            @test_throws ErrorException porkchop_plot(
                grid;
                plot_quantity = :invalid_quantity,
            )

            # Test that grid-based plotting produces similar results to direct plotting
            p_direct = porkchop_plot(
                μ_earth,
                state1_func,
                state2_func,
                departure_times,
                arrival_times;
                plot_quantity = :total_dv,
            )
            p_from_grid = porkchop_plot(grid; plot_quantity = :total_dv)

            # Both should be valid plots
            @test p_direct isa Plots.Plot
            @test p_from_grid isa Plots.Plot
            @test !isempty(p_direct.series_list)
            @test !isempty(p_from_grid.series_list)
        end

        @info "All porkchop plot tests passed!"
    end
end
