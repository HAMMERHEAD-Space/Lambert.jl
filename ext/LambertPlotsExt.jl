module LambertPlotsExt

using Lambert
using Plots
using LinearAlgebra
using SciMLBase
using StaticArraysCore

import Lambert: porkchop_plot

# Helper function to get colorbar title based on plot quantity
function _get_colorbar_title(plot_quantity::Symbol)
    if plot_quantity == :total_dv
        return "\nTotal ΔV [km/s]"
    elseif plot_quantity == :total_excess_velocity
        return "\nTotal V∞ [km/s]"
    elseif plot_quantity == :dv_departure || plot_quantity == :excess_velocity_departure
        return "\nDeparture ΔV [km/s]"
    elseif plot_quantity == :dv_arrival || plot_quantity == :excess_velocity_arrival
        return "\nArrival ΔV [km/s]"
    else
        return "\nΔV [km/s]"
    end
end
"""
    porkchop_plot(μ, state1_func, state2_func, departure_times, arrival_times;
                  solver=Lambert.IzzoSolver(), kwargs...)

Create a porkchop plot showing ΔV contours for different departure and arrival times.

# Arguments
- `μ`: Gravitational parameter of the central body [km³/s²]
- `state1_func`: Function that takes time and returns [r; v] (6-element vector) for departure body
- `state2_func`: Function that takes time and returns [r; v] (6-element vector) for arrival body
- `departure_times`: Vector of departure times to evaluate (in any units)
- `arrival_times`: Vector of arrival times to evaluate (in any units matching departure_times)

# Notes
- State functions should return 6-element vectors: [rx, ry, rz, vx, vy, vz]
- Time-of-flight is computed as `arrival_time - departure_time` and passed to LambertProblem
- Axes will display times in whatever units you provide

# Keyword Arguments
- `solver`: Lambert solver algorithm to use (default: IzzoSolver())
- `ensemble_alg`: Ensemble algorithm for parallel execution (default: EnsembleSerial()).
   Options: EnsembleSerial(), EnsembleThreads(), EnsembleDistributed()
- `max_deltav`: Maximum ΔV to display [km/s] (default: Inf)
- `levels`: Contour levels to plot (default: automatically determined)
- `tof_contours`: Whether to overlay time-of-flight contour lines (default: true)
- `tof_levels`: Specific TOF contour levels, or nothing for automatic (default: nothing)
- `tof_spacing`: Spacing between TOF lines when auto-generating, or nothing for automatic (default: nothing)
- `time_scale`: Scale factor to divide axis values by for display (default: 1.0).
   For example, if times are in seconds but you want axes in days, use 86400.0.
   This ONLY affects axis display, not the calculations.
- `plot_quantity`: Quantity to plot as contours (default: :total_dv).
   Options: :total_dv, :total_excess_velocity, :dv_departure, :dv_arrival, :excess_velocity_departure, :excess_velocity_arrival
   Note: :dv_* and :excess_velocity_* are aliases for the same quantity.
   Can also pass a vector of symbols (e.g., [:total_dv, :total_excess_velocity, :dv_departure]) to generate
   vertical subplots from a single Lambert solution run
- `title`: Plot title (default: "Porkchop Plot")
- `xlabel`: X-axis label (default: "Departure Time")
- `ylabel`: Y-axis label (default: "Arrival Time")
- `color`: Colormap (default: :turbo). Good options include :turbo, :jet, :plasma
- Additional keyword arguments are passed to Plots.contourf

# Returns
- A Plots.jl plot object
"""
function Lambert.porkchop_plot(
    μ::Number,
    state1_func::Function,
    state2_func::Function,
    departure_times::AbstractVector,
    arrival_times::AbstractVector;
    solver = Lambert.IzzoSolver(),
    ensemble_alg = SciMLBase.EnsembleSerial(),
    max_deltav::Number = Inf,
    levels::Union{Nothing,AbstractVector} = nothing,
    tof_contours::Bool = true,
    tof_levels::Union{Nothing,AbstractVector} = nothing,
    tof_spacing::Union{Nothing,Number} = nothing,
    time_scale::Number = 1.0,
    plot_quantity::Union{Symbol,AbstractVector{Symbol}} = :total_dv,
    title::String = "Porkchop Plot",
    xlabel::String = "Departure Time",
    ylabel::String = "Arrival Time",
    color = :turbo,
    kwargs...,
)
    n_dep = length(departure_times)
    n_arr = length(arrival_times)

    # Normalize plot_quantity to a vector
    plot_quantities = plot_quantity isa Symbol ? [plot_quantity] : plot_quantity
    n_quantities = length(plot_quantities)

    # Preallocate matrices for each quantity
    result_matrices = [fill(NaN, n_arr, n_dep) for _ = 1:n_quantities]

    # Build a list of all problem parameters
    problem_params = []
    indices = []

    for (i, t_dep) in enumerate(departure_times)
        state1 = state1_func(t_dep)
        r1 = SVector{3}(state1[1], state1[2], state1[3])
        v1 = SVector{3}(state1[4], state1[5], state1[6])

        for (j, t_arr) in enumerate(arrival_times)
            # Time of flight must be positive
            tof = t_arr - t_dep
            if tof <= 0
                continue
            end

            state2 = state2_func(t_arr)
            r2 = SVector{3}(state2[1], state2[2], state2[3])
            v2 = SVector{3}(state2[4], state2[5], state2[6])

            push!(problem_params, (r1, r2, tof, v1, v2))
            push!(indices, (j, i))
        end
    end

    # Create a base problem (will be remade for each trajectory)
    if isempty(problem_params)
        return nothing
    end

    r1_base, r2_base, tof_base, _, _ = problem_params[1]
    base_prob = Lambert.LambertProblem(μ, r1_base, r2_base, tof_base)

    # prob_func remakes the problem for each trajectory
    function prob_func(prob, i, repeat)
        r1, r2, tof, _, _ = problem_params[i]
        remake(prob; r1 = r1, r2 = r2, tof = tof)
    end

    # output_func computes all requested quantities from the solution
    function output_func(sol, i)
        _, _, _, v1, v2 = problem_params[i]

        if sol.retcode == :SUCCESS
            # Calculate ΔV at departure and arrival
            Δv1 = norm(sol.v1 .- v1)
            Δv2 = norm(sol.v2 .- v2)
            total_deltav = Δv1 + Δv2

            # Apply maximum ΔV filter
            if total_deltav <= max_deltav
                # Return all three quantities as a tuple
                return ((Δv1, Δv2, total_deltav), false)
            end
        end

        return ((NaN, NaN, NaN), false)
    end

    # Create and solve ensemble problem
    ensemble_prob = SciMLBase.EnsembleProblem(
        base_prob,
        prob_func = prob_func,
        output_func = output_func,
    )

    sim = SciMLBase.solve(
        ensemble_prob,
        solver,
        ensemble_alg;
        trajectories = length(problem_params),
    )

    # Fill in the result matrices
    for (idx, (j, i)) in enumerate(indices)
        Δv1, Δv2, total_dv = sim.u[idx]

        # Fill matrices based on which quantities were requested
        for (q_idx, qty) in enumerate(plot_quantities)
            result_matrices[q_idx][j, i] = if qty == :total_dv
                total_dv
            elseif qty == :total_excess_velocity
                total_dv  # Same as total_dv (sum of departure and arrival excess velocities)
            elseif qty == :dv_departure || qty == :excess_velocity_departure
                Δv1
            elseif qty == :dv_arrival || qty == :excess_velocity_arrival
                Δv2
            else
                error(
                    "Invalid plot_quantity: $qty. Must be :total_dv, :total_excess_velocity, :dv_departure, :dv_arrival, :excess_velocity_departure, or :excess_velocity_arrival",
                )
            end
        end
    end

    # Scale axes for display (does not affect calculations)
    departure_times_scaled = departure_times ./ time_scale
    arrival_times_scaled = arrival_times ./ time_scale

    # Create plots for each quantity
    plots = []
    for (q_idx, qty) in enumerate(plot_quantities)
        deltav_matrix = result_matrices[q_idx]

        # Determine contour levels for this quantity
        qty_levels = if levels === nothing
            # Automatically determine contour levels
            valid_dv = filter(!isnan, deltav_matrix)
            if !isempty(valid_dv)
                min_dv = minimum(valid_dv)
                max_dv = min(maximum(valid_dv), max_deltav)
                range(min_dv, max_dv, length = 20)
            else
                0:0.5:10
            end
        else
            levels
        end

        # Create the plot
        p = contourf(
            departure_times_scaled,
            arrival_times_scaled,
            deltav_matrix,
            levels = qty_levels,
            xlabel = xlabel,
            ylabel = ylabel,
            title = title,
            color = color,
            clims = (minimum(qty_levels), maximum(qty_levels)),
            colorbar_title = _get_colorbar_title(qty),
            size = (950, 700),
            titlefontsize = 20,
            guidefontsize = 16,
            tickfontsize = 10,
            colorbar_tickfontsize = 11,
            colorbar_titlefontsize = 14,
            left_margin = 8Plots.mm,
            right_margin = 8Plots.mm,
            top_margin = 5Plots.mm,
            bottom_margin = 5Plots.mm;
            kwargs...,
        )

        push!(plots, p)
    end

    # Add time-of-flight lines to all plots
    if tof_contours
        for p in plots
            _add_tof_lines!(
                p,
                departure_times_scaled,
                arrival_times_scaled,
                tof_levels,
                tof_spacing,
                time_scale,
            )
        end
    end

    # Return single plot or vertical subplots
    if length(plots) == 1
        return plots[1]
    else
        # Combine plots into vertical subplots with better spacing
        return plot(
            plots...,
            layout = (length(plots), 1),
            size = (1150, 700 * length(plots)),
            left_margin = 20Plots.mm,
            right_margin = 18Plots.mm,
            top_margin = 8Plots.mm,
            bottom_margin = 8Plots.mm,
        )
    end
end

# Helper function to add TOF lines to a plot
function _add_tof_lines!(
    p,
    departure_times_scaled,
    arrival_times_scaled,
    tof_levels,
    tof_spacing,
    time_scale,
)
    # Determine TOF levels
    min_dep = minimum(departure_times_scaled)
    max_dep = maximum(departure_times_scaled)
    min_arr = minimum(arrival_times_scaled)
    max_arr = maximum(arrival_times_scaled)

    if tof_levels === nothing
        # TOF range that covers the plot area
        min_tof = max(0, min_arr - max_dep)  # Minimum TOF (must be positive)
        max_tof = max_arr - min_dep  # Maximum TOF (latest arr - earliest dep)

        if tof_spacing === nothing
            # Auto-generate 5 evenly spaced lines
            tof_levels_scaled = range(min_tof, max_tof, length = 5)
        else
            # Use specified spacing
            spacing_scaled = tof_spacing / time_scale
            tof_levels_scaled = collect(min_tof:spacing_scaled:max_tof)
        end
    else
        tof_levels_scaled = filter(x -> x >= 0, tof_levels ./ time_scale)  # Filter negatives
    end

    # Calculate text rotation based on plot dimensions and data ranges
    # A line with slope=1 in data space appears at angle = atan((dy_pixels/dy_data) / (dx_pixels/dx_data))
    # Approximating: plot area is roughly 950*0.8 = 760 wide, 700*0.8 = 560 tall (accounting for margins/colorbar)
    plot_width_pixels = 760
    plot_height_pixels = 560
    data_width = max_dep - min_dep
    data_height = max_arr - min_arr

    # For a line with slope=1 in data: dy_data = dx_data
    # On screen: dy_pixels = dx_data * (plot_height_pixels / data_height)
    #            dx_pixels = dx_data * (plot_width_pixels / data_width)
    # Screen angle = atan(dy_pixels / dx_pixels)
    text_rotation =
        atand((plot_height_pixels / data_height) / (plot_width_pixels / data_width))

    # Draw diagonal TOF lines (arrival = departure + TOF)
    for tof in tof_levels_scaled
        # Line MUST have slope = 1 in data coordinates (arrival = departure + tof)
        dep_start = min_dep
        dep_end = max_dep
        arr_start = dep_start + tof
        arr_end = dep_end + tof

        # Clip to plot bounds so axes don't expand
        if arr_start < min_arr
            dep_start = min_arr - tof
            arr_start = min_arr
        end
        if arr_end > max_arr
            dep_end = max_arr - tof
            arr_end = max_arr
        end

        # Skip if line is completely outside bounds
        if dep_start > max_dep || dep_end < min_dep
            continue
        end

        # Calculate midpoint for label
        mid_dep = (dep_start + dep_end) / 2
        mid_arr = mid_dep + tof

        # Create label text first to determine gap size
        label_text = string(round(Int, tof))

        # Gap size based on text length (approximate character width in data units)
        # Rough estimate: each character is ~10 pixels at font size 14
        # Convert to data units considering the screen scaling
        char_width_pixels = 10
        num_chars = length(label_text)
        gap_data_units = (num_chars * char_width_pixels / plot_width_pixels) * data_width
        gap = gap_data_units * 0.75  # Add margin/padding around text

        # Draw line in two segments with gap in middle - thicker lines
        plot!(
            p,
            [dep_start, mid_dep - gap],
            [arr_start, mid_arr - gap],
            color = :black,
            linewidth = 4.5,
            linealpha = 1.0,
            label = "",
            linestyle = :dash,
            xlims = (min_dep, max_dep),
            ylims = (min_arr, max_arr),
        )

        plot!(
            p,
            [mid_dep + gap, dep_end],
            [mid_arr + gap, arr_end],
            color = :black,
            linewidth = 4.5,
            linealpha = 1.0,
            label = "",
            linestyle = :dash,
            xlims = (min_dep, max_dep),
            ylims = (min_arr, max_arr),
        )

        # Note: Plots.jl only supports :bold, no heavier weights
        # Using smaller size with bold for heavier appearance
        annotate!(
            p,
            mid_dep,
            mid_arr,
            text(label_text, 16, :center, :black, :bold, rotation = text_rotation),
        )
    end
end

end # module
