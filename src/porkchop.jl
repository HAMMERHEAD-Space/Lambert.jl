export PorkchopGrid, porkchop_grid
"""
    PorkchopGrid

A struct containing the results of a porkchop plot grid calculation.

# Fields
- `departure_times`: Vector of departure times used in the calculation
- `arrival_times`: Vector of arrival times used in the calculation
- `quantities`: Vector of quantity symbols that were computed
- `data`: Dictionary mapping quantity symbols to their corresponding matrices

# Usage
Access data for a specific quantity using `grid.data[:total_dv]`, `grid.data[:dv_departure]`, etc.
"""
struct PorkchopGrid
    departure_times::AbstractVector
    arrival_times::AbstractVector
    quantities::Vector{Symbol}
    data::Dict{Symbol,Matrix{Float64}}
end

# Convenience methods for accessing data
Base.getindex(grid::PorkchopGrid, quantity::Symbol) = grid.data[quantity]
Base.haskey(grid::PorkchopGrid, quantity::Symbol) = haskey(grid.data, quantity)
Base.keys(grid::PorkchopGrid) = keys(grid.data)

"""
    porkchop_grid(μ, state1_func, state2_func, departure_times, arrival_times;
                  solver=IzzoSolver(), kwargs...)

Generate a porkchop grid of ΔV values for different departure and arrival times.

# Arguments
- `μ`: Gravitational parameter of the central body [km³/s²]
- `state1_func`: Function that takes time and returns [r; v] (6-element vector) for departure body
- `state2_func`: Function that takes time and returns [r; v] (6-element vector) for arrival body
- `departure_times`: Vector of departure times to evaluate (in any units)
- `arrival_times`: Vector of arrival times to evaluate (in any units matching departure_times)

# Keyword Arguments
- `solver`: Lambert solver algorithm to use (default: IzzoSolver())
- `ensemble_alg`: Ensemble algorithm for parallel execution (default: EnsembleSerial()).
  Options: EnsembleSerial(), EnsembleThreads(), EnsembleDistributed()
- `max_deltav`: Maximum ΔV to display [km/s] (default: Inf)
- `quantities`: Vector of quantities to compute (default: [:total_dv]).
  Options: :total_dv, :total_excess_velocity, :dv_departure, :dv_arrival, :excess_velocity_departure, :excess_velocity_arrival

# Returns
- `PorkchopGrid`: A struct containing departure_times, arrival_times, quantities, and data dictionary

# Notes
- State functions should return 6-element vectors: [rx, ry, rz, vx, vy, vz]
- Time-of-flight is computed as `arrival_time - departure_time` and passed to LambertProblem
- Result matrices have dimensions (length(arrival_times), length(departure_times))
- Invalid solutions (negative TOF, failed convergence, etc.) are filled with NaN
- Access data using `grid[:total_dv]` or `grid.data[:total_dv]`
"""
function porkchop_grid(
    μ::Number,
    state1_func::Function,
    state2_func::Function,
    departure_times::AbstractVector,
    arrival_times::AbstractVector;
    solver = IzzoSolver(),
    ensemble_alg = SciMLBase.EnsembleSerial(),
    max_deltav::Number = Inf,
    quantities::AbstractVector{Symbol} = [:total_dv],
)
    n_dep = length(departure_times)
    n_arr = length(arrival_times)
    n_quantities = length(quantities)

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
        # Create empty data dictionary for all requested quantities
        data_dict = Dict{Symbol,Matrix{Float64}}()
        for qty in quantities
            data_dict[qty] = result_matrices[1]  # All NaN matrix
        end
        return PorkchopGrid(departure_times, arrival_times, quantities, data_dict)
    end

    r1_base, r2_base, tof_base, _, _ = problem_params[1]
    base_prob = LambertProblem(μ, r1_base, r2_base, tof_base)

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
        for (q_idx, qty) in enumerate(quantities)
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
                    "Invalid quantity: $qty. Must be :total_dv, :total_excess_velocity, :dv_departure, :dv_arrival, :excess_velocity_departure, or :excess_velocity_arrival",
                )
            end
        end
    end

    # Create data dictionary mapping quantities to their matrices
    data_dict = Dict{Symbol,Matrix{Float64}}()
    for (q_idx, qty) in enumerate(quantities)
        data_dict[qty] = result_matrices[q_idx]
    end

    return PorkchopGrid(departure_times, arrival_times, quantities, data_dict)
end
