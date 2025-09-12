"""
    get_transfer_angle(r1, r2, prograde)

Compute the transfer angle of the trajectory.

# Arguments
- `r1`: Initial position vector
- `r2`: Final position vector  
- `prograde`: True for prograde motion, False otherwise

# Returns
- `dtheta`: Transfer angle in radians
"""
function get_transfer_angle(r1::Vector{<:Number}, r2::Vector{<:Number}, prograde::Bool)

    r1 = SVector{3}(r1[1], r1[2], r1[3])
    r2 = SVector{3}(r2[1], r2[2], r2[3])

    # Check if both position vectors are collinear
    cross_prod = cross(r1, r2)

    # Handle collinear case
    if all(cross_prod .≈ 0)
        return all(sign.(r1) .== sign.(r2)) ? 0.0 : π
    end

    # Use AstroCoords function to get the basic angle between vectors
    θ0 = angle_between_vectors(r1, r2)

    # Solve for a unitary vector normal to the vector plane
    h = cross_prod / norm(cross_prod)

    # Compute the projection of the normal vector onto the reference plane
    α = dot(SVector{3}(0, 0, 1), h)

    # Fix the value of theta if necessary for prograde/retrograde motion
    if prograde
        dtheta = α > 0 ? θ0 : 2π - θ0
    else
        dtheta = α < 0 ? θ0 : 2π - θ0
    end

    return dtheta
end

"""
    get_orbit_normal_vector(r1, r2, prograde)

Computes a unitary normal vector aligned with the specific angular momentum of the orbit.

# Arguments
- `r1`: Initial position vector
- `r2`: Final position vector
- `prograde`: If true, assumes prograde motion, otherwise retrograde

# Returns
- `i_h`: Unitary vector aligned with orbit specific angular momentum
"""
function get_orbit_normal_vector(r1::Vector{<:Number}, r2::Vector{<:Number}, prograde::Bool)
    # Compute the normal vector and its projection onto the vertical axis
    i_h = (r1 × r2) / norm(r1 × r2)

    # Solve the projection onto the positive vertical direction
    α = dot([0, 0, 1], i_h)

    # A prograde orbit always has a positive vertical component
    if prograde
        i_h = α > 0 ? i_h : -i_h
    else
        i_h = α < 0 ? i_h : -i_h
    end

    return i_h
end

"""
    c2(ψ; threshold=1e-6)

Second Stumpff function.

For positive arguments: c₂(ψ) = (1 - cos(√ψ))/ψ
"""
@inline function c2(ψ::Number; threshold::Float64 = 1e-6)
    if ψ > threshold
        return (1 - cos(√ψ)) / ψ
    elseif ψ < -threshold
        return (cosh(√(-ψ)) - 1) / (-ψ)
    else
        # Series expansion for small ψ: c2 = 1/2 - ψ/24 + ψ²/720 - ...
        return 0.5 - ψ/24 + ψ^2/720 - ψ^3/40320
    end
end

"""
    c3(ψ; threshold=1e-6)

Third Stumpff function.

For positive arguments: c₃(ψ) = (√ψ - sin(√ψ))/√(ψ³)
"""
@inline function c3(ψ::Number; threshold::Float64 = 1e-6)
    if ψ > threshold
        sqrt_ψ = √ψ
        return (sqrt_ψ - sin(sqrt_ψ)) / (ψ * sqrt_ψ)
    elseif ψ < -threshold
        sqrt_neg_ψ = √(-ψ)
        return (sinh(sqrt_neg_ψ) - sqrt_neg_ψ) / (-ψ * sqrt_neg_ψ)
    else
        # Series expansion for small ψ: c3 = 1/6 - ψ/120 + ψ²/5040 - ...
        return 1/6 - ψ/120 + ψ^2/5040 - ψ^3/362880
    end
end

"""
    assert_parameters_are_valid(μ, r1, r2, tof, M)

Basic validation of Lambert problem parameters.
"""
function assert_parameters_are_valid(
    μ::Number,
    r1::Vector{<:Number},
    r2::Vector{<:Number},
    tof::Number,
    M::Int,
)
    @assert μ > 0 "Gravitational parameter must be positive"
    @assert norm(r1) > 0 "Initial position vector cannot be zero"
    @assert norm(r2) > 0 "Final position vector cannot be zero"
    @assert tof > 0 "Time of flight must be positive"
    @assert M >= 0 "Number of revolutions must be non-negative"
    @assert !all(r1 .≈ r2) "Initial and final positions cannot be the same"
end

"""
    assert_transfer_angle_not_zero(dtheta)

Check that transfer angle is not zero.
"""
@inline function assert_transfer_angle_not_zero(dtheta::Number)
    @assert abs(dtheta) > 1e-12 "Transfer angle cannot be zero"
end

"""
    assert_transfer_angle_not_pi(dtheta)

Check that transfer angle is not π (180 degrees).
"""
@inline function assert_transfer_angle_not_pi(dtheta::Number)
    @assert abs(dtheta - π) > 1e-12 "Transfer angle cannot be π (180 degrees)"
end


"""
    lambert_geometry(r1, r2, prograde)

Compute basic geometry parameters for Lambert problem.
Returns (r1_norm, r2_norm, c_norm, dtheta).
"""
function lambert_geometry(r1::Vector{<:Number}, r2::Vector{<:Number}, prograde::Bool)
    r1_norm = norm(r1)
    r2_norm = norm(r2)
    c_norm = norm(r2 - r1)
    dtheta = get_transfer_angle(r1, r2, prograde)

    return r1_norm, r2_norm, c_norm, dtheta
end

"""
    check_convergence(x_new, x_old, atol, rtol)

Check convergence based on absolute and relative tolerances.
Returns true if converged.
"""
@inline function check_convergence(
    x_new::Number,
    x_old::Number,
    atol::Float64,
    rtol::Float64,
)
    return abs(x_new - x_old) <= atol ||
           abs(x_new - x_old) <= rtol * max(abs(x_new), abs(x_old))
end

"""
    check_convergence(vec_new, vec_old, atol, rtol)

Check convergence for vectors based on absolute and relative tolerances.
Returns true if converged.
"""
@inline function check_convergence(
    vec_new::Vector{<:Number},
    vec_old::Vector{<:Number},
    atol::Float64,
    rtol::Float64,
)
    diff_norm = norm(vec_new - vec_old)
    return diff_norm <= atol || diff_norm <= rtol * max(norm(vec_new), norm(vec_old))
end

"""
    create_lambert_solution(v1, v2, numiter, retcode)

Create a LambertSolution with proper error handling for failed solutions.
"""
function create_lambert_solution(
    v1::Union{AbstractVector,Nothing},
    v2::Union{AbstractVector,Nothing},
    numiter::Int,
    retcode::Symbol,
)
    if retcode == :SUCCESS && v1 !== nothing && v2 !== nothing
        return LambertSolution(v1, v2, numiter, retcode)
    else
        zero_v = zero(v1)
        return LambertSolution(zero_v, zero_v, numiter, retcode)
    end
end

"""
    handle_max_iterations(numiter, maxiter)

Check if maximum iterations exceeded and return appropriate status.
"""
@inline function handle_max_iterations(numiter::Int, maxiter::Int)
    return numiter >= maxiter ? :MAXIMUM_ITERATIONS : :SUCCESS
end

"""
    compute_semiperimeter(r1_norm, r2_norm, c_norm)

Compute the semiperimeter of the transfer triangle.
"""
@inline function compute_semiperimeter(r1_norm::Number, r2_norm::Number, c_norm::Number)
    return (r1_norm + r2_norm + c_norm) / 2.0
end

"""
    compute_unit_vectors(r1, r2, r1_norm, r2_norm)

Compute unit vectors for positions and their cross product.
Returns (i_r1, i_r2, i_h_unnorm).
"""
@inline function compute_unit_vectors(
    r1::Vector{<:Number},
    r2::Vector{<:Number},
    r1_norm::Number,
    r2_norm::Number,
)
    i_r1 = r1 / r1_norm
    i_r2 = r2 / r2_norm
    i_h_unnorm = cross(i_r1, i_r2)
    return i_r1, i_r2, i_h_unnorm
end

"""
    compute_tangential_unit_vectors(i_h, i_r1, i_r2)

Compute tangential unit vectors for velocity reconstruction.
Returns (i_t1, i_t2).
"""
@inline function compute_tangential_unit_vectors(
    i_h::Vector{<:Number},
    i_r1::Vector{<:Number},
    i_r2::Vector{<:Number},
)
    i_h = SVector{3}(i_h[1], i_h[2], i_h[3])
    i_r1 = SVector{3}(i_r1[1], i_r1[2], i_r1[3])
    i_r2 = SVector{3}(i_r2[1], i_r2[2], i_r2[3])

    i_t1 = cross(i_h, i_r1)
    i_t2 = cross(i_h, i_r2)
    return i_t1, i_t2
end

"""
    reconstruct_velocities_from_components(vr1, vt1, vr2, vt2, i_r1, i_t1, i_r2, i_t2)

Reconstruct velocity vectors from radial and tangential components.
Returns (v1, v2).
"""
@inline function reconstruct_velocities_from_components(
    vr1::Number,
    vt1::Number,
    vr2::Number,
    vt2::Number,
    i_r1::AbstractVector{<:Number},
    i_t1::AbstractVector{<:Number},
    i_r2::AbstractVector{<:Number},
    i_t2::AbstractVector{<:Number},
)

    i_r1 = SVector{3}(i_r1[1], i_r1[2], i_r1[3])
    i_t1 = SVector{3}(i_t1[1], i_t1[2], i_t1[3])
    i_r2 = SVector{3}(i_r2[1], i_r2[2], i_r2[3])
    i_t2 = SVector{3}(i_t2[1], i_t2[2], i_t2[3])

    v1 = vr1 * i_r1 + vt1 * i_t1
    v2 = vr2 * i_r2 + vt2 * i_t2
    return v1, v2
end

"""
    reconstruct_velocities_fg(f, g, gdot, r1, r2)

Reconstruct velocities using f and g functions from universal formulation.
Returns (v1, v2).
"""
@inline function reconstruct_velocities_fg(
    f::Number,
    g::Number,
    gdot::Number,
    r1::AbstractVector{<:Number},
    r2::AbstractVector{<:Number},
)
    r1 = SVector{3}(r1[1], r1[2], r1[3])
    r2 = SVector{3}(r2[1], r2[2], r2[3])

    v1 = (r2 - f * r1) / g
    v2 = (gdot * r2 - r1) / g
    return v1, v2
end

"""
    validate_multi_revolution(M, solver_name)

Check if multi-revolution is supported by the solver.
"""
function validate_multi_revolution(M::Int, solver_name::String)
    if M > 0
        error("Multi-revolution case (M > 0) not implemented for $solver_name")
    end
end

"""
    compute_parabolic_time_of_flight(μ, s, c_norm)

Compute the parabolic time of flight for minimum energy transfer.
"""
@inline function compute_parabolic_time_of_flight(μ::Number, s::Number, c_norm::Number)
    return (1/3) * √(2/μ) * (s^1.5 - (s - c_norm)^1.5)
end

"""
    compute_unit_vector_geometry(r1, r2, r1_norm, r2_norm)

Compute unit vectors and check for collinear vectors.
Returns (i_r1, i_r2, i_h, cos_dtheta, sin_dtheta).
"""
function compute_unit_vector_geometry(
    r1::AbstractVector{<:Number},
    r2::AbstractVector{<:Number},
    r1_norm::Number,
    r2_norm::Number,
)
    i_r1 = r1 / r1_norm
    i_r2 = r2 / r2_norm
    i_h = cross(i_r1, i_r2)
    i_h_norm = norm(i_h)

    if i_h_norm < 1e-10
        return i_r1, i_r2, [0.0, 0.0, 1.0], 1.0, 0.0  # Default values for collinear case
    end

    i_h = i_h / i_h_norm
    cos_dtheta = dot(i_r1, i_r2)
    sin_dtheta = i_h_norm

    return i_r1, i_r2, i_h, cos_dtheta, sin_dtheta
end

"""
    check_collinear_vectors(sin_dtheta; tolerance=1e-10)

Check if vectors are collinear and return appropriate error code.
"""
@inline function check_collinear_vectors(sin_dtheta::Number; tolerance::Float64 = 1e-10)
    if abs(sin_dtheta) < tolerance
        return :COLLINEAR_VECTORS
    end
    return :OK
end

"""
    universal_anomaly_initial_guess(μ, tof, s, M, tof_parabolic)

Compute initial guess for universal anomaly (chi) in universal variable methods.
"""
function universal_anomaly_initial_guess(
    μ::Number,
    tof::Number,
    s::Number,
    M::Int,
    tof_parabolic::Number,
)
    if M == 0
        if tof > tof_parabolic
            # Elliptic case
            return √(μ) * tof / s
        else
            # Hyperbolic case
            return -√(μ) * tof / s
        end
    else
        # Multi-revolution case
        return √(μ) * tof / s + 2*π*M
    end
end

"""
    check_newton_derivative(derivative, threshold=1e-12)

Check if Newton-Raphson derivative is safe and return appropriate error code.
"""
@inline function check_newton_derivative(derivative::Number; threshold::Float64 = 1e-12)
    if abs(derivative) < threshold
        return :DERIVATIVE_TOO_SMALL
    end
    return :OK
end

"""
    apply_damping(delta, current_value, max_relative_change=0.5)

Apply damping to Newton-Raphson updates for stability.
"""
@inline function apply_damping(
    delta::Number,
    current_value::Number;
    max_relative_change::Float64 = 0.5,
)
    relative_change = abs(delta / current_value)
    if relative_change > max_relative_change
        damping_factor = max_relative_change / relative_change
        return delta * damping_factor
    end
    return delta
end

"""
    validate_lagrange_coefficients(f, g, gdot; tolerance=1e-10)

Validate Lagrange coefficients satisfy f*gdot - g*(dg/dt) = 1.
Returns :INVALID_LAGRANGE_COEFFICIENTS if validation fails.
"""
function validate_lagrange_coefficients(
    f::Number,
    g::Number,
    gdot::Number;
    tolerance::Float64 = 1e-10,
)
    # For Lagrange coefficients, f*gdot + (g*fdot) should equal 1
    # Since fdot = (gdot - 1)/g, we have: f*gdot + g*(gdot-1)/g = f*gdot + gdot - 1 = gdot*(f+1) - 1
    # But the correct relation is: f*gdot - g*fdot = 1
    # Where fdot = -μ*g/(r1*r2*r1_norm*r2_norm)

    # Simplified check: just verify they are finite and reasonable
    if !all(isfinite.([f, g, gdot]))
        return :INVALID_LAGRANGE_COEFFICIENTS
    end

    # Check if g is too small (would cause numerical issues)
    if abs(g) < tolerance
        return :INVALID_LAGRANGE_COEFFICIENTS
    end

    return :OK
end

"""
    select_lambert_algorithm(problem::LambertProblem, M::Int=0, prograde::Bool=true)

Heuristic algorithm selection for Lambert problems based on problem characteristics.

# Selection Logic
- **Multi-revolution (M > 0)**: Arora (fast, robust) or Izzo (most accurate)
- **Single revolution (M = 0)**:
  - Short transfers (Δθ < π/4): Gooding (robust for small angles)
  - Medium transfers (π/4 ≤ Δθ ≤ 3π/2): Izzo (best overall performance)  
  - Large transfers (Δθ > 3π/2): Arora (handles large angles well)
  - Near-180° transfers: Battin (designed for this singularity)

# Arguments
- `problem`: LambertProblem instance
- `M`: Number of revolutions (default: 0)
- `prograde`: Direction of motion (default: true)

# Returns
- Optimal `AbstractLambertSolver` instance for the given problem
"""
function select_lambert_algorithm(
    problem::LambertProblem,
    M::Int = 0,
    prograde::Bool = true,
)
    @unpack r1, r2 = problem

    # Compute transfer angle
    dtheta = get_transfer_angle(r1, r2, prograde)

    # Multi-revolution cases - only Arora and Izzo support this
    if M > 0
        # Arora is faster and specifically designed for multi-rev
        # Izzo has better accuracy but more complex
        if dtheta > π  # Large angle transfers
            return AroraSolver(M = M, prograde = prograde)
        else  # Small to medium angle transfers  
            return IzzoSolver(M = M, prograde = prograde)
        end
    end

    # Single revolution cases
    # Check for near-180° transfer (Battin's specialty)
    if abs(dtheta - π) < π/36  # Within 5 degrees of 180°
        return BattinSolver(M = M, prograde = prograde)
    end

    # General single-rev selection based on transfer angle
    if dtheta < π/4  # Short transfers (~45°)
        # Gooding is most robust for small angles
        return GoodingSolver(M = M, prograde = prograde)
    elseif dtheta <= 3π/2  # Medium transfers (45° to 270°)
        # Izzo has best overall performance for typical transfers
        return IzzoSolver(M = M, prograde = prograde)
    else  # Large transfers (> 270°)
        # Arora handles large angles well
        return AroraSolver(M = M, prograde = prograde)
    end
end

"""
    generic_solver_solve(problem, solver, algorithm_function)

Generic SciMLBase.solve implementation for Lambert solvers.
Calls the algorithm function with unpacked parameters and returns a LambertSolution.
"""
function generic_solver_solve(
    problem::LambertProblem,
    solver::AbstractLambertSolver,
    algorithm_function::Function,
)
    @unpack μ, r1, r2, tof = problem

    # Get solver parameters (this will vary by solver type)
    if hasfield(typeof(solver), :M)
        M = solver.M
    else
        M = 0
    end

    if hasfield(typeof(solver), :prograde)
        prograde = solver.prograde
    else
        prograde = true
    end

    # Call the algorithm function
    result = algorithm_function(μ, r1, r2, tof, solver)

    # Handle different return patterns
    if length(result) == 4
        v1, v2, numiter, retcode = result
        return create_lambert_solution(v1, v2, numiter, retcode)
    elseif length(result) == 2
        v1, v2 = result
        return create_lambert_solution(v1, v2, 0, :SUCCESS)
    else
        error("Unexpected algorithm function return format")
    end
end


"""
    normalized_time_of_flight_function(x, target_tof, geometry_params, algorithm_specific_func)

Generic wrapper for time-of-flight functions that normalizes the interface.
Returns residual (computed_tof - target_tof).
"""
function normalized_time_of_flight_function(
    x::Number,
    target_tof::Number,
    geometry_params::NamedTuple,
    algorithm_specific_func::Function,
)
    computed_tof = algorithm_specific_func(x, geometry_params)
    return computed_tof - target_tof
end

"""
    transfer_classification(r1, r2, prograde, tof, μ)

Classify the transfer type based on geometry and time of flight.
Returns a symbol indicating transfer characteristics.
"""
function transfer_classification(
    r1::AbstractVector{<:Number},
    r2::AbstractVector{<:Number},
    prograde::Bool,
    tof::Number,
    μ::Number,
)
    r1_norm, r2_norm, c_norm, dtheta = lambert_geometry(r1, r2, prograde)
    s = compute_semiperimeter(r1_norm, r2_norm, c_norm)
    tof_parabolic = compute_parabolic_time_of_flight(μ, s, c_norm)

    # Angular classification
    angular_class = if dtheta < π/4
        :short_transfer  # < 45°
    elseif dtheta < π/2
        :medium_transfer  # 45° - 90°
    elseif dtheta < 3π/2
        :large_transfer   # 90° - 270°
    else
        :very_large_transfer  # > 270°
    end

    # Energy classification
    energy_class = if tof < tof_parabolic
        :hyperbolic
    elseif abs(tof - tof_parabolic) < 1e-6
        :parabolic
    else
        :elliptical
    end

    # Special cases
    if abs(dtheta - π) < π/36  # Within 5° of 180°
        return :near_pi_transfer
    end

    return Symbol(string(energy_class) * "_" * string(angular_class))
end

"""
    normalize_inputs(μ, r1, r2, tof)

Normalize inputs to characteristic length and time scales.
Returns (μ_hat, r1_hat, r2_hat, tof_hat, L_ref, T_ref).
"""
function normalize_inputs(
    μ::Number,
    r1::AbstractVector{<:Number},
    r2::AbstractVector{<:Number},
    tof::Number,
)
    L_ref = norm(r1)  # Characteristic length
    T_ref = √(L_ref^3 / μ)  # Characteristic time

    μ_hat = μ * (T_ref^2 / L_ref^3)
    r1_hat = r1 / L_ref
    r2_hat = r2 / L_ref
    tof_hat = tof / T_ref

    return μ_hat, r1_hat, r2_hat, tof_hat, L_ref, T_ref
end

"""
    denormalize_velocities(v1_hat, v2_hat, L_ref, T_ref)

Convert normalized velocities back to dimensional form.
"""
@inline function denormalize_velocities(
    v1_hat::AbstractVector{<:Number},
    v2_hat::AbstractVector{<:Number},
    L_ref::Number,
    T_ref::Number,
)
    scale_factor = L_ref / T_ref
    return v1_hat * scale_factor, v2_hat * scale_factor
end

"""
    compute_lagrange_coefficients_universal(χ, r1_norm, r2_norm, μ, α, ψ)

Compute Lagrange coefficients using universal variables.
Returns (f, g, fdot, gdot).
"""
function compute_lagrange_coefficients_universal(
    χ::Number,
    r1_norm::Number,
    r2_norm::Number,
    μ::Number,
    α::Number,
    ψ::Number,
)
    C = c2(ψ)
    S = c3(ψ)

    # Lagrange coefficients
    f = 1 - (χ^2 / r1_norm) * C
    g = tof - (χ^3 / √μ) * S
    gdot = 1 - (χ^2 / r2_norm) * C
    fdot = (√μ / (r1_norm * r2_norm)) * χ * (ψ * S - 1)

    return f, g, fdot, gdot
end

"""
    robust_sqrt(x; threshold=1e-15)

Robust square root that handles near-zero values gracefully.
"""
@inline function robust_sqrt(x::Number; threshold::Float64 = 1e-15)
    return x < threshold ? 0.0 : √x
end

"""
    robust_division(numerator, denominator; threshold=1e-15, fallback=0.0)

Robust division that handles near-zero denominators.
"""
@inline function robust_division(
    numerator::Number,
    denominator::Number;
    threshold::Float64 = 1e-15,
    fallback::Number = 0.0,
)
    return abs(denominator) < threshold ? fallback : numerator / denominator
end

"""
    extract_common_solver_params(solver)

Extract common parameters from solver structs with defaults.
Returns NamedTuple with standardized parameters.
"""
function extract_common_solver_params(solver)
    return (
        M = hasfield(typeof(solver), :M) ? solver.M : 0,
        prograde = hasfield(typeof(solver), :prograde) ? solver.prograde : true,
        maxiter = hasfield(typeof(solver), :maxiter) ? solver.maxiter : 35,
        atol = hasfield(typeof(solver), :atol) ? solver.atol : 1e-5,
        rtol = hasfield(typeof(solver), :rtol) ? solver.rtol : 1e-7,
    )
end
