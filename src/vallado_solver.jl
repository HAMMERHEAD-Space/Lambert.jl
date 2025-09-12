export vallado2013, ValladoSolver

"""
    ValladoSolver

Vallado's algorithm uses the universal formulation with a bisection method.

This guarantees convergence but requires more iterations than other methods.

# Fields
- `M::Int`: Number of revolutions (default: 0)
- `prograde::Bool`: Direction of motion - true for prograde, false for retrograde (default: true)
- `maxiter::Int`: Maximum number of iterations (default: 100)
- `rtol::Float64`: Relative tolerance (default: 1e-7)

# References
Vallado, D. A. (2001). Fundamentals of astrodynamics and applications 
(2nd ed.). Space Technology Library. Springer Science & Business Media.
"""
@with_kw struct ValladoSolver <: AbstractLambertSolver
    M::Int = 0
    prograde::Bool = true
    maxiter::Int = 100
    rtol::Float64 = 1e-7
    stumpff_threshold::Float64 = 1e-6  # Threshold for Stumpff function series expansion
end

function SciMLBase.solve(problem::LambertProblem, solver::ValladoSolver)
    @unpack μ, r1, r2, tof = problem
    @unpack M, prograde, maxiter, rtol, stumpff_threshold = solver

    # Call the direct algorithm function
    v1, v2, numiter, retcode = vallado2013(
        μ,
        r1,
        r2,
        tof;
        M = M,
        prograde = prograde,
        maxiter = maxiter,
        rtol = rtol,
        stumpff_threshold = stumpff_threshold,
    )

    return LambertSolution(v1, v2, numiter, retcode)
end

"""
    vallado2013(μ, r1, r2, tof, M=0, prograde=true, low_path=true, maxiter=100, atol=1e-5, rtol=1e-7)

Vallado's algorithm makes use of the universal formulation to solve for the Lambert's problem.

# Arguments
- `μ`: Gravitational parameter, equivalent to GM of attractor body
- `r1`: Initial position vector
- `r2`: Final position vector
- `tof`: Time of flight
- `M`: Number of revolutions (default: 0)
- `prograde`: If true, specifies prograde motion (default: true)
- `low_path`: If two solutions are available, selects between high or low path (default: true) - Not used in implementation
- `maxiter`: Maximum number of iterations (default: 100)
- `atol`: Absolute tolerance (default: 1e-5) - Not used in implementation
- `rtol`: Relative tolerance (default: 1e-7)

# Returns
- `v1`: Initial velocity vector
- `v2`: Final velocity vector
- `numiter`: Number of iterations
- `retcode`: Return code

# References
Vallado, D. A. (2001). Fundamentals of astrodynamics and applications 
(2nd ed.). Space Technology Library. Springer Science & Business Media.
"""
function vallado2013(
    μ::Number,
    r1::Vector{<:Number},
    r2::Vector{<:Number},
    tof::Number;
    M::Int = 0,
    prograde::Bool = true,
    maxiter::Int = 100,
    rtol::Float64 = 1e-7,
    stumpff_threshold::Float64 = 1e-6,
)
    #TODO: Multi-revolution extension requires significant algorithm modifications:
    # 1. Modifying the time-of-flight equation to include 2πM√(a³/μ) term
    # 2. Adjusting initial bounds for ψ to handle multi-rev cases
    # 3. Handling the constraint that only elliptical orbits (ψ > 0) allow multi-rev
    # 4. Adding proper initial guess logic for multi-rev transfers
    # 5. Managing multiple solution branches (2M + 1 solutions for M revolutions)
    # Reference: Universal variables need extension beyond standard Vallado formulation
    # See recent universal approaches to multi-revolution Lambert problems
    (M > 0) && error(
        "Multi-revolution case not implemented - requires universal variable extension",
    )

    # Check that input parameters are safe
    assert_parameters_are_valid(μ, r1, r2, tof, M)

    # Retrieve the fundamental geometry of the problem
    r1_norm, r2_norm, c_norm, dtheta = lambert_geometry(r1, r2, prograde)

    # Compute Vallado's transfer angle parameter
    A = get_A(r1_norm, r2_norm, dtheta)
    (A == 0.0) && error("Cannot compute orbit, phase angle is 180 degrees")

    # The initial guess and limits for the bisection method
    ψ = 0.0
    ψ_low = -4 * π^2
    ψ_up = 4 * π^2

    numiter = 0
    for iter = 1:maxiter
        numiter += 1

        # Evaluate the value of y at a given psi
        y_val = y_at_psi(ψ, r1_norm, r2_norm, A, stumpff_threshold)

        if A > 0.0 && y_val < 0.0
            # Readjust psi_low until y > 0.0
            while y_val < 0.0
                ψ_low = ψ
                C2 = c2(ψ; threshold = stumpff_threshold)
                C3 = c3(ψ; threshold = stumpff_threshold)
                ψ = 0.8 * (1.0 / C3) * (1.0 - (r1_norm * r2_norm) * √C2 / A)
                y_val = y_at_psi(ψ, r1_norm, r2_norm, A, stumpff_threshold)
            end
        end

        X = X_at_psi(ψ, y_val, stumpff_threshold)
        tof_new = tof_vallado(μ, ψ, X, A, y_val, stumpff_threshold)

        # Convergence check
        if abs((tof_new - tof) / tof) < rtol
            break
        end

        # Bisection check - narrow the bounds
        condition = tof_new <= tof
        ψ_low = ψ_low + (ψ - ψ_low) * condition      # Update lower bound if tof_new too small
        ψ_up = ψ_up + (ψ - ψ_up) * (!condition)      # Update upper bound if tof_new too large

        ψ = (ψ_up + ψ_low) / 2  # New midpoint for next iteration
    end

    retcode = handle_max_iterations(numiter, maxiter)

    if retcode != :SUCCESS
        return nothing, nothing, numiter, retcode
    end

    y_val = y_at_psi(ψ, r1_norm, r2_norm, A, stumpff_threshold)
    f = 1 - y_val / r1_norm          # Lagrange coefficient f
    g = A * √(y_val / μ)             # Lagrange coefficient g  
    gdot = 1 - y_val / r2_norm       # Lagrange coefficient ġ

    v1, v2 = reconstruct_velocities_fg(f, g, gdot, r1, r2)

    return v1, v2, numiter, :SUCCESS
end

"""
    tof_vallado(μ, ψ, X, A, y)

Evaluates universal Kepler's equation.
"""
@inline function tof_vallado(
    μ::Number,
    ψ::Number,
    X::Number,
    A::Number,
    y::Number,
    stumpff_threshold::Float64,
)
    C3_val = c3(ψ; threshold = stumpff_threshold)
    tof = (X^3 * C3_val + A * √y) / √μ
    return tof
end

"""
    X_at_psi(ψ, y)

Computes the value of X at given psi.
"""
@inline function X_at_psi(ψ::Number, y::Number, stumpff_threshold::Float64)
    C2_val = c2(ψ; threshold = stumpff_threshold)
    if C2_val == 0.0
        X = 0.0
    else
        X = √(y / C2_val)
    end
    return X
end

"""
    get_A(r1_norm, r2_norm, dtheta)

Computes the value of the A constant.
"""
@inline function get_A(r1_norm::Number, r2_norm::Number, dtheta::Number)
    t_m = dtheta < π ? 1 : -1
    A = t_m * √(r1_norm * r2_norm * (1 + cos(dtheta)))
    return A
end

"""
    y_at_psi(ψ, r1_norm, r2_norm, A)

Evaluates the value of y at given psi.
"""
@inline function y_at_psi(
    ψ::Number,
    r1_norm::Number,
    r2_norm::Number,
    A::Number,
    stumpff_threshold::Float64,
)
    C2_val = c2(ψ; threshold = stumpff_threshold)
    C3_val = c3(ψ; threshold = stumpff_threshold)
    y = (r1_norm + r2_norm) + A * (ψ * C3_val - 1) / (C2_val > 0 ? √C2_val : 1.0)
    return y
end
