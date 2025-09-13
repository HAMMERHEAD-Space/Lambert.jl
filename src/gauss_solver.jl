export gauss1809, GaussSolver

"""
    GaussSolver

Lambert's problem solver devised by Carl Friedrich Gauss in 1809.

This is the original method that exploits the ratio of sector to triangle area.
It shows poor accuracy and is only suitable for low transfer angles.

# Fields
- `M::Int`: Number of revolutions (default: 0)
- `prograde::Bool`: Direction of motion - true for prograde, false for retrograde (default: true)
- `maxiter::Int`: Maximum number of iterations (default: 250)
- `atol::Float64`: Absolute tolerance (default: 1e-5)
- `rtol::Float64`: Relative tolerance (default: 1e-7)

# References
Gauss, C. F. (1809). Theoria motus corporum coelestium in sectionibus conicis 
solem ambientium [Theory of the motion of the heavenly bodies moving about 
the sun in conic sections]. Hamburg: Friedrich Perthes and I. H. Besser.
"""
@with_kw struct GaussSolver <: AbstractLambertSolver
    M::Int = 0
    prograde::Bool = true
    maxiter::Int = 250
    atol::Float64 = 1e-5
    rtol::Float64 = 1e-2  # Relative tolerance for time of flight verification
    y_init::Float64 = 1.0  # Initial guess for Gauss variable y
    y_min::Float64 = 0.1   # Minimum bound for y
    y_max::Float64 = 6.283185307179586  # Maximum bound for y (2π)
end

function solve(problem::LambertProblem, solver::GaussSolver)
    @unpack μ, r1, r2, tof = problem
    @unpack M, prograde, maxiter, atol, rtol, y_init, y_min, y_max = solver

    v1, v2, numiter, status = gauss1809(
        μ,
        r1,
        r2,
        tof;
        M = M,
        prograde = prograde,
        maxiter = maxiter,
        atol = atol,
        rtol = rtol,
        y_init = y_init,
        y_min = y_min,
        y_max = y_max,
    )

    return LambertSolution(v1, v2, numiter, status)
end

"""
    gauss1809(μ, r1, r2, tof; kwargs...)

Solve Lambert's problem using Gauss' 1809 method following Bate's book.

# Arguments
- `μ`: Gravitational parameter
- `r1`: Initial position vector  
- `r2`: Final position vector
- `tof`: Time of flight

# Keyword Arguments
- `M`: Number of revolutions (default: 0)
- `prograde`: Prograde motion flag (default: true)
- `maxiter`: Maximum iterations (default: 250)
- `atol`: Absolute tolerance (default: 1e-5)
- `rtol`: Relative tolerance (default: 1e-2)
- `y_init`: Initial guess for y variable (default: 1.0)
- `y_min`: Minimum bound for y (default: 0.1)
- `y_max`: Maximum bound for y (default: 2π)

# Returns
- `v1`: Initial velocity vector
- `v2`: Final velocity vector
- `numiter`: Number of iterations
- `retcode`: Return code
"""
function gauss1809(
    μ::Number,
    r1::Vector{<:Number},
    r2::Vector{<:Number},
    tof::Number;
    M::Int = 0,
    prograde::Bool = true,
    maxiter::Int = 250,
    atol::Float64 = 1e-5,
    rtol::Float64 = 1e-2,
    y_init::Float64 = 1.0,
    y_min::Float64 = 0.1,
    y_max::Float64 = 6.283185307179586,
)
    # Sanity checks
    assert_parameters_are_valid(μ, r1, r2, tof, M)

    # Proper Gauss 1809 method following Python lamberthub implementation
    r1_norm, r2_norm, _, dtheta = lambert_geometry(r1, r2, prograde)

    # Compute auxiliary constants s and w (Bate's book equations 5.6-2 and 5.6-3)
    s = (r1_norm + r2_norm) / (4 * sqrt(r1_norm * r2_norm) * cos(dtheta / 2)) - 0.5
    w = (μ * tof^2) / (2 * sqrt(r1_norm * r2_norm) * cos(dtheta / 2))^3

    # Initial guess for y (dependent variable)
    y0 = 1.0

    # Gauss iterative procedure
    for numiter = 1:maxiter
        # Compute x using Gauss first equation (5.6-13)
        x = w / y0^2 - s

        # Compute new y using Gauss second equation (5.6-14)
        y = 1 + X_at_x(x) * (s + x)

        # Check convergence
        if abs(y - y0) <= atol
            # Determine orbital type and compute delta anomaly
            deltaAnomaly = if x > 0
                # Elliptical orbit
                2 * acos(1 - 2 * x)
            elseif x == 0
                # Parabolic orbit
                π  # This case needs special handling
            else
                # Hyperbolic orbit
                2 * acosh(1 - 2 * x)
            end

            # Compute orbital parameter p
            p =
                (r1_norm * r2_norm * (1 - cos(dtheta))) / (
                    r1_norm + r2_norm -
                    2 * sqrt(r1_norm * r2_norm) * cos(dtheta / 2) * cos(deltaAnomaly / 2)
                )

            # Compute Lagrange coefficients f, g, f_dot, g_dot
            f = 1 - r2_norm * (1 - cos(dtheta)) / p
            g = (r1_norm * r2_norm * sin(dtheta)) / sqrt(μ * p)
            f_dot =
                sqrt(μ / p) *
                tan(dtheta / 2) *
                ((1 - cos(dtheta)) / p - 1 / r1_norm - 1 / r2_norm)
            g_dot = 1 - r1_norm * (1 - cos(dtheta)) / p

            # Compute velocity vectors
            v1 = (r2 - f * r1) / g
            v2 = f_dot * r1 + g_dot * v1

            return v1, v2, numiter, :SUCCESS
        else
            # Update y0 for next iteration
            y0 = y
        end
    end

    return zeros(3), zeros(3), maxiter, :MAXIMUM_ITERATIONS
end

# Helper function for X series (equation 5.6-15 from Bate's book)
function X_at_x(x::Number; order::Int = 50)
    coefficients = [1.0]
    for n = 3:(3+order-1)
        coeff = (2 * n) / (2 * n - 1)
        push!(coefficients, coefficients[end] * coeff * x)
    end
    X = (4 / 3) * sum(coefficients)
    return X
end
