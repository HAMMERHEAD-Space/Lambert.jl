export izzo2015, IzzoSolver

"""
    IzzoSolver

Lambert's problem solver using Izzo's devised algorithm from 2015.

# Fields
- `M::Int`: Number of revolutions (default: 0)
- `prograde::Bool`: Direction of motion - true for prograde, false for retrograde (default: true)
- `low_path::Bool`: Path selection when multiple solutions exist (default: true)
- `maxiter::Int`: Maximum number of iterations (default: 35)
- `atol::Float64`: Absolute tolerance (default: 1e-5)
- `rtol::Float64`: Relative tolerance (default: 1e-7)

# References
Izzo, D. (2015). Revisiting Lambert's problem. Celestial Mechanics
and Dynamical Astronomy, 121(1), 1-15.
"""
@with_kw struct IzzoSolver <: AbstractLambertSolver
    M::Int = 0
    prograde::Bool = true
    low_path::Bool = true
    maxiter::Int = 35
    atol::Float64 = 1e-5
    rtol::Float64 = 1e-7
end

function SciMLBase.solve(problem::LambertProblem, solver::IzzoSolver)
    @unpack μ, r1, r2, tof = problem
    @unpack M, prograde, low_path, maxiter, atol, rtol = solver
    
    # Call the direct algorithm function
    v1, v2 = izzo2015(μ, r1, r2, tof; M=M, prograde=prograde, low_path=low_path, maxiter=maxiter, atol=atol, rtol=rtol)
    
    return LambertSolution(v1, v2, 0, :SUCCESS)
end

"""
    izzo2015(μ, r1, r2, tof, M=0, prograde=true, low_path=true, maxiter=35, atol=1e-5, rtol=1e-7)

Solves Lambert problem using Izzo's devised algorithm.

# Arguments
- `μ`: Gravitational parameter, equivalent to GM of attractor body
- `r1`: Initial position vector
- `r2`: Final position vector
- `tof`: Time of flight
- `M`: Number of revolutions (default: 0)
- `prograde`: If true, specifies prograde motion (default: true)
- `low_path`: If two solutions are available, selects between high or low path (default: true)
- `maxiter`: Maximum number of iterations (default: 35)
- `atol`: Absolute tolerance (default: 1e-5)
- `rtol`: Relative tolerance (default: 1e-7)

# Returns
- `v1`: Initial velocity vector
- `v2`: Final velocity vector

# References
Izzo, D. (2015). Revisiting Lambert's problem. Celestial Mechanics
and Dynamical Astronomy, 121(1), 1-15.
"""
function izzo2015(
    μ::Number, r1::Vector{<:Number}, r2::Vector{<:Number}, tof::Number;
    M::Int=0, prograde::Bool=true, low_path::Bool=true, maxiter::Int=35, atol::Float64=1e-5, rtol::Float64=1e-7
)
    # Check that input parameters are safe
    assert_parameters_are_valid(μ, r1, r2, tof, M)
    
    # Compute basic geometry
    r1_norm, r2_norm, c_norm, dtheta = lambert_geometry(r1, r2, prograde)
    
    # Semiperimeter
    s = compute_semiperimeter(r1_norm, r2_norm, c_norm)
    
    # Unit vectors
    i_r1, i_r2, i_h_unnorm = compute_unit_vectors(r1, r2, r1_norm, r2_norm)
    i_h = i_h_unnorm / norm(i_h_unnorm)
    
    # Geometry of the problem
    ll = √(1 - min(1.0, c_norm / s))
    
    # Compute the fundamental tangential directions
    if i_h[3] < 0
        ll = -ll
        i_t1 = cross(i_r1, i_h)
        i_t2 = cross(i_r2, i_h)
    else
        i_t1, i_t2 = compute_tangential_unit_vectors(i_h, i_r1, i_r2)
    end
    
    # Correct transfer angle parameter and tangential vectors
    if !prograde
        ll = -ll
        i_t1 = -i_t1
        i_t2 = -i_t2
    end
    
    # Non dimensional time of flight
    T = √(2 * μ / s^3) * tof
    
    # Find solutions and filter them
    x, y = find_xy(ll, T, M, maxiter, atol, rtol, low_path)
    
    # Reconstruct
    γ = √(μ * s / 2)
    ρ = (r1_norm - r2_norm) / c_norm
    σ = √(1 - ρ^2)
    
    # Compute the radial and tangential components
    V_r1, V_r2, V_t1, V_t2 = reconstruct(x, y, r1_norm, r2_norm, ll, γ, ρ, σ)
    
    # Solve for the initial and final velocity
    v1 = V_r1 * (r1 / r1_norm) + V_t1 * i_t1
    v2 = V_r2 * (r2 / r2_norm) + V_t2 * i_t2
    
    return v1, v2
end

"""
    reconstruct(x, y, r1, r2, ll, γ, ρ, σ)

Reconstruct solution velocity vectors.
"""
@inline function reconstruct(x::Number, y::Number, r1::Number, r2::Number, ll::Number, γ::Number, ρ::Number, σ::Number)
    V_r1 = γ * ((ll * y - x) - ρ * (ll * y + x)) / r1
    V_r2 = -γ * ((ll * y - x) + ρ * (ll * y + x)) / r2
    V_t1 = γ * σ * (y + ll * x) / r1
    V_t2 = γ * σ * (y + ll * x) / r2
    return V_r1, V_r2, V_t1, V_t2
end

"""
    find_xy(ll, T, M, maxiter, atol, rtol, low_path)

Computes all x, y for given number of revolutions.
"""
function find_xy(ll::Number, T::Number, M::Int, maxiter::Int, atol::Float64, rtol::Float64, low_path::Bool)
    # For abs(ll) == 1 the derivative is not continuous
    @assert abs(ll) < 1 "abs(ll) must be < 1"
    
    M_max = Int(floor(T / π))                     # Maximum possible revolutions from time
    T_00 = acos(ll) + ll * √(1 - ll^2)            # Single-rev minimum time (T_0M)
    
    # Refine maximum number of revolutions if necessary
    if T < T_00 + M_max * π && M_max > 0          # Check if time is sufficient
        _, T_min = compute_T_min(ll, M_max, maxiter, atol, rtol)
        if T < T_min                              # Time too short for M_max revs
            M_max -= 1
        end
    end
    
    # Check if a feasible solution exist for the given number of revolutions
    if M > M_max
        error("No feasible solution, try lower M!")
    end
    
    # Initial guess
    x_0 = initial_guess(T, ll, M, low_path)
    
    # Start Householder iterations from x_0 and find x, y
    x = householder(x_0, T, ll, M, atol, rtol, maxiter)
    y = compute_y(x, ll)
    
    return x, y
end

"""
    compute_y(x, ll)

Computes y.
"""
@inline function compute_y(x::Number, ll::Number)
    return √(1 - ll^2 * (1 - x^2))
end

"""
    compute_psi(x, y, ll)

Computes psi using the appropriate inverse function.
"""
function compute_psi(x::Number, y::Number, ll::Number)
    if -1 <= x < 1
        # Elliptic motion - use arc cosine to avoid numerical errors
        return acos(x * y + ll * (1 - x^2))
    elseif x > 1
        # Hyperbolic motion - the hyperbolic sine is bijective
        return asinh((y - x * ll) * √(x^2 - 1))
    else
        # Parabolic motion
        return 0.0
    end
end

"""
    tof_equation(x, T0, ll, M)

Time of flight equation.
"""
function tof_equation(x::Number, T0::Number, ll::Number, M::Int)
    return tof_equation_y(x, compute_y(x, ll), T0, ll, M)
end

"""
    tof_equation_y(x, y, T0, ll, M)

Time of flight equation with externally computed y.
"""
function tof_equation_y(x::Number, y::Number, T0::Number, ll::Number, M::Int)
    if M == 0 && √(0.6) < x < √(1.4)
        η = y - ll * x
        S_1 = (1 - ll - x * η) * 0.5
        Q = 4 / 3 * hyp2f1b(S_1)
        T_ = (η^3 * Q + 4 * ll * η) * 0.5
    else
        ψ = compute_psi(x, y, ll)
        T_ = (ψ + M * π) / √(abs(1 - x^2)) - x + ll * y
        T_ = T_ / (1 - x^2)
    end
    
    return T_ - T0
end

"""
    tof_equation_p(x, y, T, ll)

First derivative of time of flight equation.
"""
function tof_equation_p(x::Number, y::Number, T::Number, ll::Number)
    return (3 * T * x - 2 + 2 * ll^3 * x / y) / (1 - x^2)
end

"""
    tof_equation_p2(x, y, T, dT, ll)

Second derivative of time of flight equation.
"""
function tof_equation_p2(x::Number, y::Number, T::Number, dT::Number, ll::Number)
    return (3 * T + 5 * x * dT + 2 * (1 - ll^2) * ll^3 / y^3) / (1 - x^2)
end

"""
    tof_equation_p3(x, y, dT, ddT, ll)

Third derivative of time of flight equation.
"""
function tof_equation_p3(x::Number, y::Number, dT::Number, ddT::Number, ll::Number)
    return (7 * x * ddT + 8 * dT - 6 * (1 - ll^2) * ll^5 * x / y^5) / (1 - x^2)
end

"""
    compute_T_min(ll, M, maxiter, atol, rtol)

Compute minimum T.
"""
function compute_T_min(ll::Number, M::Int, maxiter::Int, atol::Float64, rtol::Float64)
    if ll == 1
        x_T_min = 0.0
        T_min = tof_equation(x_T_min, 0.0, ll, M)
    else
        if M == 0
            x_T_min = Inf
            T_min = 0.0
        else
            # Set x_i > 0 to avoid problems at ll = -1
            x_i = 0.1
            T_i = tof_equation(x_i, 0.0, ll, M)
            x_T_min = halley(x_i, T_i, ll, atol, rtol, maxiter)
            T_min = tof_equation(x_T_min, 0.0, ll, M)
        end
    end
    
    return x_T_min, T_min
end

"""
    initial_guess(T, ll, M, low_path)

Compute initial guess for the solution.
"""
function initial_guess(T::Number, ll::Number, M::Int, low_path::Bool)
    if M == 0
        # Single revolution
        T_0 = acos(ll) + ll * √(1 - ll^2) + M * π  # Equation 19
        T_1 = 2 * (1 - ll^3) / 3  # Equation 21
        
        if T >= T_0
            x_0 = (T_0 / T)^(2/3) - 1
        elseif T < T_1
            x_0 = 5 / 2 * T_1 / T * (T_1 - T) / (1 - ll^5) + 1
        else
            # This is the condition T_1 < T < T_0
            x_0 = exp(log(2) * log(T / T_0) / log(T_1 / T_0)) - 1
        end
        
        return x_0
    else
        # Multiple revolution
        x_0l = (((M * π + π) / (8 * T))^(2/3) - 1) / (((M * π + π) / (8 * T))^(2/3) + 1)
        x_0r = (((8 * T) / (M * π))^(2/3) - 1) / (((8 * T) / (M * π))^(2/3) + 1)
        
        # Filter out the solution
        x_0 = low_path ? max(x_0l, x_0r) : min(x_0l, x_0r)
        
        return x_0
    end
end

"""
    halley(p0, T0, ll, atol, rtol, maxiter)

Find a minimum of time of flight equation using the Halley method.
"""
function halley(p0::Number, T0::Number, ll::Number, atol::Float64, rtol::Float64, maxiter::Int)
    for iter in 1:maxiter
        y = compute_y(p0, ll)
        fder = tof_equation_p(p0, y, T0, ll)
        fder2 = tof_equation_p2(p0, y, T0, fder, ll)
        if fder2 == 0
            error("Derivative was zero")
        end
        fder3 = tof_equation_p3(p0, y, fder, fder2, ll)
        
        # Halley step (cubic)
        p = p0 - 2 * fder * fder2 / (2 * fder2^2 - fder * fder3)
        
        if abs(p - p0) < rtol * abs(p0) + atol
            return p
        end
        p0 = p
    end
    
    error("Halley: Failed to converge")
end

"""
    householder(p0, T0, ll, M, atol, rtol, maxiter)

Find a zero of time of flight equation using the Householder method.
"""
function householder(p0::Number, T0::Number, ll::Number, M::Int, atol::Float64, rtol::Float64, maxiter::Int)
    for iter in 1:maxiter
        y = compute_y(p0, ll)
        fval = tof_equation_y(p0, y, T0, ll, M)
        T = fval + T0
        fder = tof_equation_p(p0, y, T, ll)
        fder2 = tof_equation_p2(p0, y, T, fder, ll)
        fder3 = tof_equation_p3(p0, y, fder, fder2, ll)
        
        # Householder step (quartic)
        p = p0 - fval * (
            (fder^2 - fval * fder2 / 2) /
            (fder * (fder^2 - fval * fder2) + fder3 * fval^2 / 6)
        )
        
        if abs(p - p0) < rtol * abs(p0) + atol
            return p
        end
        p0 = p
    end
    
    error("Householder: Failed to converge")
end

"""
    hyp2f1b(x)

Hypergeometric function 2F1(3, 1, 5/2, x).
"""
function hyp2f1b(x::Number)
    if x >= 1.0
        return Inf
    else
        res = 1.0
        term = 1.0
        ii = 0
        while true
            term = term * (3 + ii) * (1 + ii) / (5/2 + ii) * x / (ii + 1)
            res_old = res
            res += term
            if res_old == res
                return res
            end
            ii += 1
        end
    end
end
