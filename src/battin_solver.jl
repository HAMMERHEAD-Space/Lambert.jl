export battin1984, BattinSolver

"""
    BattinSolver

Battin's elegant algorithm for solving Lambert's problem.

This algorithm improves Gauss's original method by removing the singularity 
for 180-degree transfer angles and increasing performance.

# Fields
- `M::Int`: Number of revolutions (default: 0)
- `prograde::Bool`: Direction of motion - true for prograde, false for retrograde (default: true)
- `low_path::Bool`: Path selection when multiple solutions exist (default: true)
- `maxiter::Int`: Maximum number of iterations (default: 35)
- `atol::Float64`: Absolute tolerance (default: 1e-5)
- `rtol::Float64`: Relative tolerance (default: 1e-7)

# References
Battin, R. H., & Vaughan, R. M. (1984). An elegant Lambert algorithm.
Journal of Guidance, Control, and Dynamics, 7(6), 662-670.
"""
@with_kw struct BattinSolver <: AbstractLambertSolver
    M::Int = 0
    prograde::Bool = true
    maxiter::Int = 35
    atol::Float64 = 1e-5
    levels_xi::Int = 125
    levels_K::Int = 1000
end

function SciMLBase.solve(problem::LambertProblem, solver::BattinSolver)
    @unpack μ, r1, r2, tof = problem
    @unpack M, prograde, maxiter, atol, levels_xi, levels_K = solver
    
    # Call the direct algorithm function
    v1, v2, numiter, retcode = battin1984(μ, r1, r2, tof; M=M, prograde=prograde, maxiter=maxiter, atol=atol, levels_xi=levels_xi, levels_K=levels_K)
    
    return LambertSolution(v1, v2, numiter, retcode)
end

# Helper functions translated from the Python reference
@inline function _get_λ(c::Number, s::Number, dtheta::Number)
    λ = √(s * (s - c)) / s
    return dtheta < π ? abs(λ) : -abs(λ)
end

@inline function _get_ll(λ::Number)
    return ((1 - λ) / (1 + λ))^2
end

@inline function _get_m(μ::Number, tof::Number, s::Number, λ::Number)
    return (8 * μ * tof^2) / (s^3 * (1 + λ)^6)
end

function _xi_at_x(x::Number; levels_xi::Int=125)
    η = x / (√(1 + x) + 1)^2
    
    δ = 0.0
    u = 1.0
    σ = 1.0
    m1 = 1
    
    while abs(u) > 1e-18 && m1 <= levels_xi
        m1 += 1
        γ = (m1 + 3)^2 / (4 * (m1 + 3)^2 - 1)
        δ = 1 / (1 + γ * η * δ)
        u = u * (δ - 1)
        σ = σ + u
    end
    
    return 8 * (√(1 + x) + 1) / (3 + 1 / (5 + η + (9 * η / 7) * σ))
end

@inline function _B_at_h(h1::Number, h2::Number)
    return (27 * h2) / (4 * (1 + h1)^3)
end

@inline function _u_at_B(B::Number)
    return -B / (2 * √(1 + B) + 1)
end

@inline function _u_at_h(h1::Number, h2::Number)
    return _u_at_B(_B_at_h(h1, h2))
end

function _get_h_coefficients(x::Number, ll::Number, m::Number; levels_xi::Int=125)
    xi = _xi_at_x(x; levels_xi=levels_xi)
    h_denominator = (1 + 2 * x + ll) * (4 * x + xi * (3 + x))
    h1 = ((ll + x)^2 * (1 + 3 * x + xi)) / h_denominator
    h2 = (m * (x - ll + xi)) / h_denominator
    return h1, h2
end

function _K_at_u(u::Number; levels_K::Int=1000)
    δ = 1.0
    u0 = 1.0
    σ = 1.0
    n1 = 0
    
    while abs(u0) > 1e-18 && n1 <= levels_K
        if n1 == 0
            γ = 4 / 27
            δ = 1 / (1 - γ * u * δ)
            u0 = u0 * (δ - 1)
            σ = σ + u0
        else
            for val in [1, 2]
                if val == 1
                    γ = 2 * (3 * n1 + 1) * (6 * n1 - 1) / (9 * (4 * n1 - 1) * (4 * n1 + 1))
                else
                    γ = 2 * (3 * n1 + 2) * (6 * n1 + 1) / (9 * (4 * n1 + 1) * (4 * n1 + 3))
                end
                δ = 1 / (1 - γ * u * δ)
                u0 = u0 * (δ - 1)
                σ = σ + u0
            end
        end
        n1 = n1 + 1
    end
    
    return (σ / 3)^2
end

function _battin_second_equation(u::Number, h1::Number, h2::Number; levels_K::Int=1000)
    B = _B_at_h(h1, h2)
    K = _K_at_u(u; levels_K=levels_K)
    y = ((1 + h1) / 3) * (2 + √(B + 1) / (1 - 2 * u * K))
    return y
end

function _battin_first_equation(y::Number, ll::Number, m::Number)
    x = √(((1 - ll) / 2)^2 + m / y^2) - (1 + ll) / 2
    return x
end

"""
    battin1984(μ, r1, r2, tof, M=0, prograde=true, maxiter=35, atol=1e-5)

Solves the Lambert's problem using the Battin's method.

# Arguments
- `μ`: Gravitational parameter
- `r1`: Initial position vector
- `r2`: Final position vector
- `tof`: Time of flight
- `M`: Number of revolutions
- `prograde`: True for prograde motion, false for retrograde
- `maxiter`: Maximum number of iterations
- `atol`: Absolute tolerance

# Returns
- `v1`: Initial velocity vector
- `v2`: Final velocity vector
- `numiter`: Number of iterations
- `retcode`: Return code
"""
function battin1984(
    μ::Number, 
    r1::Vector{<:Number}, 
    r2::Vector{<:Number}, 
    tof::Number; 
    M::Int=0, 
    prograde::Bool=true, 
    maxiter::Int=35, 
    atol::Float64=1e-5,
    levels_xi::Int=125,
    levels_K::Int=1000
)
    # Sanity checks
    assert_parameters_are_valid(μ, r1, r2, tof, M)
    
    # Geometry of the problem
    r1_norm, r2_norm, c_norm, dtheta = lambert_geometry(r1, r2, prograde)
    s = compute_semiperimeter(r1_norm, r2_norm, c_norm)

    # Auxiliary parameters
    λ = _get_λ(c_norm, s, dtheta)
    ll = _get_ll(λ)
    m = _get_m(μ, tof, s, λ)

    # Initial guess
    T = √(8 * μ / s^3) * tof
    T_p = (4 / 3) * (1 - λ^3)
    x0 = T > T_p ? ll : 0.0

    iter = 0
    converged = false
    x_final = x0
    
    for i in 1:maxiter
        iter = i
        
        h1, h2 = _get_h_coefficients(x0, ll, m; levels_xi=levels_xi)
        u = _u_at_h(h1, h2)
        y = _battin_second_equation(u, h1, h2; levels_K=levels_K)
        x_final = _battin_first_equation(y, ll, m)

        if check_convergence(x_final, x0, atol, 0.0)  # Only use atol for Battin
            converged = true
            break
        end
        x0 = x_final
    end
    
    retcode = converged ? :SUCCESS : handle_max_iterations(iter, maxiter)
    
    if retcode != :SUCCESS
        return nothing, nothing, iter, retcode
    end

    # Compute final velocities using converged x value
    h1, h2 = _get_h_coefficients(x_final, ll, m; levels_xi=levels_xi)
    u = _u_at_h(h1, h2)
    y = _battin_second_equation(u, h1, h2; levels_K=levels_K)
    
    r11 = (1 + λ)^2 / (4 * tof * λ)
    s11 = y * (1 + x_final)
    t11 = (m * s * (1 + λ)^2) / s11
    
    v1 = -r11 * (s11 * (r1 - r2) - t11 * r1 / r1_norm)
    v2 = -r11 * (s11 * (r1 - r2) + t11 * r2 / r2_norm)
    
    return v1, v2, iter, :SUCCESS
end
