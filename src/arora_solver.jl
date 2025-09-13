export arora2013, AroraSolver

"""
    AroraSolver

Lambert's problem solver using Arora's cosine transformation method from 2013.

This algorithm exploits the universal formulae by defining a new cosine-based 
transformation and developing a robust initial guess.

# Fields
- `M::Int`: Number of revolutions (default: 0) - Multi-revolution supported
- `prograde::Bool`: Direction of motion - true for prograde, false for retrograde (default: true)
- `maxiter::Int`: Maximum number of iterations (default: 35)
- `atol::Float64`: Absolute tolerance (default: 1e-5)

# References
Arora, N., & Russell, R. P. (2013). A fast and robust multiple revolution Lambert 
algorithm using a cosine transformation. In AAS/AIAA Astrodynamics Specialist 
Conference, Hilton Head, SC. Paper AAS 13-728.
"""
@with_kw struct AroraSolver <: AbstractLambertSolver
    M::Int = 0
    prograde::Bool = true
    maxiter::Int = 35
    atol::Float64 = 1e-5
    ε::Float64 = 2e-2
end

function SciMLBase.solve(problem::LambertProblem, solver::AroraSolver)
    @unpack μ, r1, r2, tof = problem
    @unpack M, prograde, maxiter, atol, ε = solver

    # Call the direct algorithm function
    v1, v2, numiter, retcode = arora2013(
        μ,
        r1,
        r2,
        tof;
        M = M,
        prograde = prograde,
        maxiter = maxiter,
        atol = atol,
        ε = ε,
    )

    return LambertSolution(v1, v2, numiter, retcode)
end

"""
    arora2013(μ, r1, r2, tof, M=0, prograde=true, maxiter=35, atol=1e-5)

Arora and Russell's algorithm for solving Lambert's problem.

# Arguments
- `μ`: Gravitational parameter, equivalent to GM of attractor body
- `r1`: Initial position vector
- `r2`: Final position vector
- `tof`: Time of flight
- `M`: Number of revolutions (default: 0) - Multi-revolution supported
- `prograde`: If true, specifies prograde motion (default: true)
- `maxiter`: Maximum number of iterations (default: 35)
- `atol`: Absolute tolerance (default: 1e-5)

# Returns
- `v1`: Initial velocity vector
- `v2`: Final velocity vector
- `numiter`: Number of iterations
- `retcode`: Return code

# References
Arora, N., & Russell, R. P. (2013). A fast and robust multiple revolution Lambert 
algorithm using a cosine transformation. In AAS/AIAA Astrodynamics Specialist 
Conference, Hilton Head, SC. Paper AAS 13-728.
"""
function arora2013(
    μ::Number,
    r1::Vector{<:Number},
    r2::Vector{<:Number},
    tof::Number;
    M::Int = 0,
    prograde::Bool = true,
    maxiter::Int = 35,
    atol::Float64 = 1e-5,
    ε::Float64 = 2e-2,
)
    # Multi-revolution support enabled - algorithm was specifically designed for this case
    # Reference: Arora & Russell (2013) "A fast and robust multiple revolution Lambert algorithm"

    # Check that input parameters are safe
    assert_parameters_are_valid(μ, r1, r2, tof, M)

    # Compute basic geometry
    r1_norm, _, _, dtheta = lambert_geometry(r1, r2, prograde)
    assert_transfer_angle_not_zero(dtheta)
    d = dtheta <= π ? 1 : -1

    # Compute the characteristic length and time variables
    L_ref = r1_norm
    T_ref = √(L_ref^3 / μ)

    # All variables named "hat" have been non-dimensionalized
    r1_hat = r1 / L_ref
    r2_hat = r2 / L_ref
    r1hat_norm = norm(r1_hat)
    r2hat_norm = norm(r2_hat)
    tof_hat = tof / T_ref
    μ_hat = μ * (T_ref^2 / L_ref^3)

    # Auxiliary variables
    S = √((r1hat_norm + r2hat_norm)^3 / μ_hat)

    # Solve Lambert's geometry parameter
    τ = d * √(r1hat_norm * r2hat_norm * (1 + cos(dtheta))) / (r1hat_norm + r2hat_norm)

    # Compute the parabolic time of flight
    tof_p = S * √(1 - √2 * τ) * (τ + √2) / 3

    # Compare actual and parabolic TOF to generate the initial guess
    if tof_hat <= tof_p
        # The final orbit is expected to be hyperbolic
        tof20 = S * √(1 - 20 * τ) * (τ + 0.04940968903 * (1 - 20 * τ))
        tof100 = S * √(1 - 100 * τ) * (τ + 0.00999209404 * (1 - 100 * τ))

        if d == 1
            # Value from the first row of (Table 2)
            k_n = √2
            k_m = 1 / τ
            k_i = (k_n + k_m) / 2
            Z = 1.0 / √2
            α = 1 / 2
            F_0 = tof_p
            F_1 = 0.0
            W = get_W(k_i, M; ε = ε)
            F_i = get_TOF(k_i, τ, S, W)
            F_star = tof_hat

            # Compute the Sundman transformation and the initial guess
            x_star = get_x(F_0, F_1, F_i, F_star, Z, α)
            k = k_n + (k_m - k_n) * x_star

        elseif (d == -1) && (tof_hat > tof20)
            # Initial guess falls in hyperbolic region H1
            k_n = √2
            k_m = 20
            k_i = (2 * k_n + k_m) / 3
            Z = 1 / 3
            α = 1
            F_0 = tof_p
            F_1 = tof20
            W = get_W(k_i, M; ε = ε)
            F_i = get_TOF(k_i, τ, S, W)
            F_star = tof_hat

            # Compute the Sundman transformation and the initial guess
            x_star = get_x(F_0, F_1, F_i, F_star, Z, α)
            k = k_n + (k_m - k_n) * x_star

        elseif (d == -1) && (tof_hat <= tof20)
            # Initial guess falls in the H2 region
            t_star = tof_hat
            t_0 = tof20
            t_1 = tof100
            k =
                (
                    (t_1 * (t_0 - t_star) * 10 - t_0 * √20 * (t_1 - t_star)) /
                    (t_star * (t_0 - t_1))
                )^2
        end
    else
        # The final orbit is expected to be an ellipse
        if M == 0
            # A collection of auxiliary independent variable values
            k_set = [-1.41, -1.38, -1.00, -1/2, 0, 1/√2]

            # Precomputed values of W(k) required for the zero-rev initial guess
            W_set = [
                4839.684497246,
                212.087279879,
                5.712388981,
                1.954946607,
                1.110720735,
                0.6686397730,
            ]

            # Time of flight for each one of the previous auxiliary values
            t_set = [get_TOF(k, τ, S, W) for (k, W) in zip(k_set, W_set)]
            tof_m141, tof_m138, tof_m1, tof_m1half, tof_0, tof_1oversq2 = t_set

            # Filter out the region
            if tof_hat <= tof_0
                # Region E1 applies
                k_n = 0
                k_m = √2
                k_i = 1/√2
                Z, α = 1/2, 1
                F_0 = tof_0
                F_1 = tof_p
                F_i = tof_1oversq2
                F_star = tof_hat
                x_star = get_x(F_0, F_1, F_i, F_star, Z, α)
                k = k_n + (k_m - k_n) * x_star

            elseif tof_0 <= tof_hat <= tof_m1
                # Region E2 applies
                k_n = 0
                k_m = -1
                k_i = -1/2
                Z, α = 1/2, 1
                F_0 = tof_0
                F_1 = tof_m1
                F_i = tof_m1half
                F_star = tof_hat
                x_star = get_x(F_0, F_1, F_i, F_star, Z, α)
                k = k_n + (k_m - k_n) * x_star

            elseif tof_m1 <= tof_hat <= tof_m138
                # Region E3 applies
                k_n = -1
                k_m = -√2
                k_i = -1.38
                c1 = 540649/3125
                c2 = 256
                c3 = 1
                c4 = 1
                α = 16
                F_n = tof_m1^(-1)
                F_i = tof_m138^(-1)
                F_star = tof_hat^(-1)
                γ1, γ2, γ3 = get_gammas(F_i, F_n, F_star)
                k =
                    -c4 *
                    (
                        ((γ1 * c1 - c3 * γ3) * c2 + c3 * c1 * γ2) /
                        (γ3 * c1 - c3 * γ1 - γ2 * c2)
                    )^(1/α)

            else
                # Region E4 applies
                k_n = -1.38
                k_m = -√2
                k_i = -1.41
                c1 = 49267/27059
                c2 = 67286/17897
                c3 = 2813/287443
                c4 = 4439/3156
                α = 243
                F_n = tof_m138^(-1)
                γ1, γ2, γ3 = get_gammas(F_i, F_n, F_star)
                k =
                    -c4 *
                    (
                        ((γ1 * c1 - c3 * γ3) * c2 + c3 * c1 * γ2) /
                        (γ3 * c1 - c3 * γ1 - γ2 * c2)
                    )^(1/α)
            end
        else
            # Multi-revolution case (M > 0)
            # Use a simple initial guess based on the number of revolutions
            k = -1.0 - 0.1 * M  # Simple heuristic for multi-rev cases
        end
    end

    # Iterative process using Halley's method
    numiter = 0
    for iter = 1:maxiter
        numiter = iter

        # Evaluate the auxiliary function and its derivatives
        W = get_W(k, M; ε = ε)
        Wp = get_Wprime(k, W; ε = ε)
        Wpp = get_W2prime(k, W, Wp; ε = ε)

        # Evaluate the predicted time of flight with the current value of k
        c = (1 - k * τ) / τ
        tofc = get_TOF(k, τ, S, W)

        # Check computed time of flight matches target one
        if abs(tof_hat - tofc) <= atol
            break
        end

        # Compute the time derivatives to proceed by using Halley's method
        tofc_p = (-tofc / (2 * c)) + S * τ * √(c * τ) * (Wp * c - W)
        tofc_pp = (-tofc / (4 * c^2)) + S * τ * √(c * τ) * (W / c + c * Wpp - 3 * Wp)

        # Solve Halley's step
        deltak = -(tofc - tof_hat) / (tofc_p - (tofc - tof_hat) * tofc_pp / (2.0 * tofc_p))

        # Update the value of the independent variable
        k += deltak

        # Bound the value of k if required
        if k < -√2
            k = -√2 + 1e-12  # Lower bound for elliptical orbits
        end

        if k < √2 && tof_hat < tof_p
            k = √2 + 1e-12  # Ensure hyperbolic regime for short transfer times
        end

        if k > √2 && tof_hat > tof_p
            k = √2 - 1e-12  # Ensure elliptical regime for long transfer times
        end

        if tof_hat < tof_p && d > 0 && (1 - τ * k) < 0
            k = 1 / τ - 1e-12  # Prevent division by zero in y calculation
        end
    end

    retcode = handle_max_iterations(numiter, maxiter)

    if retcode != :SUCCESS
        return nothing, nothing, numiter, retcode
    end

    # Evaluate f and g functions
    f = 1 - (1 - k * τ) * (r1hat_norm + r2hat_norm) / r1hat_norm
    g = S * τ * √((1 - k * τ) * μ_hat)
    g_dot = 1 - (1 - k * τ) * (r1hat_norm + r2hat_norm) / r2hat_norm

    # Compute the initial and final velocity vectors using f&g reconstruction
    v1_hat, v2_hat = reconstruct_velocities_fg(f, g, g_dot, r1_hat, r2_hat)
    v1 = v1_hat * (L_ref / T_ref)
    v2 = v2_hat * (L_ref / T_ref)

    return v1, v2, numiter, :SUCCESS
end

# Helper functions for Arora solver

@inline function get_gammas(F_i::Number, F_n::Number, F_star::Number)
    γ1 = F_i * (F_star - F_n)
    γ2 = F_star * (F_n - F_i)
    γ3 = F_n * (F_star - F_i)
    return γ1, γ2, γ3
end

@inline function get_x(
    F_0::Number,
    F_1::Number,
    F_i::Number,
    F_star::Number,
    Z::Number,
    α::Number,
)
    x =
        (
            (Z * (F_0 - F_star) * (F_1 - F_i)) /
            ((F_i - F_star) * (F_1 - F_0) * Z + (F_0 - F_i) * (F_1 - F_star))
        )^(1/α)
    return x
end

function get_W(k::Number, M::Int; ε::Float64 = 2e-2)
    # Evaluate the sign of k
    m = 2 - k^2
    sgn_k = sign(k)
    sq2 = √2

    if -sq2 <= k < (sq2 - ε)
        # Elliptical orbits
        W = ((1 - sgn_k) * π + sgn_k * acos(1 - m) + 2 * π * M) / √(m^3) - k / m
    elseif k > sq2 + ε
        # Hyperbolic orbits
        W = -acosh(1 - m) / √(-(m^3)) - k / m
    elseif sq2 - ε <= k <= sq2 + ε
        # Direct transfer, no complete revolutions (M = 0)
        v = k - sq2
        v2, v3, v4, v5, v6, v7, v8 = [v^i for i = 2:8]

        W =
            (√2 / 3) - (1 / 5) * v + (2 / 35) * sq2 * v2 - (2 / 63) * v3 +
            (2 / 231) * sq2 * v4 - (2 / 429) * v5 + (8 / 6435) * sq2 * v6 -
            (8 / 12155) * v7 + (8 / 46189) * sq2 * v8
    else
        throw(ValueError("Did not find a suitable equation to find W!"))
    end

    return W
end

function get_Wprime(k::Number, W::Number; ε::Float64 = 2e-2)
    # Evaluate m
    m = 2 - k^2

    if k < √2 - ε
        W_prime = (-2 + 3 * W * k) / m
    elseif √2 - ε < k < √2 + ε
        # Series derivative
        sq2 = √2
        v = k - sq2
        v2, v3, v4, v5, v6, v7 = [v^i for i = 2:7]

        W_prime =
            -1/5 + sq2 * v * (4/35) - v2 * (6/63) + sq2 * v3 * (8/231) - v4 * (10/429) +
            sq2 * v5 * (48/6435) - v6 * (56/12155) + sq2 * v7 * (64/46189)
    else
        W_prime = (-2 + 3 * W * k) / m
    end

    return W_prime
end

function get_W2prime(k::Number, W::Number, W_prime::Number; ε::Float64 = 2e-2)
    # Evaluate m
    m = 2 - k^2

    if k < √2 - ε
        W_2prime = (5 * W_prime * k + 3 * W) / m
    elseif √2 - ε < k < √2 + ε
        # Series second derivative
        sq2 = √2
        v = k - sq2
        v2, v3, v4, v5, v6 = [v^i for i = 2:6]

        W_2prime =
            sq2 * (4/35) - v * (12/63) + sq2 * v2 * (24/231) - v3 * (40/429) +
            sq2 * v4 * (240/6435) - v5 * (336/12155) + sq2 * v6 * (448/46189)
    else
        W_2prime = (5 * W_prime * k + 3 * W) / m
    end

    return W_2prime
end

@inline function get_TOF(k::Number, τ::Number, S::Number, W::Number)
    TOF = S * √(1 - k * τ) * (τ + (1 - k * τ) * W)
    return TOF
end
