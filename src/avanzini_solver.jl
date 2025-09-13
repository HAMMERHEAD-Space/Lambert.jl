export avanzini2008, AvanziniSolver

"""
    AvanziniSolver

Lambert's problem solver using Avanzini's eccentricity-based method from 2008.

This solver takes advantage of the conservation of the eccentricity projection 
along the chord direction. Note: This solver does not support multi-revolutions.

# Fields
- `M::Int`: Number of revolutions (default: 0) - Only M=0 supported
- `prograde::Bool`: Direction of motion - true for prograde, false for retrograde (default: true)
- `maxiter::Int`: Maximum number of iterations (default: 35)
- `atol::Float64`: Absolute tolerance (default: 1e-5)
- `rtol::Float64`: Relative tolerance (default: 1e-7)

# References
Avanzini, G. (2008). A simple Lambert algorithm. Journal of Guidance, 
Control, and Dynamics, 31(6), 1587-1594.
"""
@with_kw struct AvanziniSolver <: AbstractLambertSolver
    M::Int = 0
    prograde::Bool = true
    maxiter::Int = 35
    atol::Float64 = 1e-5
    dx::Float64 = 1e-3  # Finite difference step for derivative
    x_init::Float64 = 0.0  # Initial guess for x parameter
end

function SciMLBase.solve(problem::LambertProblem, solver::AvanziniSolver)
    @unpack μ, r1, r2, tof = problem
    @unpack M, prograde, maxiter, atol, dx, x_init = solver

    # Call the direct algorithm function
    v1, v2 = avanzini2008(
        μ,
        r1,
        r2,
        tof;
        M = M,
        prograde = prograde,
        maxiter = maxiter,
        atol = atol,
        dx = dx,
        x_init = x_init,
    )

    return LambertSolution(v1, v2, 0, :SUCCESS)
end

"""
    avanzini2008(μ, r1, r2, tof, M=0, prograde=true, maxiter=35, atol=1e-5)

Solves the Lambert's problem using the Avanzini's method.

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
"""
function avanzini2008(
    μ::Number,
    r1::Vector{<:Number},
    r2::Vector{<:Number},
    tof::Number;
    M::Int = 0,
    prograde::Bool = true,
    maxiter::Int = 35,
    atol::Float64 = 1e-5,
    dx::Float64 = 1e-3,
    x_init::Float64 = 0.0,
)
    # Check multi-revolution case
    # Multi-revolution extension possible but complex - requires:
    # 1. Modifying the eccentricity parameterization to handle transfer angles > 2π
    # 2. Extending the transverse eccentricity method to multi-rev cases (see He et al.)
    # 3. Adding proper initial guess logic for multi-rev transfers
    # 4. Handling the multiplicity of solutions (2M + 1 solutions for M revolutions)
    # Reference: He et al. extension of Avanzini's transverse eccentricity approach
    (M > 0) && error("Avanzini solver does not support multi-revolution scenarios")

    # Get geometry following Python reference exactly
    r1_norm = norm(r1)
    r2_norm = norm(r2)
    c_norm = norm(r2 - r1)
    dtheta = get_transfer_angle(r1, r2, prograde)
    w_c = get_transfer_angle(r1, (r2 - r1), prograde)

    # Get fundamental ellipse properties (Python reference)
    ecc_F = (r1_norm - r2_norm) / c_norm
    a_F = (r1_norm + r2_norm) / 2
    p_F = a_F * (1 - ecc_F^2)

    # Get eccentricity limits (Python reference)
    ecc_max = -1 / cos(dtheta / 2)
    ecc_H = dtheta > π ? sqrt(ecc_max^2 - ecc_F^2) : Inf
    ecc_P = sqrt(1 - ecc_F^2)

    # Define the function to find zero of (Python _f function exactly)
    function f(x)
        # Compute ecc_T at this x based on transfer angle (Python _get_eccT_at_x)
        ecc_T = if dtheta > π
            X = exp((1/ecc_H + 1/ecc_P) * x)
            ecc_P * ecc_H * (X - 1) / (ecc_P + ecc_H * X)
        else
            ecc_P * (1 - exp(-x / ecc_P))
        end

        # Compute transfer orbit parameters (Python eap_from_eccT)
        ecc = sqrt(ecc_T^2 + ecc_F^2)
        p = p_F - ecc_T * r1_norm * r2_norm * sin(dtheta) / c_norm
        a = p / (1 - ecc_F^2 - ecc_T^2)  # Python formula

        # Compute argument of periapsis (Python w_at_eccT)
        y = ecc_F * sin(w_c) + ecc_T * cos(w_c)
        x_coord = ecc_F * cos(w_c) - ecc_T * sin(w_c)
        w = atan(y, x_coord)

        # True anomalies (Python get_true_anomalies)
        nu_1 = -w
        nu_2 = nu_1 + dtheta

        # Mean anomalies and time of flight (Python kepler_tof_at_eccT)
        M_1 = trueAnomaly2MeanAnomaly(nu_1, ecc)
        M_2 = trueAnomaly2MeanAnomaly(nu_2, ecc)
        deltaM = M_2 - M_1

        # Compute time of flight (Python logic)
        if 0 <= ecc < 1
            tof_computed = sqrt(a^3 / μ) * deltaM
            if tof_computed < 0
                tof_computed += 2π * sqrt(a^3 / μ)
            end
        elseif ecc == 1
            tof_computed = 0.5 * sqrt(p^3 / μ) * deltaM
            if tof_computed < 0
                return Inf
            end
        else
            tof_computed = sqrt(-a^3 / μ) * deltaM
        end

        # Return the error exactly as Python _f function
        t_ref = sqrt(r1_norm^3 / μ)
        tau_computed = tof_computed / t_ref
        tau_desired = tof / t_ref

        # Compute y values as in Python
        y_computed = log(tau_computed)
        y_desired = log(tau_desired)

        return y_computed - y_desired
    end

    # Use Roots.jl to find the solution (like Python scipy.optimize.newton)
    # Use Bisection method to avoid automatic differentiation issues
    x_sol = find_zero(f, (-5.0, 5.0), Roots.Bisection(), atol = atol, maxevals = maxiter)

    # Compute final solution using solved x (Python coe_at_eccT)
    ecc_T = if dtheta > π
        X = exp((1/ecc_H + 1/ecc_P) * x_sol)
        ecc_P * ecc_H * (X - 1) / (ecc_P + ecc_H * X)
    else
        ecc_P * (1 - exp(-x_sol / ecc_P))
    end
    ecc = sqrt(ecc_T^2 + ecc_F^2)
    p = p_F - ecc_T * r1_norm * r2_norm * sin(dtheta) / c_norm

    # Compute argument of periapsis and true anomalies
    y = ecc_F * sin(w_c) + ecc_T * cos(w_c)
    x_coord = ecc_F * cos(w_c) - ecc_T * sin(w_c)
    w = atan(y, x_coord)
    nu_1 = -w
    nu_2 = nu_1 + dtheta

    # Convert to velocity vectors using AstroCoords directly
    # Convert semi-latus rectum to semi-major axis: a = p / (1 - e²)
    a = p / (1 - ecc^2)

    # Create Keplerian coordinates and convert to Cartesian using AstroCoords
    kepler1 = Keplerian(a, ecc, 0.0, 0.0, w, nu_1)
    kepler2 = Keplerian(a, ecc, 0.0, 0.0, w, nu_2)

    cart1 = Cartesian(kepler1, μ)
    cart2 = Cartesian(kepler2, μ)

    v1 = [cart1.ẋ, cart1.ẏ, cart1.ż]
    v2 = [cart2.ẋ, cart2.ẏ, cart2.ż]

    return v1, v2, maxiter, :SUCCESS
end

# Note: Full Avanzini implementation moved to utils.jl for reusability
