export gooding1990, GoodingSolver

"""
    GoodingSolver

Lambert's problem solver using the method proposed by R. H. Gooding in 1990.

# Fields
- `M::Int`: Number of revolutions (default: 0)
- `prograde::Bool`: Direction of motion - true for prograde, false for retrograde (default: true)
- `low_path::Bool`: Path selection when multiple solutions exist (default: true)
- `maxiter::Int`: Maximum number of iterations (default: 35)
- `atol::Float64`: Absolute tolerance (default: 1e-5)
- `rtol::Float64`: Relative tolerance (default: 1e-7)

# References
Gooding, R. H. (1990). A procedure for the solution of Lambert's orbital
boundary-value problem. Celestial Mechanics and Dynamical Astronomy, 48(2), 145-165.
"""
@with_kw struct GoodingSolver <: AbstractLambertSolver
    M::Int = 0
    prograde::Bool = true
    low_path::Bool = true
    maxiter::Int = 35
    atol::Float64 = 1e-5
    rtol::Float64 = 1e-7
end

function SciMLBase.solve(problem::LambertProblem, solver::GoodingSolver)
    @unpack μ, r1, r2, tof = problem
    @unpack M, prograde, low_path, maxiter, atol, rtol = solver

    # Call the direct algorithm function
    v1, v2, numiter, retcode = gooding1990(
        μ,
        r1,
        r2,
        tof;
        M = M,
        prograde = prograde,
        low_path = low_path,
        maxiter = maxiter,
        atol = atol,
        rtol = rtol,
    )

    return LambertSolution(v1, v2, numiter, retcode)
end

"""
    gooding1990(μ, r1, r2, tof, M=0, prograde=true, low_path=true, maxiter=35, atol=1e-5, rtol=1e-7)

Lambert's problem solver using the method proposed by R. H. Gooding in 1990.

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
- `numiter`: Number of iterations
- `retcode`: Return code

# References
Gooding, R. H. (1990). A procedure for the solution of Lambert's orbital
boundary-value problem. Celestial Mechanics and Dynamical Astronomy, 48(2), 145-165.
"""
function gooding1990(
    μ::Number,
    r1::Vector{<:Number},
    r2::Vector{<:Number},
    tof::Number;
    M::Int = 0,
    prograde::Bool = true,
    low_path::Bool = true,
    maxiter::Int = 35,
    atol::Float64 = 1e-5,
    rtol::Float64 = 1e-7,
)
    #TODO: Multi-revolution extension requires complex modifications:
    # 1. Extending tlamb function to handle M > 0 properly in series computations
    # 2. Modifying initial guess in xlamb for multi-rev cases (currently only handles M=0)
    # 3. Adding proper handling of multiple solutions (2M + 1 solutions exist)
    # 4. Implementing solution selection logic for high/low path in multi-rev cases  
    # 5. Extending vlamb function to compute correct velocity components
    # Reference: Gooding (1990) focused on single-rev, extensions needed for multi-rev
    # See recent work on universal approaches and Izzo's multi-rev handling for guidance
    (M > 0) && error(
        "Multi-revolution case not implemented - requires significant algorithm extension",
    )

    # Check that input parameters are safe
    assert_parameters_are_valid(μ, r1, r2, tof, M)

    # Compute basic geometry
    r1_norm, r2_norm, c_norm, θ = lambert_geometry(r1, r2, prograde)
    i_r1, i_r2, i_h_unnorm = compute_unit_vectors(r1, r2, r1_norm, r2_norm)
    assert_transfer_angle_not_zero(θ)

    # Include additional revolutions if necessary
    dtheta = 2.0π * M + θ

    # Compute a vector normal to orbit plane
    i_h = get_orbit_normal_vector(r1, r2, prograde)

    # Compute the tangential unitary vectors
    i_t1, i_t2 = compute_tangential_unit_vectors(i_h, i_r1, i_r2)

    # Compute velocity components
    vr1, vt1, vr2, vt2, numiter =
        vlamb(μ, r1_norm, r2_norm, dtheta, tof, low_path, maxiter, atol, rtol)

    retcode = handle_max_iterations(numiter, maxiter)

    if retcode != :SUCCESS
        return nothing, nothing, numiter, retcode
    end

    # Final velocity vectors
    v1, v2 =
        reconstruct_velocities_from_components(vr1, vt1, vr2, vt2, i_r1, i_t1, i_r2, i_t2)

    return v1, v2, numiter, :SUCCESS
end

"""
    tlamb(m, q, qsqfm1, x, n)

Auxiliary routine for computing the non-dimensional time of flight.
"""
function tlamb(m::Int, q::Number, qsqfm1::Number, x::Number, n::Int)
    sw = 0.4
    lm1 = n == -1
    l1 = n >= 1
    l2 = n >= 2
    l3 = n == 3
    qsq = q * q
    xsq = x * x
    u = (1.0 - x) * (1.0 + x)

    if !lm1
        dt, d2t, d3t = 0.0, 0.0, 0.0
    end

    if lm1 || m > 0 || x < 0.0 || abs(u) > sw
        # Direct computation
        y = √(abs(u))
        z = √(qsqfm1 + qsq * xsq)
        qx = q * x

        if qx <= 0.0
            a = z - qx
            b = q * z - x
        end

        if qx <= 0.0 && lm1
            aa = qsqfm1 / a
            bb = qsqfm1 * (qsq * u - xsq) / b
        end

        if (qx == 0.0 && lm1) || (qx > 0.0)
            aa = z + qx
            bb = q * z + x
        end

        if qx > 0.0
            a = qsqfm1 / aa
            b = qsqfm1 * (qsq * u - xsq) / bb
        end

        if lm1
            t, dt, d2t, d3t = 0, b, bb, aa
        else
            if qx * u >= 0.0
                g = x * z + q * u
            else
                g = (xsq - qsq * u) / (x * z - q * u)
            end

            f = a * y

            if x <= 1.0
                t = m * π + atan(f, g)
            else
                if f > sw
                    t = log(f + g)
                else
                    fg1 = f / (g + 1.0)
                    term = 2.0 * fg1
                    fg1sq = fg1 * fg1
                    t = term
                    twoi1 = 1.0

                    # Series computation
                    while true
                        twoi1 = twoi1 + 2.0
                        term = term * fg1sq
                        told = t
                        t = t + term / twoi1
                        if t == told
                            break
                        end
                    end
                end
            end

            t = 2.0 * (t / y + b) / u

            if l1 && z != 0.0
                qz = q / z
                qz2 = qz * qz
                qz = qz * qz2
                dt = (3.0 * x * t - 4.0 * (a + qx * qsqfm1) / z) / u
                if l2
                    d2t = (3.0 * t + 5.0 * x * dt + 4.0 * qz * qsqfm1) / u
                end
                if l3
                    d3t = (8.0 * dt + 7.0 * x * d2t - 12.0 * qz * qz2 * x * qsqfm1) / u
                end
            end
        end
    else
        # Compute by series
        u0i = 1.0

        if l1
            ;
            u1i = 1.0;
        end
        if l2
            ;
            u2i = 1.0;
        end
        if l3
            ;
            u3i = 1.0;
        end

        term = 4.0
        tq = q * qsqfm1
        i = 0

        if q < 0.5
            tqsum = 1.0 - q * qsq
        else
            tqsum = (1.0 / (1.0 + q) + q) * qsqfm1
        end

        ttmold = term / 3.0
        t = ttmold * tqsum

        # Series loop
        told = t
        while i < n || t != told
            i = i + 1
            p = i
            u0i = u0i * u

            if l1 && i > 1
                ;
                u1i = u1i * u;
            end
            if l2 && i > 2
                ;
                u2i = u2i * u;
            end
            if l3 && i > 3
                ;
                u3i = u3i * u;
            end

            term = term * (p - 0.5) / p
            tq = tq * qsq
            tqsum = tqsum + tq
            told = t
            tterm = term / (2.0 * p + 3.0)
            tqterm = tterm * tqsum
            t = t - u0i * ((1.5 * p + 0.25) * tqterm / (p * p - 0.25) - ttmold * tq)
            ttmold = tterm
            tqterm = tqterm * p

            if l1
                ;
                dt = dt + tqterm * u1i;
            end
            if l2
                ;
                d2t = d2t + tqterm * u2i * (p - 1.0);
            end
            if l3
                ;
                d3t = d3t + tqterm * u3i * (p - 1.0) * (p - 2.0);
            end
        end

        if l3
            ;
            d3t = 8.0 * x * (1.5 * d2t - xsq * d3t);
        end
        if l2
            ;
            d2t = 2.0 * (2.0 * xsq * d2t - dt);
        end
        if l1
            ;
            dt = -2.0 * x * dt;
        end

        t = t / xsq
    end

    return t, dt, d2t, d3t
end

"""
    vlamb(μ, r1_norm, r2_norm, dtheta, tof, low_path, maxiter, atol, rtol)

Auxiliary routine for computing the velocity vector components.
"""
function vlamb(
    μ::Number,
    r1_norm::Number,
    r2_norm::Number,
    dtheta::Number,
    tof::Number,
    low_path::Bool,
    maxiter::Int,
    atol::Float64,
    rtol::Float64,
)
    # The following yields m = 0 when th = 2π exactly
    thr2 = dtheta
    m = 0

    while thr2 > 2π
        thr2 = thr2 - 2π
        m = m + 1
    end
    thr2 = thr2 / 2.0

    # Compute auxiliary parameters
    dr = r1_norm - r2_norm
    r1r2 = r1_norm * r2_norm
    r1r2th = 4.0 * r1r2 * sin(thr2)^2
    csq = dr * dr + r1r2th
    c = √csq
    s = (r1_norm + r2_norm + c) / 2.0
    mus = √(μ * s / 2.0)
    qsqfm1 = c / s
    q = √r1r2 * cos(thr2) / s

    if c != 0.0
        ρ = dr / c
        σ = r1r2th / csq
    else
        ρ = 0.0
        σ = 1.0
    end

    t = 4.0 * mus * tof / s^2

    # Compute the number of solutions
    n_sol, x1_sol, x2_sol, numiter = xlamb(m, q, qsqfm1, t, maxiter, atol, rtol)

    # Filter the solution
    if n_sol > 1
        if low_path
            x_sol = max(x1_sol, x2_sol)
        else
            x_sol = min(x1_sol, x2_sol)
        end
    else
        x_sol = x1_sol
    end

    # Compute radial and tangential velocity components
    _, qzminx, qzplx, zplqx = tlamb(m, q, qsqfm1, x_sol, -1)
    vt2 = mus * zplqx * √σ
    vr1 = mus * (qzminx - qzplx * ρ) / r1_norm
    vt1 = vt2 / r1_norm
    vr2 = -mus * (qzminx + qzplx * ρ) / r2_norm
    vt2 = vt2 / r2_norm

    return vr1, vt1, vr2, vt2, numiter
end

"""
    xlamb(m, q, qsqfm1, tin, maxiter, atol, rtol)

Auxiliary routine for finding the independent variable.
"""
function xlamb(
    m::Int,
    q::Number,
    qsqfm1::Number,
    tin::Number,
    maxiter::Int,
    atol::Float64,
    rtol::Float64,
)
    # Declare auxiliary parameters
    xpl = 0
    c0 = 1.7
    c1 = 0.5
    c2 = 0.03
    thr2 = atan(qsqfm1, 2.0 * q) / π

    # Auxiliary function for 8th root
    d8rt(x) = x^(1/8)

    if m == 0
        # Single-rev starter
        n = 1

        # Call TLAMB routine
        t0, dt, d2t, d3t = tlamb(m, q, qsqfm1, 0.0, 0)
        tdiff = tin - t0

        if tdiff <= 0.0
            x = t0 * tdiff / (-4.0 * tin)
        else
            x = -tdiff / (tdiff + 4.0)
            w = x + c0 * √(2.0 * (1.0 - thr2))

            if w < 0.0
                x = x - √(d8rt(-w)) * (x + √(tdiff / (tdiff + 1.5 * t0)))
            end

            w = 4.0 / (4.0 + tdiff)
            x = x * (1.0 + x * (c1 * w - c2 * x * √w))
        end
    end

    # Newton-Raphson iterations
    numiter = 0

    for iter = 1:maxiter
        numiter = iter
        t, dt, d2t, d3t = tlamb(m, q, qsqfm1, x, 2)
        t = tin - t

        if dt != 0.0
            xold = x
            x = x + t * dt / (dt * dt + t * d2t / 2.0)

            # Check convergence
            x_atol, x_rtol = abs(x - xold), abs(x / xold - 1)
            if x_atol <= atol && x_rtol <= rtol
                break
            end
        end
    end

    return n, x, xpl, numiter
end
