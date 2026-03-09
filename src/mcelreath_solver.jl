export mcelreath2025, McElreathSolver

"""
    McElreathSolver

Lambert's problem solver using McElreath, Down, and Majji's universal approach (2025).

Uses matrix exponential solutions to the Sundman transform, evaluating the transfer
time equation through exponential/trigonometric universal functions. Provides a unified
formulation across elliptic and hyperbolic regimes with no separate code paths for the
transfer time evaluation. Multi-revolution solutions are handled by proper bounding of
the universal variable, removing the explicit multi-revolution term from the transfer
time equation for improved numerical stability at ultra-high revolution counts.

# Fields
- `M::Int`: Number of revolutions (default: 0)
- `prograde::Bool`: Direction of motion (default: true)
- `low_path::Bool`: For multi-rev, select low-energy path (default: true)
- `maxiter::Int`: Maximum number of iterations (default: 25)
- `rtol::Float64`: Relative tolerance for convergence (default: 1e-10)

# References
[1] McElreath, J., Down, I. M., & Majji, M. (2025). A universal approach for solving
    the multi-revolution Lambert's problem. *Celestial Mechanics and Dynamical Astronomy*,
    137, 22. DOI: [10.1007/s10569-025-10251-5](https://doi.org/10.1007/s10569-025-10251-5)
"""
@with_kw struct McElreathSolver <: AbstractLambertSolver
    M::Int = 0
    prograde::Bool = true
    low_path::Bool = true
    maxiter::Int = 25
    rtol::Float64 = 1e-10
end

function SciMLBase.solve(problem::LambertProblem, solver::McElreathSolver)
    @unpack μ, r1, r2, tof = problem
    @unpack M, prograde, low_path, maxiter, rtol = solver

    v1, v2, numiter, retcode = mcelreath2025(
        μ,
        r1,
        r2,
        tof;
        M = M,
        prograde = prograde,
        low_path = low_path,
        maxiter = maxiter,
        rtol = rtol,
    )

    return LambertSolution(v1, v2, numiter, retcode)
end

# ============================================================================
# Universal functions via trigonometric evaluation (Appendix C of the paper)
# ============================================================================

"""
Evaluate the universal functions U0, U1, U2, U3, U0*, U1* from the universal
variable x and the inverse semi-major axis α.

For elliptic orbits (α > 0, x real): uses standard trig functions.
For hyperbolic orbits (α < 0, x_hyp real positive): uses hyperbolic trig functions.
"""
@inline function _mcelreath_universal_elliptic(x::Number, α::Number)
    sqα = sqrt(α)
    χ = x / sqα
    halfx = x / 2

    U0 = cos(x)
    U1 = sin(x) / sqα
    U2 = (1 - U0) / α
    U3 = (χ - U1) / α
    U0s = cos(halfx)
    U1s = sin(halfx) / sqα

    return U0, U1, U2, U3, U0s, U1s
end

@inline function _mcelreath_universal_hyperbolic(x_hyp::Number, α::Number)
    neg_α = -α
    sq_neg_α = sqrt(neg_α)
    χ = x_hyp / sq_neg_α
    halfx = x_hyp / 2

    U0 = cosh(x_hyp)
    U1 = sinh(x_hyp) / sq_neg_α
    U2 = (U0 - 1) / neg_α
    U3 = (U1 - χ) / neg_α
    U0s = cosh(halfx)
    U1s = sinh(halfx) / sq_neg_α

    return U0, U1, U2, U3, U0s, U1s
end

"""
Compute α from x using Eq. 37: α = (1 - U0) / (rf + r0 - 2m*U0*)
For elliptic: x is the universal variable directly.
For hyperbolic: x_hyp is the magnitude of the hyperbolic anomaly difference.
"""
@inline function _mcelreath_alpha_elliptic(x::Number, r0::Number, rf::Number, m::Number)
    U0 = cos(x)
    U0s = cos(x / 2)
    return (1 - U0) / (rf + r0 - 2 * m * U0s)
end

@inline function _mcelreath_alpha_hyperbolic(
    x_hyp::Number,
    r0::Number,
    rf::Number,
    m::Number,
)
    U0 = cosh(x_hyp)
    U0s = cosh(x_hyp / 2)
    return (1 - U0) / (rf + r0 - 2 * m * U0s)
end

# ============================================================================
# Transfer time and derivatives (Section 3.1 and Appendix A)
# ============================================================================

"""
Compute the transfer time T = 2*U1s*m + U3 and its first three derivatives
with respect to x for elliptic orbits.
"""
function _mcelreath_tof_elliptic(x::Number, r0::Number, rf::Number, m::Number, μ::Number)
    α = _mcelreath_alpha_elliptic(x, r0, rf, m)

    if abs(α) < 1e-14
        χ0 = sqrt(2 * (rf + r0 - 2 * m))
        T = m * χ0 + χ0^3 / 6
        return T, 0.0, 0.0, 0.0, α
    end

    U0, U1, U2, U3, U0s, U1s = _mcelreath_universal_elliptic(x, α)
    sqα = sqrt(α)

    T = 2 * U1s * m + U3

    # Derivatives of α (Appendix A, Eqs. A4-A6)
    α′ = sqα * (2 * U0s - m * α) / (2 * U1s)
    α′′ = -α / 2 + α′ * (U0s - 2 * m * α) / (2 * sqα * U1s)

    denom_a3 = 4 * U2
    numer_a3_t1 = (2 * m^2 * α - m * U0s) / denom_a3
    numer_a3_t2 = -α′ * (U0s + 2 * α * m) / (4 * α^(3 / 2) * U1s)
    numer_a3_t3 = -3 / 4
    α′′′ =
        α′ * (numer_a3_t1 + numer_a3_t2 + numer_a3_t3) +
        α′′ * (U0s - 2 * m * α) / (2 * sqα * U1s)

    # Derivatives of T (Appendix A, Eqs. A1-A3)
    T′ = ((m * U0s + U2) * sqα - (T / 2 + U3) * α′) / α
    T′′ =
        (
            U1 - m * α * U1s / 2 + α′ * (2 * m * U0s / sqα - 3 * T′) -
            3 * T * (α′)^2 / (4 * α) - (T / 2 + U3) * α′′
        ) / α
    T′′′ =
        (
            U0 / sqα - m * sqα * U0s / 2 -
            α′ *
            (9 * T′′ / 2 + 3 * m * U1s / 2 + α′ * (9 * T′ / 4 - 3 * T * α′ / (8 * α^2))) +
            α′′ *
            (2 * m * U0s / sqα - U2 / sqα - 7 * T′ / 2 + α′ * (U3 / α - 7 * T / (4 * α))) -
            α′′′ * (T / 2 + U3)
        ) / α

    return T, T′, T′′, T′′′, α
end

"""
Compute T and derivatives for hyperbolic orbits, where x_hyp = ΔH > 0.
"""
function _mcelreath_tof_hyperbolic(
    x_hyp::Number,
    r0::Number,
    rf::Number,
    m::Number,
    μ::Number,
)
    α = _mcelreath_alpha_hyperbolic(x_hyp, r0, rf, m)
    neg_α = -α
    sq_neg_α = sqrt(neg_α)

    U0, U1, U2, U3, U0s, U1s = _mcelreath_universal_hyperbolic(x_hyp, α)

    T = 2 * U1s * m + U3

    # For hyperbolic, derivatives w.r.t. x_hyp need a sign flip because
    # x = j*x_hyp, so d/dx = (1/j)*d/dx_hyp = -j*d/dx_hyp
    # We compute derivatives w.r.t. the real variable x_hyp and account
    # for the chain rule: the Householder update works on x_hyp directly.

    # α and derivatives w.r.t. x_hyp (chain rule: d/dx_hyp of cos(jx_hyp)=cosh(x_hyp) → sinh(x_hyp))
    # Re-derive for hyperbolic: U0 = cosh(xh), U0s = cosh(xh/2)
    # dU0/dxh = sinh(xh), dU0s/dxh = sinh(xh/2)/2

    sxh = sinh(x_hyp)
    sxh2 = sinh(x_hyp / 2)
    cxh = U0    # cosh(x_hyp)
    cxh2 = U0s  # cosh(x_hyp/2)

    denom = rf + r0 - 2 * m * cxh2
    α_h = (1 - cxh) / denom  # This is negative for hyperbolic

    # dα/dxh via quotient rule
    d_num = -sxh
    d_den = -m * sxh2
    α′_h = (d_num * denom - (1 - cxh) * d_den) / denom^2

    # dT/dxh: T = 2*m*sinh(xh/2)/sqrt(-α) + (xh/sqrt(-α) - sinh(xh)/sqrt(-α))/(-α)
    # This is complex to differentiate analytically. Use the chain rule approach.
    # Since the Householder update works on the real variable x_hyp, we need
    # dT/dxh, d²T/dxh², d³T/dxh³.

    # Numerical derivatives for robustness in the hyperbolic regime
    eps_h = max(abs(x_hyp) * 1e-7, 1e-12)

    Tp, _, _, _, _ = _mcelreath_tof_hyperbolic_val(x_hyp + eps_h, r0, rf, m)
    Tm_v, _, _, _, _ = _mcelreath_tof_hyperbolic_val(x_hyp - eps_h, r0, rf, m)
    Tp2, _, _, _, _ = _mcelreath_tof_hyperbolic_val(x_hyp + 2 * eps_h, r0, rf, m)
    Tm2, _, _, _, _ = _mcelreath_tof_hyperbolic_val(x_hyp - 2 * eps_h, r0, rf, m)

    T′ = (Tp - Tm_v) / (2 * eps_h)
    T′′ = (Tp - 2 * T + Tm_v) / (eps_h^2)
    T′′′ = (Tp2 - 2 * Tp + 2 * Tm_v - Tm2) / (2 * eps_h^3)

    return T, T′, T′′, T′′′, α
end

"""
Compute only the transfer time value for hyperbolic (no derivatives).
"""
@inline function _mcelreath_tof_hyperbolic_val(
    x_hyp::Number,
    r0::Number,
    rf::Number,
    m::Number,
)
    α = _mcelreath_alpha_hyperbolic(x_hyp, r0, rf, m)
    _, _, _, U3, _, U1s = _mcelreath_universal_hyperbolic(x_hyp, α)
    T = 2 * U1s * m + U3
    return T, 0.0, 0.0, 0.0, α
end

# ============================================================================
# Parabolic transfer time (Eq. 49)
# ============================================================================

@inline function _mcelreath_parabolic_time(r0::Number, rf::Number, m::Number)
    χ0 = sqrt(2 * (rf + r0 - 2 * m))
    return m * χ0 + χ0^3 / 6
end

# ============================================================================
# Initial guess schemes (Section 3.2)
# ============================================================================

"""
Single-revolution elliptic initial guess (Eq. 58).
"""
function _mcelreath_guess_elliptic_singlerev(
    T_tilde::Number,
    T0::Number,
    θ::Number,
    r0::Number,
    rf::Number,
    m::Number,
)
    xθ = θ
    Tθ = _mcelreath_eval_T_at_x(xθ, r0, rf, m)

    if abs(Tθ - T0) < 1e-15
        return xθ
    end

    β4 = (xθ / (2π - xθ)) * sqrt((T_tilde - T0) / (Tθ - T0))
    return 2π * β4 / (1 + β4)
end

"""
Maximum valid xh value where α transitions from negative to positive.
The denominator of α is rf+r0-2m*cosh(xh/2), which is zero when
cosh(xh/2) = (rf+r0)/(2|m|), giving the upper bound for the hyperbolic domain.
"""
@inline function _mcelreath_xh_max(r0::Number, rf::Number, m::Number)
    ratio = (r0 + rf) / (2 * abs(m))
    if ratio <= 1
        return 0.0
    end
    return 2 * acosh(ratio)
end

"""
Hyperbolic initial guess for θ ≥ π (Eq. 52).
"""
function _mcelreath_guess_hyperbolic_large_theta(
    T_tilde::Number,
    T0::Number,
    θ::Number,
    r0::Number,
    rf::Number,
    m::Number,
    xh_max::Number,
)
    xjθ = min(θ, xh_max * 0.8)
    Tjθ = _mcelreath_eval_T_hyp_at_xh(xjθ, r0, rf, m)

    log_ratio_star = log(T_tilde / T0)
    log_ratio_jθ = log(Tjθ / T0)

    if abs(log_ratio_jθ) < 1e-15
        return xjθ
    end

    xh0 = xjθ * sqrt(log_ratio_star / log_ratio_jθ)
    return min(xh0, xh_max * 0.99)
end

"""
Hyperbolic initial guess for θ < π (Eqs. 53-56).
τ = (r0+rf)/(2m) corresponds to the rectilinear limit where α→-∞.
"""
function _mcelreath_guess_hyperbolic_small_theta(
    T_tilde::Number,
    T0::Number,
    r0::Number,
    rf::Number,
    m::Number,
    xh_max::Number,
)
    τ_val = (r0 + rf) / (2 * abs(m))
    if τ_val <= 1
        return xh_max * 0.5
    end
    xτ = log(τ_val + sqrt(τ_val^2 - 1))

    xτ_eval = min(xτ * exp(-2.5), xh_max * 0.8)
    T_at_eval = _mcelreath_eval_T_hyp_at_xh(xτ_eval, r0, rf, m)

    log_num = -T0^(8 / 5)
    log_den = T_at_eval^(8 / 5) - T0^(8 / 5)
    if abs(log_den) < 1e-15
        return xτ * 0.5
    end
    β2 = (2 / 5) * log(log_num / log_den)

    if abs(β2) < 1e-15
        return xτ * 0.5
    end

    ratio = T_tilde / T0
    if ratio >= 1
        return xτ * 0.01
    end
    xh0 = xτ * (1 - ratio)^(1 / β2)
    return min(xh0, xh_max * 0.99)
end

"""
Multi-revolution initial guess for minimum time xm (Eq. 59).
"""
function _mcelreath_guess_xm(θ::Number, N::Int)
    sqN = sqrt(N)
    if θ < π
        return 2 * (θ + (1 / 6)^(2.5 / sqN)^(1 / 2.5) - (1 / 3)^sqN) + 2π * N
    else
        return -2 * (θ + 2π + (1 / 20)^sqN^(1 / 2.5) - (1 / 4)^sqN) + 2π * (N + 1)
    end
end

"""
Multi-revolution initial guesses for the two solutions (Eqs. 60-64).
"""
function _mcelreath_guess_multirev(
    T_tilde::Number,
    Tm::Number,
    xm::Number,
    θ::Number,
    N::Int,
    r0::Number,
    rf::Number,
    m_geom::Number,
    low_path::Bool,
)
    xθ = 2π * N + θ
    Tθ = _mcelreath_eval_T_at_x(xθ, r0, rf, m_geom)

    β5_num = abs(Tθ - Tm) * (2π * N - xθ)^2 * (2π * (N + 1) - xθ)^2
    β5_den = (xm - xθ)^2
    if abs(β5_den) < 1e-30
        β5 = 1e30
    else
        β5 = β5_num / β5_den
    end

    diff = T_tilde - Tm
    if diff < 0
        return NaN
    end

    γ = low_path ? sqrt(diff / β5) : -sqrt(diff / β5)

    c1 = γ
    c2 = 1 - 2π * (2 * N + 1) * γ
    c3 = 4π^2 * N * (N + 1) * γ - xm

    disc = c2^2 - 4 * c1 * c3
    if disc < 0
        return xm
    end

    return (-c2 + sqrt(disc)) / (2 * c1)
end

"""
Evaluate scaled transfer time at a given x (for elliptic initial guess computation).
"""
function _mcelreath_eval_T_at_x(x::Number, r0::Number, rf::Number, m::Number)
    α = _mcelreath_alpha_elliptic(x, r0, rf, m)

    if abs(α) < 1e-14
        χ0 = sqrt(2 * (rf + r0 - 2 * m))
        return m * χ0 + χ0^3 / 6
    end

    _, _, _, U3, _, U1s = _mcelreath_universal_elliptic(x, α)
    return 2 * U1s * m + U3
end

"""
Evaluate scaled transfer time at a given x_hyp (for hyperbolic initial guess computation).
"""
function _mcelreath_eval_T_hyp_at_xh(x_hyp::Number, r0::Number, rf::Number, m::Number)
    α = _mcelreath_alpha_hyperbolic(x_hyp, r0, rf, m)
    _, _, _, U3, _, U1s = _mcelreath_universal_hyperbolic(x_hyp, α)
    return 2 * U1s * m + U3
end

# ============================================================================
# Root solvers
# ============================================================================

"""
Find xm (minimum transfer time) using Halley's method (Eq. 66).
"""
function _mcelreath_find_xm(
    xm0::Number,
    r0::Number,
    rf::Number,
    m::Number,
    μ::Number,
    T_tilde::Number,
    maxiter::Int,
)
    xm = xm0

    for iter = 1:maxiter
        _, T′, T′′, T′′′, _ = _mcelreath_tof_elliptic(xm, r0, rf, m, μ)

        g = T′
        g′ = T′′
        g′′ = T′′′

        denom = 2 * g′^2 - g * g′′
        if abs(denom) < 1e-30
            break
        end

        dx = -2 * g * g′ / denom
        xm_new = xm + dx

        # Adaptive tolerance (Eq. 76)
        Tm_curr = _mcelreath_eval_T_at_x(xm, r0, rf, m)
        ratio = abs(1 - Tm_curr / T_tilde)
        ϵ = ratio > 1e-4 ? 1e-2 : 1e-8

        ref = xm_new - 2π * max(0, round(Int, xm_new / (2π)) - 1)
        if abs(ref) < 1e-30
            ref = 1.0
        end

        if abs(dx / ref) < ϵ
            xm = xm_new
            break
        end
        xm = xm_new
    end

    return xm
end

"""
3rd-order Householder iteration for the main solution (Eq. 70).
Returns (x_converged, numiter, converged).
"""
function _mcelreath_householder_elliptic(
    x0::Number,
    T_tilde::Number,
    r0::Number,
    rf::Number,
    m::Number,
    μ::Number,
    N::Int,
    maxiter::Int,
    rtol::Float64,
    x_lo::Number,
    x_hi::Number,
)
    x = x0
    numiter = 0

    for iter = 1:maxiter
        numiter = iter

        T, T′, T′′, T′′′, _ = _mcelreath_tof_elliptic(x, r0, rf, m, μ)
        h = T - T_tilde
        h′ = T′
        h′′ = T′′
        h′′′ = T′′′

        if abs(h′) < 1e-30
            break
        end

        # 3rd-order Householder (Eq. 70)
        num = 6 * h * h′^2 - 3 * h′′^2 * h  # simplified: keep full form
        num = 6 * h * (h′^2 - h * h′′ / 2)
        den = 6 * h′^3 - 6 * h * h′ * h′′ + h^2 * h′′′

        if abs(den) < 1e-30
            # Fall back to Newton
            dx = -h / h′
        else
            dx = -num / den
        end

        x_new = x + dx

        # Bound checking - fall back to lower order if out of bounds
        if x_new <= x_lo || x_new >= x_hi
            # Try 2nd-order Householder (Halley)
            den2 = 2 * h′^2 - h * h′′
            if abs(den2) > 1e-30
                dx = -2 * h * h′ / den2
                x_new = x + dx
            end
            if x_new <= x_lo || x_new >= x_hi
                # Fall back to Newton
                dx = -h / h′
                x_new = x + dx
            end
            if x_new <= x_lo || x_new >= x_hi
                # Bisect toward center
                x_new = (x_lo + x_hi) / 2
            end
        end

        ref = x_new - 2π * N
        if abs(ref) < 1e-30
            ref = 1.0
        end

        if abs(dx / ref) < rtol
            return x_new, numiter, true
        end

        x = x_new
    end

    return x, numiter, false
end

"""
Bracket the hyperbolic root using bisection within the valid domain [0, xh_max).
T(xh) is monotonically decreasing from T0 at xh=0 toward α→-∞ near xh_max.
"""
function _mcelreath_bracket_hyperbolic(
    T_tilde::Number,
    r0::Number,
    rf::Number,
    m::Number,
    xh_max::Number,
    maxbisect::Int = 50,
)
    xh_lo = 1e-10
    xh_hi = xh_max * (1 - 1e-8)

    for _ = 1:maxbisect
        xh_mid = (xh_lo + xh_hi) / 2
        T_mid = _mcelreath_eval_T_hyp_at_xh(xh_mid, r0, rf, m)
        if T_mid > T_tilde
            xh_lo = xh_mid
        else
            xh_hi = xh_mid
        end
        if (xh_hi - xh_lo) / max(xh_hi, 1.0) < 1e-10
            break
        end
    end

    return (xh_lo + xh_hi) / 2
end

"""
Householder iteration for hyperbolic solutions (x_hyp > 0), bounded in [0, xh_max).
"""
function _mcelreath_householder_hyperbolic(
    xh0::Number,
    T_tilde::Number,
    r0::Number,
    rf::Number,
    m::Number,
    μ::Number,
    maxiter::Int,
    rtol::Float64,
    xh_max::Number,
)
    xh = xh0
    numiter = 0
    xh_lo = 1e-14
    xh_hi = xh_max * (1 - 1e-8)

    for iter = 1:maxiter
        numiter = iter

        T, T′, T′′, T′′′, _ = _mcelreath_tof_hyperbolic(xh, r0, rf, m, μ)
        h = T - T_tilde
        h′ = T′

        if abs(h′) < 1e-30
            break
        end

        # 3rd-order Householder
        num = 6 * h * (h′^2 - h * T′′ / 2)
        den = 6 * h′^3 - 6 * h * h′ * T′′ + h^2 * T′′′

        if abs(den) < 1e-30
            dx = -h / h′
        else
            dx = -num / den
        end

        xh_new = xh + dx

        # Bound enforcement with fallback to lower-order methods
        if xh_new <= xh_lo || xh_new >= xh_hi
            den2 = 2 * h′^2 - h * T′′
            if abs(den2) > 1e-30
                dx = -2 * h * h′ / den2
                xh_new = xh + dx
            end
            if xh_new <= xh_lo || xh_new >= xh_hi
                dx = -h / h′
                xh_new = xh + dx
            end
            if xh_new <= xh_lo || xh_new >= xh_hi
                xh_new = (xh_lo + xh_hi) / 2
            end
        end

        if abs(xh) > 1e-30 && abs(dx / xh) < rtol
            return xh_new, numiter, true
        end

        xh = xh_new
    end

    return xh, numiter, false
end

# ============================================================================
# Velocity reconstruction (Stark transformation, Eqs. 72-75)
# ============================================================================

function _mcelreath_velocity(
    x::Number,
    is_hyp::Bool,
    r0_vec::AbstractVector{<:Number},
    rf_vec::AbstractVector{<:Number},
    r0::Number,
    rf::Number,
    θ::Number,
    m::Number,
    μ::Number,
)
    if is_hyp
        α = _mcelreath_alpha_hyperbolic(x, r0, rf, m)
        U0s = cosh(x / 2)
    else
        α = _mcelreath_alpha_elliptic(x, r0, rf, m)
        U0s = cos(x / 2)
    end

    sinθ = sin(θ)
    cosθ = cos(θ)

    p_semi = r0 * rf * (1 - cosθ) / (r0 + rf - 2 * m * U0s)

    r0s = SVector{3}(r0_vec[1], r0_vec[2], r0_vec[3])
    rfs = SVector{3}(rf_vec[1], rf_vec[2], rf_vec[3])

    c_vec = rfs - r0s
    c_mag = norm(c_vec)
    ic = c_vec / c_mag
    ir0 = r0s / r0
    irf = rfs / rf

    vc = c_mag * sqrt(μ * p_semi) / (r0 * rf * sinθ)
    vρ = sqrt(μ / p_semi) * (1 - cosθ) / sinθ

    v1 = vc * ic + vρ * ir0
    v2 = vc * ic - vρ * irf

    return v1, v2
end

# ============================================================================
# Main solver
# ============================================================================

"""
    mcelreath2025(μ, r1, r2, tof; M=0, prograde=true, low_path=true, maxiter=25, rtol=1e-10)

Solve Lambert's problem using McElreath, Down, and Majji's universal approach
with Sundman transform-based universal functions and Householder iteration.

# Arguments
- `μ`: Gravitational parameter [L³/T²]
- `r1`: Initial position vector [L]
- `r2`: Final position vector [L]
- `tof`: Time of flight [T]
- `M`: Number of complete revolutions (default: 0)
- `prograde`: Direction of motion (default: true)
- `low_path`: For multi-rev, select low-energy path (default: true)
- `maxiter`: Maximum iterations (default: 25)
- `rtol`: Relative tolerance (default: 1e-10)

# Returns
- `v1`: Departure velocity vector [L/T]
- `v2`: Arrival velocity vector [L/T]
- `numiter`: Number of iterations used
- `retcode`: `:SUCCESS`, `:MAXIMUM_ITERATIONS`, or `:NO_SOLUTION`

# References
[1] McElreath, Down, Majji (2025), Celestial Mechanics and Dynamical Astronomy, 137:22
"""
function mcelreath2025(
    μ::Number,
    r1::AbstractVector{<:Number},
    r2::AbstractVector{<:Number},
    tof::Number;
    M::Int = 0,
    prograde::Bool = true,
    low_path::Bool = true,
    maxiter::Int = 25,
    rtol::Float64 = 1e-10,
)
    assert_parameters_are_valid(μ, r1, r2, tof, M)

    r0_norm = norm(r1)
    rf_norm = norm(r2)
    _, _, _, dtheta = lambert_geometry(r1, r2, prograde)
    assert_transfer_angle_not_zero(dtheta)

    # Transfer angle for this formulation (0, 2π)
    θ = dtheta

    # Geometric constant m (Eq. 41) — sign alternates with N for multi-rev
    m_geom = (-1)^M * sqrt(r0_norm * rf_norm) * cos(θ / 2)

    # Scaled transfer time
    T_tilde = sqrt(μ) * tof

    N = M  # revolution count

    if N == 0
        # Parabolic time (Eq. 49)
        T0 = _mcelreath_parabolic_time(r0_norm, rf_norm, m_geom)

        if abs(T_tilde - T0) < 1e-12 * T0
            # Parabolic solution — reconstruct directly
            v1, v2 = _mcelreath_velocity(0.0, false, r1, r2, r0_norm, rf_norm, θ, m_geom, μ)
            return v1, v2, 0, :SUCCESS
        elseif T_tilde > T0
            # Elliptic single-rev
            x0 = _mcelreath_guess_elliptic_singlerev(
                T_tilde,
                T0,
                θ,
                r0_norm,
                rf_norm,
                m_geom,
            )
            x0 = clamp(x0, 1e-10, 2π - 1e-10)

            x, numiter, converged = _mcelreath_householder_elliptic(
                x0,
                T_tilde,
                r0_norm,
                rf_norm,
                m_geom,
                μ,
                0,
                maxiter,
                rtol,
                0.0,
                2π,
            )

            retcode = converged ? :SUCCESS : handle_max_iterations(numiter, maxiter)
            if retcode != :SUCCESS
                return zero(SVector{3,Float64}), zero(SVector{3,Float64}), numiter, retcode
            end

            v1, v2 = _mcelreath_velocity(x, false, r1, r2, r0_norm, rf_norm, θ, m_geom, μ)
            return v1, v2, numiter, :SUCCESS
        else
            # Hyperbolic — compute valid domain bound and use robust initial guess
            xh_max = _mcelreath_xh_max(r0_norm, rf_norm, m_geom)

            if xh_max < 1e-12
                return zero(SVector{3,Float64}), zero(SVector{3,Float64}), 0, :NO_SOLUTION
            end

            # Try paper's initial guess first, fall back to bracketing
            if θ >= π
                xh0 = _mcelreath_guess_hyperbolic_large_theta(
                    T_tilde,
                    T0,
                    θ,
                    r0_norm,
                    rf_norm,
                    m_geom,
                    xh_max,
                )
            else
                xh0 = _mcelreath_guess_hyperbolic_small_theta(
                    T_tilde,
                    T0,
                    r0_norm,
                    rf_norm,
                    m_geom,
                    xh_max,
                )
            end

            # Validate initial guess; fall back to bisection if it produces invalid values
            if isnan(xh0) || xh0 <= 0 || xh0 >= xh_max
                xh0 =
                    _mcelreath_bracket_hyperbolic(T_tilde, r0_norm, rf_norm, m_geom, xh_max)
            end

            xh0 = clamp(xh0, 1e-10, xh_max * (1 - 1e-8))

            xh, numiter, converged = _mcelreath_householder_hyperbolic(
                xh0,
                T_tilde,
                r0_norm,
                rf_norm,
                m_geom,
                μ,
                maxiter,
                rtol,
                xh_max,
            )

            retcode = converged ? :SUCCESS : handle_max_iterations(numiter, maxiter)
            if retcode != :SUCCESS
                return zero(SVector{3,Float64}), zero(SVector{3,Float64}), numiter, retcode
            end

            v1, v2 = _mcelreath_velocity(xh, true, r1, r2, r0_norm, rf_norm, θ, m_geom, μ)
            return v1, v2, numiter, :SUCCESS
        end
    else
        # Multi-revolution (N ≥ 1)
        xm0 = _mcelreath_guess_xm(θ, N)
        xm0 = clamp(xm0, 2π * N + 1e-6, 2π * (N + 1) - 1e-6)

        xm = _mcelreath_find_xm(xm0, r0_norm, rf_norm, m_geom, μ, T_tilde, maxiter)
        xm = clamp(xm, 2π * N + 1e-10, 2π * (N + 1) - 1e-10)

        Tm = _mcelreath_eval_T_at_x(xm, r0_norm, rf_norm, m_geom)

        if T_tilde < Tm
            return zero(SVector{3,Float64}), zero(SVector{3,Float64}), 0, :NO_SOLUTION
        end

        if abs(T_tilde - Tm) < 1e-12 * Tm
            v1, v2 = _mcelreath_velocity(xm, false, r1, r2, r0_norm, rf_norm, θ, m_geom, μ)
            return v1, v2, 0, :SUCCESS
        end

        x0 = _mcelreath_guess_multirev(
            T_tilde,
            Tm,
            xm,
            θ,
            N,
            r0_norm,
            rf_norm,
            m_geom,
            low_path,
        )

        if isnan(x0)
            return zero(SVector{3,Float64}), zero(SVector{3,Float64}), 0, :NO_SOLUTION
        end

        x0 = clamp(x0, 2π * N + 1e-10, 2π * (N + 1) - 1e-10)

        x, numiter, converged = _mcelreath_householder_elliptic(
            x0,
            T_tilde,
            r0_norm,
            rf_norm,
            m_geom,
            μ,
            N,
            maxiter,
            rtol,
            2π * N,
            2π * (N + 1),
        )

        retcode = converged ? :SUCCESS : handle_max_iterations(numiter, maxiter)
        if retcode != :SUCCESS
            return zero(SVector{3,Float64}), zero(SVector{3,Float64}), numiter, retcode
        end

        v1, v2 = _mcelreath_velocity(x, false, r1, r2, r0_norm, rf_norm, θ, m_geom, μ)
        return v1, v2, numiter, :SUCCESS
    end
end
