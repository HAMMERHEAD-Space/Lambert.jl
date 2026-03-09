export russell2021, RussellSolver

"""
    RussellSolver

Lambert's problem solver using Russell's vercosine formulation (2019, 2021).

This algorithm uses the vercosine Lambert formulation with interpolated initial
guesses and up to third-order iteration corrections for robust, high-precision
convergence across all orbit types (elliptic, parabolic, hyperbolic) and
revolution counts.

# Fields
- `M::Int`: Number of revolutions (default: 0)
- `prograde::Bool`: Direction of motion (default: true)
- `low_path::Bool`: For multi-rev, select low-energy path (default: true)
- `maxiter::Int`: Maximum number of iterations (default: 50)
- `rtol::Float64`: Relative tolerance for convergence (default: 1e-14)
- `order::Int`: Correction order: 1=Newton, 2=Halley, 3=third-order (default: 3)

# References
[1] Russell, R. P. (2019). On the Solution to Every Lambert Problem. *Celestial
    Mechanics and Dynamical Astronomy*, 131, Article 50.
    DOI: [10.1007/s10569-019-9927-z](https://doi.org/10.1007/s10569-019-9927-z)

[2] Russell, R. P. (2021). Complete Lambert Solver Including Second-Order
    Sensitivities. *Journal of Guidance, Control, and Dynamics*, 45(2).
    DOI: [10.2514/1.G006089](https://doi.org/10.2514/1.G006089)
"""
@with_kw struct RussellSolver <: AbstractLambertSolver
    M::Int = 0
    prograde::Bool = true
    low_path::Bool = true
    maxiter::Int = 50
    rtol::Float64 = 1e-14
    order::Int = 3
end

function SciMLBase.solve(problem::LambertProblem, solver::RussellSolver)
    @unpack μ, r1, r2, tof = problem
    @unpack M, prograde, low_path, maxiter, rtol, order = solver

    v1, v2, numiter, retcode = russell2021(
        μ,
        r1,
        r2,
        tof;
        M = M,
        prograde = prograde,
        low_path = low_path,
        maxiter = maxiter,
        rtol = rtol,
        order = order,
    )

    return LambertSolution(v1, v2, numiter, retcode)
end

# ============================================================================
# Thresholds from Russell's ivLam implementation (tuned for Float64)
# ============================================================================
const _RUSSELL_SQRT2 = sqrt(2.0)
const _RUSSELL_ZERO_BAND = 0.02
const _RUSSELL_PARABOLA_BAND = 0.02
const _RUSSELL_BIG_K = 1000.0
const _RUSSELL_LITTLE_P = 0.1
const _RUSSELL_HUGE_TOF = 1e4
const _RUSSELL_ALT_TAU_THRESH = 1e-2
const _RUSSELL_K_CLOSE_MSQRT2 = 1e-2

const _RUSSELL_MAX_K_STEP_MREV = 2e-4
const _RUSSELL_SERIES_CONVERGE = 1.01
const _RUSSELL_SERIES_TAMP = 0.75
const _RUSSELL_SERIES_TAMP_THRESH = 2e-7
const _RUSSELL_K_BUMP_WEIGHT = 0.925
const _RUSSELL_K_MARGIN_EPS = 50.0 * eps(Float64)
const _RUSSELL_K_BUMP_INIT = 7e-7
const _RUSSELL_DELTA_VAR_TOL = 0.5

# ============================================================================
# W(k) function and derivatives — the core of the vercosine formulation
# ============================================================================

"""
    _russell_W!(dW, k, M, order)

Compute the vercosine W(k) function and its derivatives up to the specified
`order` (1-3). Results stored in `dW[1:order+1]` where `dW[1] = W(k)`.

The W(k) function plays the role of Stumpff functions in the vercosine
formulation. Four evaluation regions handle numerical precision:
- Region 0: k near zero (Taylor series)
- Region 1: Normal ellipse (closed-form with acos)
- Region 2: k near √2 — parabola band (Taylor series in ν = k - √2)
- Region 3: Hyperbola (closed-form with log)
"""
@inline function _russell_W(k::T, M::Int, order::Int) where {T<:Number}
    k² = k * k
    ν = k - _RUSSELL_SQRT2
    is_zero_rev = (M == 0)
    twoπN = 2π * M

    if k² <= _RUSSELL_ZERO_BAND^2
        kregion = 0
    elseif is_zero_rev && abs(ν) < _RUSSELL_PARABOLA_BAND
        kregion = 2
    elseif k² > 2.0
        kregion = 3
    else
        kregion = 1
    end

    if kregion == 2
        return _russell_W_parabola_series(ν, order)
    end

    m = 2.0 - k²
    onebym = 1.0 / m
    tnpp = twoπN + T(π)

    W = 0.0
    if kregion == 0
        # Taylor series near k=0 for precision
        k³ = k² * k
        k⁴ = k² * k²
        k⁵ = k⁴ * k
        k⁶ = k⁴ * k²
        k⁷ = k⁴ * k³
        k⁸ = k⁴ * k⁴
        W =
            0.3535533905932738 * tnpp - k + 0.2651650429449553 * tnpp * k² - 2.0/3.0 * k³ +
            0.1657281518405971 * tnpp * k⁴ - 0.4 * k⁵ + 0.09667475524034829 * tnpp * k⁶ -
            2.0/7.0 * k⁷ + 0.05437954982269591 * tnpp * k⁸
    elseif kregion == 1
        # Normal ellipse
        k²m1 = k² - 1.0
        if k > 0.0
            W = (twoπN + acos(k²m1)) * sqrt(onebym^3) - k * onebym
        else
            kps2 = k + _RUSSELL_SQRT2
            if kps2 < _RUSSELL_K_CLOSE_MSQRT2
                tNp = 2 * T(π) + twoπN
                t1 = kps2^2
                t2 = sqrt(kps2)
                t3 = t2 * t1
                t9 = t1^2
                t10 = t2 * t9
                W =
                    0.12110150049603174603174603174603e-7 / t3 * (
                        -0.38926398009946925989672338336519e8 * t3 -
                        0.16515072e8 * t2 * kps2 * t1 -
                        0.1976320e7 * t10 * (kps2 + 2.4532575164897006338798206469711) -
                        0.18246749067162621557658908595243e7 * t10 +
                        0.25959796716951899525909665607350e6 *
                        (
                            t9 +
                            6.4646464646464646464646464646465 * t1 +
                            35.463203463203463203463203463203
                        ) *
                        tNp *
                        t1 +
                        0.66750357442839860425810740303391e6 *
                        tNp *
                        (
                            t9 +
                            6.0952380952380952380952380952381 * t1 +
                            26.006349206349206349206349206349
                        ) *
                        kps2 - 645120.0 * t2 * kps2 * t9
                    )
            else
                W = (2 * T(π) + twoπN - acos(k²m1)) * sqrt(onebym^3) - k * onebym
            end
        end
    else
        # Hyperbola (kregion == 3)
        k²m1 = k² - 1.0
        W = -log(k²m1 + sqrt(k²m1^2 - 1.0)) * sqrt((-onebym)^3) - k * onebym
    end

    t2 = 3.0 * W
    dW1 = order >= 1 ? (t2 * k - 2.0) * onebym : 0.0
    dW2 = order >= 2 ? (5.0 * dW1 * k + t2) * onebym : 0.0
    dW3 = order >= 3 ? (7.0 * dW2 * k + 8.0 * dW1) * onebym : 0.0

    return (W, dW1, dW2, dW3)
end

"""
    _russell_W_parabola_series!(dW, ν, order)

Evaluate W(k) and derivatives using Taylor series near the parabola (k ≈ √2),
where ν = k - √2. Coefficients from Russell's ivLam (32-digit precision).
"""
@inline function _russell_W_parabola_series(ν::Number, order::Int)
    ν² = ν^2
    ν³ = ν² * ν
    ν⁴ = ν²^2
    ν⁵ = ν⁴ * ν
    ν⁶ = ν⁴ * ν²
    ν⁷ = ν⁴ * ν³
    ν⁸ = ν⁴^2

    W =
        0.47140452079103168 - 0.2 * ν + 0.080812203564176860 * ν² -
        0.031746031746031746 * ν³ + 0.012244273267299524 * ν⁴ - 0.0046620046620046620 * ν⁵ +
        0.0017581520588942907 * ν⁶ - 0.00065816536404771699 * ν⁷ +
        0.00024494378529487022 * ν⁸

    dW1 =
        order >= 1 ?
        (
            -0.2 + 0.16162440712835372 * ν - 0.095238095238095238 * ν² +
            0.048977093069198097 * ν³ - 0.023310023310023310 * ν⁴ +
            0.010548912353365744 * ν⁵ - 0.0046071575483340189 * ν⁶ +
            0.0019595502823589617 * ν⁷
        ) : 0.0

    dW2 =
        order >= 2 ?
        (
            0.16162440712835372 - 0.19047619047619048 * ν + 0.14693127920759429 * ν² -
            0.093240093240093240 * ν³ + 0.052744561766828720 * ν⁴ -
            0.027642945290004114 * ν⁵ + 0.013716851976512732 * ν⁶
        ) : 0.0

    dW3 =
        order >= 3 ?
        (
            -0.19047619047619048 + 0.29386255841518858 * ν - 0.27972027972027972 * ν² +
            0.21097824706731488 * ν³ - 0.13821472645002057 * ν⁴ + 0.082301111859076392 * ν⁵
        ) : 0.0

    return (W, dW1, dW2, dW3)
end

# ============================================================================
# TOF residual and derivatives (including big-k series for huge hyperbolas)
# ============================================================================

"""
    _russell_tof_derivs!(dfunc, k, p, τ, τ², τ³, S, tofbyS, logTofbyS,
                         dW, M, order, huge_tof, huge_k)

Compute the normalized TOF/S residual and its derivatives wrt k. Stores
results in `dfunc[1:order+1]` where `dfunc[1]` is the residual.

Uses a log-form root-solve for very large TOF/S to maintain precision,
and a big-k series expansion for extreme hyperbolic cases.
"""
@inline function _russell_tof_derivs(
    k::Number,
    p::Number,
    τ::Number,
    τ²::Number,
    τ³::Number,
    S::Number,
    tofbyS::Number,
    logTofbyS::Number,
    W::Number,
    dW1::Number,
    dW2::Number,
    dW3::Number,
    M::Int,
    order::Int,
    huge_tof::Bool,
    huge_k::Bool,
)
    if huge_k
        return _russell_tof_bigk(k, p, τ, τ², S, tofbyS, order)
    end

    sqrtp = sqrt(p)

    LeftSide = sqrtp * (p * W + τ)

    F1 = 0.0
    F2 = 0.0
    F3 = 0.0

    if order >= 1
        t7 = p^2
        onebyrootp = 1.0 / sqrtp
        F1 = (-3.0 * p * τ * W + 2.0 * t7 * dW1 - τ²) * onebyrootp * 0.5
    end
    if order >= 2
        onebyp = onebyrootp^2
        onebyp32 = onebyrootp * onebyp
        p3dw0 = 3.0 * p * W
        t7τ = t7 * τ
        F2 = (p3dw0 * τ² + 4.0 * t7 * p * dW2 - 12.0 * t7τ * dW1 - τ³) * onebyp32 * 0.25
    end
    if order >= 3
        F3 =
            (
                p3dw0 * τ³ + 18.0 * t7 * τ² * dW1 - 36.0 * t7τ * p * dW2 +
                8.0 * t7^2 * dW3 - 3.0 * τ²^2
            ) *
            onebyp32 *
            onebyp *
            0.125
    end

    F0 = 0.0
    if huge_tof
        RH = logTofbyS
        F0 = log(LeftSide) - RH
        t1 = 1.0 / LeftSide
        t10 = F1
        t20 = F2
        t30 = F3
        F1 = t10 * t1
        t4 = F1^2
        if order >= 2
            F2 = t20 * t1 - t4
        end
        if order >= 3
            t2 = t1^2
            F3 = t30 * t1 - 3.0 * t20 * t10 * t2 + 2.0 * t4 * F1
        end
    else
        F0 = LeftSide - tofbyS
    end

    return (F0, F1, F2, F3)
end

"""
    _russell_tof_bigk!(dfunc, k, p, τ, τ², S, tofbyS, order)

Series expansion of TOF/S residual and derivatives for very large k
(extreme hyperbolic cases) to maintain numerical precision.
"""
@inline function _russell_tof_bigk(
    k::Number,
    p::Number,
    τ::Number,
    τ²::Number,
    S::Number,
    tofbyS::Number,
    order::Int,
)
    log2 = log(2.0)
    minp = -p
    t1 = k^2
    t3 = minp + 1.0
    t5 = minp * (t1 + 3.0)
    t4 = 1.0 / k
    t8 = log(t4)
    t13 = t1^2
    t14 = t1 * k
    t15 = t14 * τ
    t16 = 2.0 * t15

    t42 = t4^2
    t72 = t42^2
    dQ0 = t72 * t4 * (log2 * t5 - 2.0 * t8 * t5 + 2.0 * t1 + t13 - t16 - 5.0 * t3 + 5.0)

    t102 = -minp
    t103 = sqrt(t102)
    F0 = t103 * dQ0 - tofbyS

    F1 = 0.0
    F2 = 0.0
    F3 = 0.0

    if order >= 1
        t32_a = 6.0 * t15
        t73 = (-t16 + 3.0 * t1 - 12.0 * t3 + 15.0)
        t70 = (-2.0 * t8 + log2)
        t43 = t72 * t42
        dQ1 = t43 * (t70 * t73 - t13 + t32_a - 8.0 * t1 + 26.0 * t3 - 31.0)
        t32_b = dQ0 * τ
        t106 = 1.0 / t103
        F1 = -t106 * t32_b * 0.5 + t103 * dQ1

        if order >= 2
            t74 = (t32_b - 12.0 * t1 + 60.0 * t3 - 90.0)
            dQ2 =
                t43 *
                t4 *
                (t70 * t74 + 2.0 * t14 - 22.0 * t15 + 38.0 * t1 - 154.0 * t3 + 216.0)
            t107 = t106^2
            t113 = t106 * t107
            t115 = τ²
            F2 = -t113 * dQ0 * t115 * 0.25 - t106 * dQ1 * τ + t103 * dQ2

            if order >= 3
                t75 = (-24.0 * t15 + 60.0 * t1 - 360.0 * t3 + 630.0)
                dQ3 =
                    t72^2 * (
                        t70 * t75 - 6.0 * t14 + 100.0 * t15 - 214.0 * t1 + 1044.0 * t3 -
                        1692.0
                    )
                F3 =
                    -0.375 * t113 * t107 * t32_b * t115 - 0.75 * t113 * dQ1 * t115 -
                    1.5 * t106 * dQ2 * τ + t103 * dQ3
            end
        end
    end

    return (F0, F1, F2, F3)
end

# ============================================================================
# Higher-order correction with convergence checks
# ============================================================================

"""
    _russell_correction!(dvars, dfunc, order, is_zero_rev)

Compute the root-solve correction up to the specified order with convergence
series checking and trust-region safeguards for multi-rev cases.

Returns `info` flag (0 = normal, nonzero = safeguard activated).
"""
@inline function _russell_correction(
    F0::Number,
    F1::Number,
    F2::Number,
    F3::Number,
    order::Int,
    is_zero_rev::Bool,
)
    dv1 = 0.0
    dv2 = 0.0
    dv3 = 0.0
    info = 0
    monebydf = 0.0

    if order >= 1
        monebydf = -1.0 / F1
        dv1 = F0 * monebydf
    end

    # Trust region for multi-rev
    if !is_zero_rev
        if dv1 > _RUSSELL_MAX_K_STEP_MREV
            return (3, _RUSSELL_MAX_K_STEP_MREV, 0.0, 0.0)
        elseif dv1 < -_RUSSELL_MAX_K_STEP_MREV
            return (3, -_RUSSELL_MAX_K_STEP_MREV, 0.0, 0.0)
        end
    end

    dvm = 0.0
    dkA = 0.0
    if order >= 2
        dvm = dv1 * monebydf
        dkA = dv1 * dvm
        dv2 = 0.5 * dkA * F2

        absdv1 = abs(dv1)
        absdv2 = abs(dv2)
        if absdv2 * _RUSSELL_SERIES_CONVERGE > absdv1
            if absdv1 > _RUSSELL_SERIES_TAMP_THRESH
                dv1 *= _RUSSELL_SERIES_TAMP
            end
            return (2, dv1, 0.0, 0.0)
        end
    end

    if order >= 3
        dv3 = dkA * dv1 * F3 / 6.0 + dv2 * F2 * dvm

        absdv3 = abs(dv3)
        if absdv3 * _RUSSELL_SERIES_CONVERGE^2 > abs(dv1)
            dv3 = 0.0
            info = 3
        end
    end

    return (info, dv1, dv2, dv3)
end

# ============================================================================
# Initial guess for zero-rev (adapted from Arora & Russell 2013)
# ============================================================================

"""
    _russell_initial_guess_zerorev(τ, S, tof_hat, d)

Compute initial guess for k in the zero-revolution case using the
piecewise rational interpolation from Arora & Russell (2013).
"""
function _russell_initial_guess_zerorev(τ::Float64, S::Float64, tof_hat::Float64, d::Int)
    sq2 = _RUSSELL_SQRT2

    # Parabolic TOF
    tof_p = S * sqrt(1.0 - sq2 * τ) * (τ + sq2) / 3.0

    if tof_hat <= tof_p
        # Hyperbolic regime
        if d == 1
            k_n = sq2
            k_m = 1.0 / τ
            k_i = (k_n + k_m) / 2.0
            Z = 1.0 / sq2
            α = 0.5
            Wi, _, _, _ = _russell_W(k_i, 0, 1)
            F_i = S * sqrt(1.0 - k_i * τ) * (τ + (1.0 - k_i * τ) * Wi)
            x_star = _arora_x(tof_p, 0.0, F_i, tof_hat, Z, α)
            return k_n + (k_m - k_n) * x_star
        else
            tof20 = S * sqrt(1.0 - 20.0 * τ) * (τ + 0.04940968903 * (1.0 - 20.0 * τ))
            tof100 = S * sqrt(1.0 - 100.0 * τ) * (τ + 0.00999209404 * (1.0 - 100.0 * τ))

            if tof_hat > tof20
                # H1 region
                k_n = sq2
                k_m = 20.0
                k_i = (2.0 * k_n + k_m) / 3.0
                Z = 1.0 / 3.0
                α = 1.0
                Wi2, _, _, _ = _russell_W(k_i, 0, 1)
                F_i = S * sqrt(1.0 - k_i * τ) * (τ + (1.0 - k_i * τ) * Wi2)
                x_star = _arora_x(tof_p, tof20, F_i, tof_hat, Z, α)
                return k_n + (k_m - k_n) * x_star
            else
                # H2 region
                return (
                    (
                        tof100 * (tof20 - tof_hat) * 10.0 -
                        tof20 * sqrt(20.0) * (tof100 - tof_hat)
                    ) / (tof_hat * (tof20 - tof100))
                )^2
            end
        end
    else
        # Elliptic regime
        k_set = SVector{6}(-1.41, -1.38, -1.0, -0.5, 0.0, 1.0/sq2)
        W_set = SVector{6}(
            4839.684497246,
            212.087279879,
            5.712388981,
            1.954946607,
            1.110720735,
            0.6686397730,
        )
        _tof(k, W) = S * sqrt(1.0 - k * τ) * (τ + (1.0 - k * τ) * W)
        tof_m141 = _tof(k_set[1], W_set[1])
        tof_m138 = _tof(k_set[2], W_set[2])
        tof_m1 = _tof(k_set[3], W_set[3])
        tof_m1half = _tof(k_set[4], W_set[4])
        tof_0 = _tof(k_set[5], W_set[5])
        tof_1oversq2 = _tof(k_set[6], W_set[6])

        if tof_hat <= tof_0
            x_star = _arora_x(tof_0, tof_p, tof_1oversq2, tof_hat, 0.5, 1.0)
            return (sq2) * x_star
        elseif tof_hat <= tof_m1
            x_star = _arora_x(tof_0, tof_m1, tof_m1half, tof_hat, 0.5, 1.0)
            return -x_star
        elseif tof_hat <= tof_m138
            c1 = 540649.0 / 3125.0
            c2 = 256.0
            c3 = 1.0
            c4 = 1.0
            α = 16.0
            F_n = 1.0 / tof_m1
            F_i = 1.0 / tof_m138
            F_star = 1.0 / tof_hat
            γ1 = F_i * (F_star - F_n)
            γ2 = F_star * (F_n - F_i)
            γ3 = F_n * (F_star - F_i)
            return -c4 *
                   (
                ((γ1 * c1 - c3 * γ3) * c2 + c3 * c1 * γ2) / (γ3 * c1 - c3 * γ1 - γ2 * c2)
            )^(1.0/α)
        else
            c1 = 49267.0 / 27059.0
            c2 = 67286.0 / 17897.0
            c3 = 2813.0 / 287443.0
            c4 = 4439.0 / 3156.0
            α = 243.0
            F_n = 1.0 / tof_m138
            F_i = 1.0 / tof_m141
            F_star = 1.0 / tof_hat
            γ1 = F_i * (F_star - F_n)
            γ2 = F_star * (F_n - F_i)
            γ3 = F_n * (F_star - F_i)
            return -c4 *
                   (
                ((γ1 * c1 - c3 * γ3) * c2 + c3 * c1 * γ2) / (γ3 * c1 - c3 * γ1 - γ2 * c2)
            )^(1.0/α)
        end
    end
end

@inline function _arora_x(F_0, F_1, F_i, F_star, Z, α)
    return (
        (Z * (F_0 - F_star) * (F_1 - F_i)) /
        ((F_i - F_star) * (F_1 - F_0) * Z + (F_0 - F_i) * (F_1 - F_star))
    )^(1.0/α)
end

# ============================================================================
# Initial guess for multi-rev
# ============================================================================

"""
    _russell_initial_guess_multirev(τ, S, tofbyS, M, low_path)

Compute initial guess for k in the multi-revolution case using the
asymptotic k_bottom formula and TOF estimation.

Returns `(k0, tofbyS_bottom)` or `(NaN, NaN)` if no solution exists.
"""
function _russell_initial_guess_multirev(
    τ::T,
    S::V,
    tofbyS::W,
    M::Int,
    low_path::Bool,
) where {T<:Number,V<:Number,W<:Number}

    RT = promote_type(T, V, W)

    abstau = abs(τ)
    τ² = τ^2

    # Compute k_bottom: asymptotic value as N→∞ (works well for all N)
    if abstau < 1e-3
        k_bottom = τ + 0.5 * τ^3
    else
        k_bottom = (1.0 - sqrt(1.0 - 2.0 * τ²)) / τ
    end

    # Improve k_bottom with one Newton-like iteration using getD4W
    twoπN = 2 * RT(π) * M
    Wb, dWb1, dWb2, dWb3 = _russell_W(k_bottom, M, 3)

    p_bot = 1.0 - k_bottom * τ
    tofbyS_bot = sqrt(p_bot) * (τ + p_bot * Wb)

    t1 = τ²
    t3 = k_bottom * t1 - τ
    t6 = k_bottom^2
    t11 = -4.0 * k_bottom * τ + 2.0 * t6 * t1 + 2.0
    y = t11 * dWb1 + 3.0 * t3 * Wb - t1
    dyp = 3.0 * t1 * Wb + t11 * dWb2 + 7.0 * t3 * dWb1

    if abs(dyp) > 1e-30
        dk = -y / dyp
        k_bottom += dk

        Wb_new, _, _, _ = _russell_W(k_bottom, M, 1)
        p_bot = 1.0 - k_bottom * τ
        if p_bot > 0.0
            tofbyS_bot = sqrt(p_bot) * (τ + p_bot * Wb_new)
        end
    end

    # Check if solution exists
    if tofbyS < tofbyS_bot
        return NaN, NaN
    end

    sq2 = _RUSSELL_SQRT2
    kmargin = _RUSSELL_K_MARGIN_EPS

    if low_path
        # k > k_bottom, toward √2
        k_right = sq2 - kmargin
        # Simple linear interpolation as initial guess
        frac = min(0.8, (tofbyS - tofbyS_bot) / (tofbyS + tofbyS_bot))
        k0 = k_bottom + frac * (k_right - k_bottom)
        k0 = clamp(k0, k_bottom + _RUSSELL_K_BUMP_INIT, k_right)
    else
        # k < k_bottom, toward -√2
        k_left = -sq2 + kmargin
        frac = min(0.8, (tofbyS - tofbyS_bot) / (tofbyS + tofbyS_bot))
        k0 = k_bottom - frac * (k_bottom - k_left)
        k0 = clamp(k0, k_left, k_bottom - _RUSSELL_K_BUMP_INIT)
    end

    return k0, tofbyS_bot
end

# ============================================================================
# Velocity reconstruction from converged p value
# ============================================================================

"""
    _russell_velocity_from_p(p, τ, S, r1vec, r2vec, r1, r2, r1pr2, addTog)

Compute departure and arrival velocity vectors from the converged semi-parameter
related variable p using Lagrange f, g, ġ coefficients.
"""
function _russell_velocity_from_p(
    p::Float64,
    τ::Float64,
    S::Float64,
    r1vec::AbstractVector{<:Number},
    r2vec::AbstractVector{<:Number},
    r1::Float64,
    r2::Float64,
    r1pr2::Float64,
    addTog::Float64,
)
    sqrtp = sqrt(p)
    pr12 = p * r1pr2
    f = 1.0 - pr12 / r1
    g = S * τ * sqrtp + addTog
    gdot = 1.0 - pr12 / r2
    onebyg = 1.0 / g

    r1s = SVector{3}(r1vec[1], r1vec[2], r1vec[3])
    r2s = SVector{3}(r2vec[1], r2vec[2], r2vec[3])

    v1 = (r2s - f * r1s) * onebyg
    v2 = (gdot * r2s - r1s) * onebyg
    return v1, v2
end

# ============================================================================
# Main solver function
# ============================================================================

"""
    russell2021(μ, r1, r2, tof; M=0, prograde=true, low_path=true, maxiter=50, rtol=1e-14, order=3)

Solve Lambert's problem using Russell's vercosine formulation with up to
third-order iteration corrections.

The algorithm normalizes to canonical units (μ=1), computes the geometry
parameter τ and characteristic scale S, then iterates on the independent
variable k using Newton/Halley/third-order corrections with safeguards.

# Arguments
- `μ`: Gravitational parameter [L³/T²]
- `r1`: Initial position vector [L]
- `r2`: Final position vector [L]
- `tof`: Time of flight [T]
- `M`: Number of complete revolutions (default: 0)
- `prograde`: Direction of motion (default: true)
- `low_path`: For multi-rev, select low-energy path (default: true)
- `maxiter`: Maximum iterations (default: 50)
- `rtol`: Relative tolerance (default: 1e-14)
- `order`: Correction order 1-3 (default: 3)

# Returns
- `v1`: Departure velocity vector [L/T]
- `v2`: Arrival velocity vector [L/T]
- `numiter`: Number of iterations used
- `retcode`: `:SUCCESS` or `:MAXIMUM_ITERATIONS`

# References
[1] Russell (2019), Celestial Mechanics and Dynamical Astronomy, 131:50
[2] Russell (2021), Journal of Guidance, Control, and Dynamics, 45(2)
"""
function russell2021(
    μ::Number,
    r1::AbstractVector{<:Number},
    r2::AbstractVector{<:Number},
    tof::Number;
    M::Int = 0,
    prograde::Bool = true,
    low_path::Bool = true,
    maxiter::Int = 50,
    rtol::Float64 = 1e-14,
    order::Int = 3,
)
    assert_parameters_are_valid(μ, r1, r2, tof, M)
    order = clamp(order, 1, 3)

    # Compute transfer angle and direction
    _, _, _, dtheta = lambert_geometry(r1, r2, prograde)
    assert_transfer_angle_not_zero(dtheta)
    d = dtheta <= π ? 1 : -1

    # Normalize to canonical units (μ = 1)
    r1_norm = norm(r1)
    Lref = r1_norm
    Tref = sqrt(Lref^3 / μ)
    r1_hat = r1 / Lref
    r2_hat = r2 / Lref
    r2_norm = norm(r2_hat)
    r1hat_norm = 1.0  # norm(r1_hat) by construction
    r2hat_norm = r2_norm
    tof_hat = Float64(tof / Tref)

    # Geometry (from Russell's getGeom)
    r1pr2 = r1hat_norm + r2hat_norm
    r1r2 = r1hat_norm * r2hat_norm
    ctheta = dot(r1_hat, r2_hat) / (r1hat_norm * r2hat_norm)
    ctheta = clamp(ctheta, -1.0, 1.0)
    onePctheta = ctheta + 1.0
    oneMctheta = 1.0 - ctheta

    # Compute τ with precision-saving alternative near θ ≈ π
    if onePctheta < _RUSSELL_ALT_TAU_THRESH
        r1xr2 = cross(
            SVector{3}(r1_hat[1], r1_hat[2], r1_hat[3]),
            SVector{3}(r2_hat[1], r2_hat[2], r2_hat[3]),
        )
        sthetar1r2 = norm(r1xr2)
        abstau = sqrt(1.0 / (oneMctheta * r1r2)) / r1pr2 * sthetar1r2
    else
        abstau = sqrt(r1r2 * onePctheta) / r1pr2
    end
    τ = Float64(d) * abstau

    τ² = τ^2
    τ³ = τ² * τ
    S = r1pr2 * sqrt(r1pr2)
    tofbyS = tof_hat / S

    # Near half-rev detection
    addTog = 0.0
    sqrttiny = sqrt(floatmin(Float64))
    if onePctheta < 3.8e-5
        if abstau < sqrttiny
            addTog = sqrttiny^0.5
        end
    end

    huge_tof = tofbyS > _RUSSELL_HUGE_TOF
    logTofbyS = huge_tof ? log(tofbyS) : 0.0

    W0 = 0.0;
    W1 = 0.0;
    W2 = 0.0;
    W3 = 0.0
    F0 = 0.0;
    F1 = 0.0;
    F2 = 0.0;
    F3 = 0.0
    dv1 = 0.0;
    dv2 = 0.0;
    dv3 = 0.0

    is_zero_rev = (M == 0)

    # ---- Initial guess ----
    if is_zero_rev
        k = _russell_initial_guess_zerorev(τ, S, tof_hat, d)
    else
        k, tofbyS_bot = _russell_initial_guess_multirev(τ, S, tofbyS, M, low_path)
        if isnan(k)
            return zero(SVector{3,Float64}), zero(SVector{3,Float64}), 0, :NO_SOLUTION
        end
    end

    p = 1.0 - k * τ

    # ---- Iteration bounds ----
    sq2 = _RUSSELL_SQRT2
    kmargin = _RUSSELL_K_MARGIN_EPS
    one_m_margin = 1.0 - kmargin

    k_left = -sq2 + kmargin
    if is_zero_rev && τ > 0.0
        k_right = (1.0 / τ) * one_m_margin
    elseif is_zero_rev
        k_right = 1e90
    else
        k_right = sq2 - kmargin
    end

    # Clamp initial guess to valid range
    k_left_start = k_left + _RUSSELL_K_BUMP_INIT
    k_right_start = k_right - _RUSSELL_K_BUMP_INIT
    k = clamp(k, k_left_start, k_right_start)
    p = 1.0 - k * τ

    iterate_on_p = is_zero_rev && p < _RUSSELL_LITTLE_P

    sign_bump =
        is_zero_rev ? 0.0 :
        (low_path ? _RUSSELL_MAX_K_STEP_MREV : -_RUSSELL_MAX_K_STEP_MREV)
    deltaVarTol = _RUSSELL_DELTA_VAR_TOL

    # ---- Main iteration loop ----
    numiter = 0
    converged = false

    for iter = 1:maxiter
        numiter = iter
        k_last = k

        huge_k = is_zero_rev && k > _RUSSELL_BIG_K
        if !huge_k
            W0, W1, W2, W3 = _russell_W(k, M, order)
        end

        F0, F1, F2, F3 = _russell_tof_derivs(
            k,
            p,
            τ,
            τ²,
            τ³,
            S,
            tofbyS,
            logTofbyS,
            W0,
            W1,
            W2,
            W3,
            M,
            order,
            huge_tof,
            huge_k,
        )

        wrong_side = !is_zero_rev && (sign_bump * F1 < 0.0)

        if wrong_side
            k = k_last + sign_bump
            p = 1.0 - k * τ
        else
            if iterate_on_p
                pF1 = order >= 1 ? F1 / (-τ) : F1
                pF2 = order >= 2 ? F2 / τ² : F2
                pF3 = order >= 3 ? F3 / (-τ³) : F3

                info_corr, dv1, dv2, dv3 =
                    _russell_correction(F0, pF1, pF2, pF3, order, is_zero_rev)
                dvar = dv1 + dv2 + dv3

                p += dvar
                comp_val = _RUSSELL_LITTLE_P
                k = (1.0 - p) / τ
            else
                info_corr, dv1, dv2, dv3 =
                    _russell_correction(F0, F1, F2, F3, order, is_zero_rev)
                dvar = dv1 + dv2 + dv3

                k += dvar
                comp_val = k
                p = 1.0 - k * τ
            end
        end

        # Bound violation check
        if k < k_left
            k = (1.0 - _RUSSELL_K_BUMP_WEIGHT) * k_last + _RUSSELL_K_BUMP_WEIGHT * k_left
            p = 1.0 - k * τ
        elseif k > k_right
            k = (1.0 - _RUSSELL_K_BUMP_WEIGHT) * k_last + _RUSSELL_K_BUMP_WEIGHT * k_right
            p = 1.0 - k * τ
        end

        if F0 == 0.0
            converged = true
            break
        end

        # Relative convergence check
        if !wrong_side && maxiter > 1
            hoc = order >= 3 ? dv3 : (order >= 2 ? dv2 : dv1)
            if hoc * deltaVarTol + comp_val == comp_val
                if abs(F0) / max(1.0, abs(huge_tof ? logTofbyS : tofbyS)) < 1.0
                    converged = true
                    break
                end
            end
        end
    end

    retcode = converged ? :SUCCESS : handle_max_iterations(numiter, maxiter)

    if retcode != :SUCCESS
        return zero(SVector{3,Float64}), zero(SVector{3,Float64}), numiter, retcode
    end

    # Reconstruct velocities
    v1_hat, v2_hat = _russell_velocity_from_p(
        p,
        τ,
        S,
        r1_hat,
        r2_hat,
        r1hat_norm,
        r2hat_norm,
        r1pr2,
        addTog,
    )

    # Denormalize
    vscale = Lref / Tref
    v1 = v1_hat * vscale
    v2 = v2_hat * vscale

    return v1, v2, numiter, :SUCCESS
end
