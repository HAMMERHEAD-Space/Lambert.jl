export lambert_jacobian

"""
    lambert_jacobian(prob::LambertProblem, sol::LambertSolution)

Compute the analytical Jacobian of the Lambert solution velocities (v₁, v₂) with respect
to the inputs (μ, r₁, r₂, tof). Returns a `NamedTuple` of sensitivity matrices:

- `dv1_dr1::SMatrix{3,3}`, `dv1_dr2::SMatrix{3,3}`, `dv1_dtof::SVector{3}`, `dv1_dmu::SVector{3}`
- `dv2_dr1::SMatrix{3,3}`, `dv2_dr2::SMatrix{3,3}`, `dv2_dtof::SVector{3}`, `dv2_dmu::SVector{3}`

Uses the two-body state transition matrix (STM) for position/tof sensitivities and
analytical Lagrange coefficient derivatives for μ sensitivity, with zero AD dependency.

# References
- Battin (1999), *An Introduction to the Methods of Astrodynamics*, §9.4
- Arora & Russell (2014), *Partial Derivatives of the Lambert Problem*
"""
function lambert_jacobian(prob::LambertProblem, sol::LambertSolution)
    @unpack μ, r1, r2, tof = prob
    v1 = sol.v1
    v2 = sol.v2

    r1_vec = SVector{3}(r1[1], r1[2], r1[3])
    r2_vec = SVector{3}(r2[1], r2[2], r2[3])
    v1_vec = SVector{3}(v1[1], v1[2], v1[3])
    v2_vec = SVector{3}(v2[1], v2[2], v2[3])

    Phi_rr, Phi_rv, Phi_vr, Phi_vv = _twobody_stm(r1_vec, v1_vec, tof, μ)

    Phi_rv_inv = inv(Phi_rv)
    Phi_rv_inv_Phi_rr = Phi_rv_inv * Phi_rr

    # Position sensitivities from STM boundary-value formulation
    dv1_dr1 = -Phi_rv_inv_Phi_rr
    dv1_dr2 = Phi_rv_inv
    dv2_dr1 = Phi_vr - Phi_vv * Phi_rv_inv_Phi_rr
    dv2_dr2 = Phi_vv * Phi_rv_inv

    # TOF sensitivities
    dv1_dtof = -Phi_rv_inv * v2_vec
    a2 = -(μ / norm(r2_vec)^3) * r2_vec
    dv2_dtof = a2 - Phi_vv * (Phi_rv_inv * v2_vec)

    # μ sensitivity via STM + variational equation for ∂r₂/∂μ at fixed (r₁,v₁)
    dr2_dmu_orbit, dv2_dmu_orbit =
        _mu_orbit_sensitivity(r1_vec, v1_vec, tof, μ, Phi_rr, Phi_rv, Phi_vr, Phi_vv)
    dv1_dmu = -Phi_rv_inv * dr2_dmu_orbit
    dv2_dmu = dv2_dmu_orbit + Phi_vv * dv1_dmu

    return (
        dv1_dr1 = dv1_dr1,
        dv1_dr2 = dv1_dr2,
        dv1_dtof = dv1_dtof,
        dv1_dmu = dv1_dmu,
        dv2_dr1 = dv2_dr1,
        dv2_dr2 = dv2_dr2,
        dv2_dtof = dv2_dtof,
        dv2_dmu = dv2_dmu,
    )
end

# ============================================================================
# Two-body state transition matrix
# ============================================================================

"""
Compute the 4 blocks (each 3×3 SMatrix) of the two-body STM using
the universal variable formulation with implicit differentiation of
the Kepler equation.
"""
function _twobody_stm(r1_vec, v1_vec, tof, μ)
    r1_mag = norm(r1_vec)
    sqrtmu = sqrt(μ)
    r1dv1 = dot(r1_vec, v1_vec)
    α = 2.0 / r1_mag - dot(v1_vec, v1_vec) / μ

    χ = _solve_kepler_universal(r1_mag, r1dv1, α, sqrtmu, tof)
    ψ = α * χ^2
    C2val = c2(ψ)
    C3val = c3(ψ)

    r2_mag =
        χ^2 * C2val + (r1dv1 / sqrtmu) * χ * (1.0 - ψ * C3val) + r1_mag * (1.0 - ψ * C2val)

    f = 1.0 - (χ^2 / r1_mag) * C2val
    g = tof - (χ^3 / sqrtmu) * C3val
    fdot = (sqrtmu / (r1_mag * r2_mag)) * χ * (ψ * C3val - 1.0)
    gdot = 1.0 - (χ^2 / r2_mag) * C2val

    r2_vec = f * r1_vec + g * v1_vec
    v2_vec = fdot * r1_vec + gdot * v1_vec

    dC2_dψ =
        abs(ψ) > 1e-6 ? (1.0 - ψ * C3val - 2.0 * C2val) / (2.0 * ψ) :
        -1.0 / 24.0 + ψ / 360.0
    dC3_dψ = abs(ψ) > 1e-6 ? (C2val - 3.0 * C3val) / (2.0 * ψ) : -1.0 / 120.0 + ψ / 2520.0

    # --- Implicit differentiation of Kepler equation for dχ/d(state) ---
    r1hat = r1_vec / r1_mag
    χ2 = χ^2

    # ∂K/∂(R, σ, α) where R=|r₁|, σ=r₁·v₁
    dK_dR = χ * (1.0 - ψ * C3val)
    dK_dσ = χ2 * C2val / sqrtmu
    dK_dα =
        (r1dv1 / sqrtmu) * χ2 * dC2_dψ * χ2 +
        r1_mag * χ * (-(C3val + ψ * dC3_dψ) * χ2) +
        χ^3 * dC3_dψ * χ2

    # Chain to state vectors (∂R/∂r₁=r̂₁, ∂σ/∂r₁=v₁, ∂α/∂r₁=-(2/R²)r̂₁)
    dK_dr1 = (dK_dR - 2.0 * dK_dα / r1_mag^2) * r1hat + dK_dσ * v1_vec
    dK_dv1 = dK_dσ * r1_vec + dK_dα * (-2.0 * v1_vec / μ)

    inv_r2 = 1.0 / r2_mag
    dchi_dr1 = -inv_r2 * dK_dr1
    dchi_dv1 = -inv_r2 * dK_dv1

    # --- Lagrange coefficient partials via χ chain rule ---
    # Stumpff identities: d(χ²C₂)/dχ = χ(1-ψC₃), d(χ³C₃)/dχ = χ²C₂
    df_dchi = -χ * (1.0 - ψ * C3val) / r1_mag
    dg_dchi = -χ2 * C2val / sqrtmu

    # ∂f/∂α, ∂g/∂α at fixed χ (through ψ = αχ²)
    df_dα = -(χ2^2 * dC2_dψ) / r1_mag
    dg_dα = -(χ^3 / sqrtmu) * dC3_dψ * χ2

    df_dr1 =
        df_dchi * dchi_dr1 + (df_dα * (-2.0 / r1_mag^2) + χ2 * C2val / r1_mag^2) * r1hat
    df_dv1 = df_dchi * dchi_dv1 + df_dα * (-2.0 * v1_vec / μ)
    dg_dr1 = dg_dchi * dchi_dr1 + dg_dα * (-2.0 / r1_mag^2) * r1hat
    dg_dv1 = dg_dchi * dchi_dv1 + dg_dα * (-2.0 * v1_vec / μ)

    # For fdot and gdot partials, also need dr₂_mag/d(state)
    dr2_dchi =
        (r1dv1 / sqrtmu) * (1.0 - ψ * C2val) + χ * (1.0 - ψ * C3val) * (1.0 - α * r1_mag)

    C3_plus_ψdC3 = (C2val - C3val) / 2.0
    C2_plus_ψdC2 = (1.0 - ψ * C3val) / 2.0
    dr2_dψ_explicit =
        -(r1dv1 / sqrtmu) * χ * C3_plus_ψdC3 - r1_mag * C2_plus_ψdC2 + χ2 * dC2_dψ
    dr2_dα_explicit = dr2_dψ_explicit * χ2

    dr2_dr1_explicit =
        (v1_vec / sqrtmu) * χ * (1.0 - ψ * C3val) +
        r1hat * (1.0 - ψ * C2val) +
        dr2_dα_explicit * (-(2.0 / r1_mag^2) * r1hat)
    dr2_dr1 = dr2_dchi * dchi_dr1 + dr2_dr1_explicit

    dr2_dv1_explicit =
        (r1_vec / sqrtmu) * χ * (1.0 - ψ * C3val) + dr2_dα_explicit * (-2.0 * v1_vec / μ)
    dr2_dv1 = dr2_dchi * dchi_dv1 + dr2_dv1_explicit

    # fdot = √μ/(R r₂) χ(ψC₃-1)
    P = χ * (ψ * C3val - 1.0)
    dP_dchi = -(1.0 - ψ * C2val)
    dP_dα = χ * C3_plus_ψdC3 * χ2
    coeff_fdot = sqrtmu / (r1_mag * r2_mag)

    dfdot_dr1 =
        coeff_fdot * (dP_dchi * dchi_dr1 + dP_dα * (-(2.0 / r1_mag^2) * r1hat)) +
        sqrtmu *
        P *
        (-(1.0 / (r1_mag^2 * r2_mag)) * r1hat - (1.0 / (r1_mag * r2_mag^2)) * dr2_dr1)
    dfdot_dv1 =
        coeff_fdot * (dP_dchi * dchi_dv1 + dP_dα * (-2.0 * v1_vec / μ)) +
        sqrtmu * P * (-(1.0 / (r1_mag * r2_mag^2)) * dr2_dv1)

    # gdot = 1 - χ²C₂/r₂
    B_val = χ2 * C2val
    dB_dchi = χ * (1.0 - ψ * C3val)
    dB_dα = χ2^2 * dC2_dψ

    dgdot_dr1 =
        -(1.0 / r2_mag) * (dB_dchi * dchi_dr1 + dB_dα * (-(2.0 / r1_mag^2) * r1hat)) +
        (B_val / r2_mag^2) * dr2_dr1
    dgdot_dv1 =
        -(1.0 / r2_mag) * (dB_dchi * dchi_dv1 + dB_dα * (-2.0 * v1_vec / μ)) +
        (B_val / r2_mag^2) * dr2_dv1

    # --- Assemble STM blocks ---
    # r₂ = f r₁ + g v₁ → ∂r₂ᵢ/∂r₁ⱼ = f δᵢⱼ + r₁ᵢ(∂f/∂r₁ⱼ) + v₁ᵢ(∂g/∂r₁ⱼ)
    I3 = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

    Phi_rr = f * I3 + r1_vec * df_dr1' + v1_vec * dg_dr1'
    Phi_rv = g * I3 + r1_vec * df_dv1' + v1_vec * dg_dv1'
    Phi_vr = fdot * I3 + r1_vec * dfdot_dr1' + v1_vec * dgdot_dr1'
    Phi_vv = gdot * I3 + r1_vec * dfdot_dv1' + v1_vec * dgdot_dv1'

    return Phi_rr, Phi_rv, Phi_vr, Phi_vv
end

"""
Solve the universal Kepler equation for χ via Newton-Raphson.
"""
function _solve_kepler_universal(
    r1_mag::Number,
    r1dv1::Number,
    α::Number,
    sqrtmu::Number,
    tof::Number;
    maxiter::Int = 50,
    tol::Float64 = 1e-14,
)
    χ = sqrtmu * tof / r1_mag
    for _ = 1:maxiter
        ψ = α * χ^2
        C2val = c2(ψ)
        C3val = c3(ψ)
        r =
            χ^2 * C2val +
            (r1dv1 / sqrtmu) * χ * (1.0 - ψ * C3val) +
            r1_mag * (1.0 - ψ * C2val)
        F =
            (r1dv1 / sqrtmu) * χ^2 * C2val + r1_mag * χ * (1.0 - ψ * C3val) + χ^3 * C3val -
            sqrtmu * tof
        δχ = -F / r
        χ += δχ
        abs(δχ) < tol * max(1.0, abs(χ)) && return χ
    end
    return χ
end

# ============================================================================
# μ sensitivity via orbit-level variational equation
# ============================================================================

"""
Compute ∂r₂/∂μ and ∂v₂/∂μ at fixed (r₁, v₁, tof) — how the propagated state
changes when μ is varied but the initial conditions are held constant.

Uses the Lagrange coefficient derivatives through the universal variable.
r₂ = f(μ) r₁ + g(μ) v₁, v₂ = fdot(μ) r₁ + gdot(μ) v₁.
"""
function _mu_orbit_sensitivity(r1_vec, v1_vec, tof, μ, Phi_rr, Phi_rv, Phi_vr, Phi_vv)
    r1_mag = norm(r1_vec)
    sqrtmu = sqrt(μ)
    r1dv1 = dot(r1_vec, v1_vec)
    v1sq = dot(v1_vec, v1_vec)
    α = 2.0 / r1_mag - v1sq / μ

    χ = _solve_kepler_universal(r1_mag, r1dv1, α, sqrtmu, tof)
    ψ = α * χ^2
    C2val = c2(ψ)
    C3val = c3(ψ)
    χ2 = χ^2
    χ3 = χ2 * χ

    r2_mag =
        χ2 * C2val + (r1dv1 / sqrtmu) * χ * (1.0 - ψ * C3val) + r1_mag * (1.0 - ψ * C2val)

    f = 1.0 - (χ2 / r1_mag) * C2val
    g = tof - (χ3 / sqrtmu) * C3val
    fdot = (sqrtmu / (r1_mag * r2_mag)) * χ * (ψ * C3val - 1.0)
    gdot = 1.0 - (χ2 / r2_mag) * C2val

    dC2_dψ =
        abs(ψ) > 1e-6 ? (1.0 - ψ * C3val - 2.0 * C2val) / (2.0 * ψ) :
        -1.0 / 24.0 + ψ / 360.0
    dC3_dψ = abs(ψ) > 1e-6 ? (C2val - 3.0 * C3val) / (2.0 * ψ) : -1.0 / 120.0 + ψ / 2520.0

    # dχ/dμ from implicit differentiation of the Kepler equation at fixed (r₁, v₁, tof)
    dK_dα =
        (r1dv1 / sqrtmu) * χ2 * dC2_dψ * χ2 +
        r1_mag * χ * (-(C3val + ψ * dC3_dψ) * χ2) +
        χ3 * dC3_dψ * χ2

    dα_dμ = v1sq / μ^2

    # Explicit μ terms in K = (σ/√μ)χ²C₂ + Rχ(1-ψC₃) + χ³C₃ - √μ tof
    dK_dμ_explicit = -r1dv1 * χ2 * C2val / (2.0 * μ * sqrtmu) - tof / (2.0 * sqrtmu)
    dK_dμ = dK_dμ_explicit + dK_dα * dα_dμ

    dchi_dμ = -dK_dμ / r2_mag

    # ∂f/∂μ at fixed (r₁,v₁,tof): f = 1 - χ²C₂/R
    df_dchi = -χ * (1.0 - ψ * C3val) / r1_mag
    df_dα = -(χ2^2 * dC2_dψ) / r1_mag
    df_dμ = df_dchi * dchi_dμ + df_dα * dα_dμ

    # ∂g/∂μ: g = tof - χ³C₃/√μ
    dg_dchi = -χ2 * C2val / sqrtmu
    dg_dα = -(χ3 / sqrtmu) * dC3_dψ * χ2
    dg_dμ_explicit = χ3 * C3val / (2.0 * μ * sqrtmu)
    dg_dμ = dg_dchi * dchi_dμ + dg_dα * dα_dμ + dg_dμ_explicit

    # ∂r₂/∂μ (vector) = (∂f/∂μ) r₁ + (∂g/∂μ) v₁
    dr2_dmu_orbit = df_dμ * r1_vec + dg_dμ * v1_vec

    # For dv₂/dμ we also need ∂fdot/∂μ and ∂gdot/∂μ, but these are complex.
    # Instead, use: dv₂/dμ = Φ_vr·0 + Φ_vv·(dv₁/dμ) + [∂v₂/∂μ]_orbit
    # The caller computes dv₁/dμ = -Φ_rv⁻¹ dr2_dmu_orbit.
    # dv₂/dμ = [∂v₂/∂μ]_orbit + Φ_vv dv₁/dμ
    # where [∂v₂/∂μ]_orbit = (∂fdot/∂μ) r₁ + (∂gdot/∂μ) v₁

    # ∂fdot/∂μ and ∂gdot/∂μ require dr₂_mag/dμ
    P = χ * (ψ * C3val - 1.0)
    dP_dchi = -(1.0 - ψ * C2val)
    C3_plus_ψdC3 = (C2val - C3val) / 2.0
    dP_dα = χ * C3_plus_ψdC3 * χ2

    dr2_dchi =
        (r1dv1 / sqrtmu) * (1.0 - ψ * C2val) + χ * (1.0 - ψ * C3val) * (1.0 - α * r1_mag)
    C2_plus_ψdC2 = (1.0 - ψ * C3val) / 2.0
    dr2_dψ_expl = -(r1dv1 / sqrtmu) * χ * C3_plus_ψdC3 - r1_mag * C2_plus_ψdC2 + χ2 * dC2_dψ
    dr2_dα = dr2_dψ_expl * χ2
    dr2_dμ_explicit_scalar = -(r1dv1 / (2.0 * μ * sqrtmu)) * χ * (1.0 - ψ * C3val)
    dr2mag_dμ = dr2_dchi * dchi_dμ + dr2_dα * dα_dμ + dr2_dμ_explicit_scalar

    coeff_fdot = sqrtmu / (r1_mag * r2_mag)
    dfdot_dμ =
        coeff_fdot * (dP_dchi * dchi_dμ + dP_dα * dα_dμ) +
        P / (2.0 * sqrtmu * r1_mag * r2_mag) - sqrtmu * P / (r1_mag * r2_mag^2) * dr2mag_dμ

    B_val = χ2 * C2val
    dB_dchi = χ * (1.0 - ψ * C3val)
    dB_dα = χ2^2 * dC2_dψ
    dgdot_dμ =
        -(1.0 / r2_mag) * (dB_dchi * dchi_dμ + dB_dα * dα_dμ) +
        (B_val / r2_mag^2) * dr2mag_dμ

    dv2_dmu_orbit = dfdot_dμ * r1_vec + dgdot_dμ * v1_vec

    return dr2_dmu_orbit, dv2_dmu_orbit
end
