# Development Protocol
## Magnetic field and gradient of an axially magnetized cylinder

This document records the engineering protocol followed to develop the
two MATLAB programs in this directory:

- `cylinder_magnet.m`           — brute-force surface integral via `integral2`
- `cylinder_magnet_analytic.m`  — analytical phi-integration + 1-D radial quadrature

---

## 1. Problem statement

Compute the magnetic field **B**(x,y,z) and its full gradient tensor
∂B_α/∂x_β outside a uniformly axially magnetized cylinder
(radius R, length L, magnetization M = M_0 ẑ with μ_0 M_0 = 1 T,
axis aligned with the Z-axis).

Required outputs:
- 1-D on-axis profile of B_z for validation;
- 2-D maps of B_x, B_y, B_z and the gradient tensor in an XY plane at
  a chosen height z = z_0 above the magnet;
- visualization of in-plane field direction and gradient direction.

---

## 2. Physics formulation

### 2.1 Scalar magnetic potential

For uniform magnetization the body is replaced by surface magnetic
charges with density σ_m = M·n̂.  Only the two flat end faces carry
charge for axial magnetization:

| face   | z'      | σ_m   |
|--------|---------|-------|
| top    | +L/2    | +M_0  |
| bottom | -L/2    | -M_0  |

The scalar potential is

    φ_m(r) = (1/4π) ∬_S σ_m / |r - r'| dS'.

### 2.2 Field and gradient

Outside the magnet **H** = -∇φ_m and **B** = μ_0 **H**.
Both **B** and ∂B_α/∂x_β are obtained by differentiating *under the
integral sign* analytically:

    H_α      = (σ/4π) ∬ (α-α')/d^3 dS'
    ∂H_α/∂β  = (σ/4π) ∬ [-δ_αβ/d^3 + 3(α-α')(β-β')/d^5] dS'.

### 2.3 Validation references

- **Near / intermediate field**: closed form on the cylinder axis,
  derived directly from the surface integral at ρ=0.
- **Far field**: point-dipole formula with moment m = M_0·πR²·L.

---

## 3. Implementation strategy

Two implementations were produced, each a logical step of the same plan.

### 3.1 Brute-force (`cylinder_magnet.m`)

- Each of the three field components and each of the five independent
  gradient-tensor components is computed by `integral2` per
  observation point on the (x,y) grid.
- Total: 8 × N² adaptive 2-D integrations.
- Pros: minimal physics derivation, robust.
- Cons: ~minutes for a 50×50 grid; not useful for design iteration.

### 3.2 Analytical reduction (`cylinder_magnet_analytic.m`)

Step-wise improvement:

1. **Analytical φ' integration.**  The 2π integral over the disc
   parameter φ' is evaluated in closed form; the result is expressed
   in terms of the complete elliptic integrals K(m) and E(m) with
   parameter

       m(r',ρ,ζ) = 4 r' ρ / [(r'+ρ)² + ζ²].

2. **Closed-form integrands.**  The remaining integrand depends only
   on r', ρ, ζ:

       H_z (ρ,ζ) = (σ/π) ζ ∫₀^R r' E(m) /
                     [ ((r'-ρ)² + ζ²) · √((r'+ρ)² + ζ²) ] dr'

       H_ρ(ρ,ζ) = σ/(2π ρ) ∫₀^R r'/[((r'-ρ)²+ζ²)·√((r'+ρ)²+ζ²)] ·
                       [ (ρ²-r'²-ζ²) E(m) + ((r'-ρ)²+ζ²) K(m) ] dr'.

3. **Vectorized 1-D quadrature.**  The radial integral on [0, R] is
   evaluated by a single 96-point Gauss-Legendre rule (Golub-Welsch),
   broadcast across the entire observation grid as one matrix
   product.

4. **Gradient tensor.**  Six perturbation evaluations of the analytic
   field at observation points x±h, y±h, z±h give the gradient tensor
   by central differences.  The field is quadrature-noiseless, so a
   step h = 10⁻⁶ m gives ~10⁻⁹ accuracy.

Speed-up vs brute force: ≈ 10³.

---

## 4. Numerical-precision protocol

| issue                              | mitigation                                                                 |
|------------------------------------|-----------------------------------------------------------------------------|
| K(m → 1) logarithmic singularity   | clamp m ≤ 1 - 10⁻¹⁴ (rim only; observation points stay outside the magnet) |
| ρ = 0 (1/ρ in H_ρ formula)         | guard `max(rho, 1e-30)`, then force H_ρ(ρ < 10⁻¹² R) = 0                     |
| 0/0 in Cartesian projection on axis| use mask `nz = rho > eps_r`, set B_x = B_y = 0 on axis                       |
| integration singularity at the rim | observation plane lifted by z0add (default 1 mm) above the magnet            |

---

## 5. Validation protocol

Performed each time the analytic core is modified.

1. **On-axis comparison** at 200 log-spaced points z ∈ [b + 10⁻⁴, 20 L]:
   - numerical (1-D Gauss-Legendre) vs exact closed form  →  must agree
     to relative 10⁻⁹ over the full range;
   - numerical vs dipole far-field  →  agree only for z ≫ L; the
     log-log error plot must asymptote to a 1/r⁴-like decay of the
     residual.
2. **Divergence check** on the 2-D grid: trace of the gradient tensor
   ( ∂B_x/∂x + ∂B_y/∂y + ∂B_z/∂z ) plotted alongside the diagonal
   components.  Outside source-free conditions ⇒ map should be
   numerical noise (≪ component magnitudes).
3. **Cross-check between programs**: at a small set of grid points,
   the brute-force script and the analytic script must agree to 10⁻⁸
   relative tolerance.

---

## 6. Visualization protocol

| figure | content                                                                         |
|--------|---------------------------------------------------------------------------------|
|   1    | on-axis B_z and its relative error vs both references (log-log, 2 panels)       |
|   2    | B_x, B_y as `pcolor`, B_z as `surf`, common color limits, axes in cm            |
|   3    | B_z color map with quiver of normalized in-plane (B_x,B_y) and magnet outline   |
|   4    | five gradient components plus div B (sanity check) on a 2×3 tiled layout        |
|   5    | three components of ∇\|B\| obtained from the gradient tensor and B itself       |

Plots are docked windows, white background, jet colormap for the
field components, default for gradients.

---

## 7. Project parameters (current run)

```
R     = 35 mm        cylinder radius
L     = 60 mm        cylinder length
B_r   = 1.0 T        remanent flux density
z0add = 1 mm         distance above the top face for the XY map
N     = 80           XY grid resolution (analytic version)
Nq    = 96           Gauss-Legendre nodes on [0, R]
hd    = 10⁻⁶ m       central-difference step for the gradient tensor
```

---

## 8. Reproducibility

Both scripts are self-contained, no toolboxes beyond base MATLAB are
required (`ellipke` ships with the core product).  Running each
script from a clean session reproduces all figures.  The brute-force
script is kept in the repository as the slow but transparent
reference implementation; the analytic script is the one used for
production runs.
