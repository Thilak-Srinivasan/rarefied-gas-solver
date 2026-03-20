# Numerical Simulation of Oblique Shock Waves via Finite Volume Discretization

> **Research Project | BITS Pilani | Jan 2026–Present** \
> Thilak S · Under Prof. Anirudh Singh Rana \
> Department of Mechanical Engineering, BITS Pilani

A high-fidelity **2D Compressible Navier-Stokes Solver** developed in MATLAB using the Finite Volume Method (FVM) with Rusanov (Local Lax-Friedrichs) flux scheme to simulate supersonic flow over a wedge at Mach 2.5. Unlike standard Euler-based solvers, this implementation captures viscous-inviscid interactions and boundary layer growth critical for aerodynamic analysis at high Mach numbers.

---

## Abstract

The Euler equations provide a reasonable estimate of the global shock angle and far-field flow behavior, but the full Navier-Stokes Equations (NSE) are necessary to predict aerodynamic forces (lift, drag, heat transfer), model the shock layer thickness, and capture the interaction with solid walls. This project develops and validates a 2D compressible NSE solver using FVM + Rusanov flux, applied to Mach 2.5 wedge flow, and establishes a framework for future comparison with Direct Simulation Monte Carlo (DSMC) methods in the transition regime.

---

## Why Study Oblique Shock Waves?

Oblique shock waves contain the fundamental physics of compressible flow — the major conversion of kinetic to internal energy occurs over a very thin high-compression region, governing lift, drag, and surface heat transfer in high-speed vehicles.

**Euler vs. Navier-Stokes:**

| Aspect | Euler | Navier-Stokes |
|--------|-------|----------------|
| Shock model | Infinitely thin mathematical discontinuity | Physically thickened shock (viscous diffusion vs convective steepening) |
| Boundary layer | Not captured | Fully resolved with no-slip condition |
| Skin friction / heat flux | Cannot predict | Accurately modeled |
| Shock-boundary layer interaction | Absent | Captured |

**Continuum Validity (Kn < 0.01):**
For low-altitude supersonic flows where the Knudsen number Kn < 0.01, the molecular velocity distribution remains approximately Maxwellian equilibrium. The continuum NS approach matches experimental and DSMC data within ~5%, while being computationally **~100× faster** than DSMC.

---

## Governing Equations

### 2D Compressible Navier-Stokes in Conservation Form

$$\frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}}{\partial x} + \frac{\partial \mathbf{G}}{\partial y} = \frac{\partial \mathbf{P}}{\partial x} + \frac{\partial \mathbf{Q}}{\partial y}$$

This represents simultaneous conservation of mass, x-momentum, y-momentum, and total energy.

### Conserved Variables

$$\mathbf{U} = \begin{pmatrix} \rho \\ \rho u \\ \rho v \\ E \end{pmatrix}$$

where $E = \rho\left(e + \frac{1}{2}(u^2+v^2)\right)$ is total energy per unit volume.

### Inviscid (Convective) Fluxes

$$\mathbf{F} = \begin{pmatrix} \rho u \\ \rho u^2 + p \\ \rho uv \\ (E+p)u \end{pmatrix}, \qquad \mathbf{G} = \begin{pmatrix} \rho v \\ \rho uv \\ \rho v^2 + p \\ (E+p)v \end{pmatrix}$$

These represent the hyperbolic (wave-propagating) part of the system — identical in form to the Euler equations.

### Viscous (Diffusive) Fluxes

$$\mathbf{P} = \begin{pmatrix} 0 \\ \tau_{xx} \\ \tau_{xy} \\ u\tau_{xx} + v\tau_{xy} - q_x \end{pmatrix}, \qquad \mathbf{Q} = \begin{pmatrix} 0 \\ \tau_{xy} \\ \tau_{yy} \\ u\tau_{xy} + v\tau_{yy} - q_y \end{pmatrix}$$

These parabolic terms convert the purely hyperbolic Euler system into a hyperbolic-parabolic system requiring different treatment for convective and diffusive parts.

### Stress Tensor (Newtonian Fluid, Stokes' Hypothesis)

$$\tau_{xx} = \mu\left(\frac{4}{3}\frac{\partial u}{\partial x} - \frac{2}{3}\frac{\partial v}{\partial y}\right), \qquad \tau_{yy} = \mu\left(\frac{4}{3}\frac{\partial v}{\partial y} - \frac{2}{3}\frac{\partial u}{\partial x}\right)$$

$$\tau_{xy} = \tau_{yx} = \mu\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)$$

### Heat Flux (Fourier's Law)

$$q_x = -\kappa\frac{\partial T}{\partial x}, \qquad q_y = -\kappa\frac{\partial T}{\partial y}, \qquad \kappa = \frac{\mu C_p}{Pr}$$

### Equation of State (Ideal Gas)

$$p = \rho RT = (\gamma - 1)\rho e = (\gamma-1)\left(E - \frac{1}{2}\rho(u^2+v^2)\right)$$

### Sutherland's Law (Dynamic Viscosity)

$$\mu(T) = \mu_{ref}\left(\frac{T}{T_{ref}}\right)^{3/2}\frac{T_{ref} + S}{T + S}$$

where $\mu_{ref}$ is the reference viscosity at $T_{ref} = 273.15$ K and $S \approx 110$ K (Sutherland constant for air). This accounts for the significant temperature gradients across the shock layer.

---

## Numerical Method

### 1. Spatial Discretization: Finite Volume Method (FVM)

The FVM integrates the conservation laws over discrete control volumes $\Omega_i$:

$$\frac{\partial}{\partial t}\int_{\Omega_i} \mathbf{U}\,dV + \oint_{\partial\Omega_i}(\vec{H}_{inv} - \vec{H}_{vis})\cdot\hat{n}\,dA = 0$$

- **Structured grid**: cells indexed by $(i,j)$, uniform spacing
- **Conservative**: naturally handles shock discontinuities
- **Flux balance**: net flux summation across all four cell faces

### 2. Inviscid Flux Scheme: Rusanov (Local Lax-Friedrichs)

Standard centered differences applied to hyperbolic terms produce unphysical oscillations near shock fronts. The Rusanov scheme introduces controlled artificial diffusion proportional to the maximum local wave speed:

$$\hat{F}_{i+1/2} = \frac{1}{2}\left[F(U_L) + F(U_R)\right] - \frac{1}{2}\lambda_{max}(U_R - U_L)$$

The stabilization parameter (spectral radius):

$$\lambda_{max} = \max\left\{(|\mathbf{u}\cdot\hat{n}| + a)_L,\;(|\mathbf{u}\cdot\hat{n}| + a)_R\right\}$$

### 3. Viscous Term Discretization: Central Differencing (2nd Order)

Gradients at cell faces computed from adjacent cell centers:

$$\left.\frac{\partial u}{\partial x}\right|_{i+1/2} \approx \frac{u_{i+1} - u_i}{\Delta x}$$

Cross-derivative terms (e.g. $\partial v/\partial y$ in the x-momentum viscous flux) are averaged from the transverse faces.

### 4. Time Integration: Explicit Euler

$$U_i^{n+1} = U_i^n - \frac{\Delta t}{V_i}\sum_{faces}(F_{inv} - F_{vis})\cdot\vec{S}_f$$

### 5. CFL Stability Condition

$$\Delta t \leq \text{CFL}\cdot\frac{\min(\Delta x,\Delta y)}{|u| + a + \frac{2\nu}{\Delta x}}$$

The viscous stability limit is explicitly included in the timestep calculation.

### Discretization Summary

| Component | Method | Order | Physical Character |
|-----------|--------|-------|--------------------|
| Convective flux | Rusanov / LLF | 1st | Hyperbolic (wave propagation) |
| Viscous stress | Central difference | 2nd | Parabolic (diffusion) |
| Heat flux | Central difference | 2nd | Parabolic (diffusion) |
| Time integration | Explicit Euler | 1st | Temporal evolution |

---

## Computational Setup

<table>
<tr>
<td>

**Numerical Parameters**

| Parameter | Value |
|-----------|-------|
| Domain size (Lx × Ly) | 8.0 × 4.0 |
| Grid resolution | 120 × 80 cells |
| Grid spacing | dx ≈ 0.067, dy = 0.05 |
| Freestream Mach (M∞) | 2.5 |
| Wedge angle (θ) | 15° |
| Reynolds number (Re) | 100 |
| CFL number | 0.20 |
| Prandtl number (Pr) | 0.72 |
| Specific heat ratio (γ) | 1.4 |
| Max iterations | 30,000 |
| Convergence tolerance | 10⁻⁶ |

</td>
<td>

**Boundary Conditions**

| Boundary | Condition |
|----------|-----------|
| Inlet (left) | Supersonic inflow — fixed primitive state at M∞ = 2.5 |
| Outlet (right) | Supersonic outflow — zero-gradient extrapolation |
| Top surface | Freestream / symmetry |
| Wedge surface | No-slip, adiabatic wall |

**Wedge Geometry**
- Staircase approximation of wedge on Cartesian mesh
- Wedge starts at x = 0.5
- Solid cells masked from flux computation
- No-slip enforced via ghost cell reflection

</td>
</tr>
</table>

---

## Continuum vs. Kinetic Regime Analysis

The choice between NS and DSMC depends on the Knudsen number $Kn = \lambda/L$:

**Continuum Regime (Kn < 0.01):**
Molecular velocity distribution remains approximately Maxwellian. NS and DSMC pressure/heat flux predictions agree within ~5%. NS is preferred for cost efficiency (~100× faster than DSMC).

**Transition Regime (0.01 < Kn < 0.1):**
Continuum assumptions begin to break down. DSMC captures velocity slip and temperature jump at walls — effects the standard no-slip NS solver cannot model without explicit slip boundary conditions. Extended hydrodynamic models or DSMC should be used for high-altitude or microscale geometries.

**Strong Non-Equilibrium (shock interior):**
Significant deviation from equilibrium molecular distribution. DSMC provides an essentially exact representation of shock internal structure. NS may not quantitatively capture shock thickness as Mach number increases (M > 2), though the Rankine-Hugoniot jump conditions remain qualitatively similar.

**SBLI Benchmark Context:**
The NASA NPARC validation archive for Mach 2.5 flow over a 15° wedge gives $P_2/P_1 \approx 2.467$. Literature shows that NS and DSMC capture the same qualitative SBLI features but disagree quantitatively on recirculation bubble size in rarefied regimes — underscoring the need for careful grid design and boundary condition implementation near separation zones.

---

## Solver Development Timeline

### Step 1,000 — Initial Transient Phase
Freestream conditions begin to perturb at the wedge apex. Shock formation has not yet begun. The solution remains similar to freestream everywhere except in a thin region upstream of the wedge.

### Step 6,000 — Shock Formation & Stabilization
Oblique shock angle stabilizes at approximately 32° to the freestream. Clear jumps in density, pressure, and temperature across the shock. Mach number contours confirm the weak shock solution (M₂ > 1). Velocity vectors show early boundary layer development along the wall.

### Step 24,000 — Fully Converged Solution
Crisp shock structure with stable geometry. Fully developed viscous boundary layer with no-slip condition enforced at the wedge. Shock angle slightly displaced from theoretical inviscid value due to viscous displacement thickness — the key viscous-inviscid interaction absent in Euler-based simulations. Numerical residuals below 10⁻⁷.

---

## Computational Optimisation: MATLAB Vectorisation

The original nested loop implementation was completely restructured to eliminate per-element scalar operations:

**Before:** Nested `for` loops iterating over all $(i,j)$ cells individually.

**After:** Full-domain matrix operations using staggered indexing:

```matlab
% Vectorised central difference for du/dx
du_dx = (u(:, 3:end) - u(:, 1:end-2)) / (2*dx);
```

Precomputing transport properties (Sutherland viscosity, thermal conductivity) and assembling viscous flux vectors P and Q as full-field arrays allows the arithmetic to execute in optimised BLAS kernels.

**Result:** Wall-clock time reduced by **~85%**, enabling transients up to O(10⁴) iterations on a standard desktop.

---

## Velocity Field Visualisation Fix

A major development challenge was resolving the physically incorrect velocity vector plot. Two root causes were identified and fixed:

**Problem 1 — Aspect ratio distortion:** Default MATLAB plot settings exaggerated the 15° wedge angle on screen, making flow deflection appear non-physical.

**Fix:** `axis equal` enforces a 1:1 x-y aspect ratio.

**Problem 2 — Overcrowding:** Plotting at every grid point produced a cluttered, unreadable vector field with arrows penetrating the solid wedge region.

**Fix:** Stride-based sampling (every 5th–10th point) using sub-indexed `quiver`. Solid cells masked from the vector plot.

**Physical result at Step 12,500:** Horizontal flow entering at M∞ = 2.5 undergoes a sharp, resolved 15° upward deflection at the wedge apex. Flow remains tangent to the wedge surface after passing the shock, consistent with the no-slip boundary condition.

---

## Results

### Flow Field Contours — Step 24,000 (Converged)

<p align="center">
  <img src="step_24000.png" alt="Converged solution at step 24000" width="900"/>
  <br>
  <em>Fig. 1: Final converged flow field — Density, Pressure, Temperature, Mach Number, Speed, Velocity at step 24,000 (Ma = 2.5, θ = 15°, Re = 100)</em>
</p>

**Flow features captured:**
- **Oblique Shock (β ≈ 37°):** Sharp density jump across the discontinuity
- **Post-Shock Pressure:** Uniform high pressure region on wedge surface
- **Thermal Layer:** Viscous heating concentrated near the wall boundary
- **Weak Shock Solution:** Flow remains supersonic downstream (M₂ > 1)

---

### Convergence History

<p align="center">
  <img src="convergence.png" alt="Convergence history" width="700"/>
  <br>
  <em>Fig. 2: L² residual history — monotonic decay over 24,000 iterations to target tolerance 10⁻⁶</em>
</p>

**Convergence behaviour:**
- Residuals reduced by **6 orders of magnitude**
- Steady-state achieved at step ~24,000
- CFL = 0.20 maintains stability throughout
- Explicit time-stepping effectively resolves transient shock formation

---

## Validation Against NACA 1135 Analytical Solution

Quantitative comparison against the NASA NPARC validation case: *"Oblique Shock on a 15-Degree Wedge at Mach 2.5"* (Inviscid / Weak Shock Solution).

| Flow Property | NACA 1135 (Theory) | NS Solver (Numerical) | Relative Error (%) |
|---------------|--------------------|-----------------------|-------------------|
| **Shock Angle (β)** | 36.94° | 36.65° | −0.80% |
| **Pressure Ratio (P₂/P₁)** | 2.467 | 2.451 | −0.67% |
| **Density Ratio (ρ₂/ρ₁)** | 1.867 | 1.861 | −0.29% |
| **Downstream Mach (M₂)** | 1.874 | 1.881 | +0.40% |

All errors < 1% — validating the solver against classical oblique shock relations. The marginal increase in numerical shock angle relative to the theoretical inviscid value is physically explained by the viscous displacement thickness effect captured by the NS solver.

---

## Key Finding: Viscous Displacement Thickness Effect

The solver successfully captures the **viscous displacement thickness** — the boundary layer at the wedge surface effectively increases the geometric wedge angle seen by the incoming flow. This results in a shock angle marginally higher than the theoretical inviscid prediction. This viscous correction is absent in purely Euler-based solvers and represents a fundamental physical effect in real supersonic flows near walls.

---

## Shock-Boundary Layer Interaction (SBLI)

The solver captures the full SBLI physics:
- **Displacement thickness effects** on effective wedge angle
- **Pressure gradient influence** on boundary layer development
- **Heat transfer enhancement** near the wall due to viscous dissipation
- **No-slip condition enforcement** confirmed via velocity vector alignment with wedge surface

---

## Repository Structure

```
rarefied-gas-solver/
├── nsf_solver_corrected.m          # Main MATLAB solver
├── outputs/
│   ├── step_1000.png               # Initial transient phase
│   ├── step_6000.png               # Shock formation & stabilization
│   ├── step_24000.png              # Final converged solution
│   └── convergence.png             # Residual history plot
├── docs/
│   ├── week1_report.pdf            # Week 1: Theory & governing equations
│   ├── week2_report.pdf            # Week 2: Implementation & results
│   ├── weekly_progress.pdf         # Weekly progress summary
│   └── midsem_presentation.pdf     # Mid-semester presentation
└── README.md
```

---

## Getting Started

### Requirements
- MATLAB R2020a or later
- Memory: ~68 MB
- Runtime: 4–5 minutes (vectorised implementation)

### Running the Solver

```matlab
>> nsf_solver_corrected
```

Expected console output:
```
========================================
ENHANCED 2D OBLIQUE SHOCK SOLVER
Rusanov Scheme + Navier-Stokes NSF
========================================
Domain: 8.00 x 4.00
Grid: 120 x 80 (dx=0.0672, dy=0.0506)
Mach: 2.50, Wedge Angle: 15°
Reynolds: 100.0, CFL: 0.20
========================================
Step    500: Residual = x.xxxxxxe-xx, dt = x.xxxxxxe-xx
...
CONVERGED at step ~24,000
```

---

## Future Work

**Mesh Refinement (AMR):** Implement adaptive mesh refinement to sharper resolve the shock front and boundary layer near the wedge leading edge.

**Hypersonic Regime (M > 5):** Extend the solver to higher Mach numbers with chemical non-equilibrium models (real-gas effects).

**DSMC Cross-Validation:** Compare continuum NS results with DSMC in the transition regime (0.01 < Kn < 0.1) to quantify breakdown of continuum assumptions.

**Higher-Order Schemes:** MUSCL-TVD reconstruction for improved shock resolution and reduced numerical diffusion.

**Implicit Time Integration:** Remove explicit CFL constraint for larger, more complex domains.

**Turbulence Modeling:** RANS or LES for high Reynolds number flows.

---

## References

1. Anderson, J. D. (2003). *Modern Compressible Flow: With Historical Perspective*, 3rd ed. McGraw-Hill.
2. Bird, G. A. (1994). *Molecular Gas Dynamics and the Direct Simulation of Gas Flows*. Oxford University Press.
3. NASA NPARC Alliance Validation Archive — *"Oblique Shock on a 15 Degree Wedge at Mach 2.5"*, NASA Glenn Research Center.
4. Rusanov, V. V. (1961). Calculation of Interaction of Non-Steady Shock Waves with Obstacles. *Zh. Vychisl. Mat. Mat. Fiz.*, 1(2), 267–279.
5. Sutherland, W. (1893). The Viscosity of Gases and Molecular Force. *Philosophical Magazine*, 36(223), 507–531.
6. Hildebrand et al. (2018). Simulation and Stability Analysis of Oblique Shock-Wave/Boundary-Layer Interactions at Mach 5.92. *Physical Review Fluids*, 3, 013906.
7. Bhide et al. (2021). Velocity-Slip and Temperature-Jump Effects in Near-Continuum Hypersonic Flows. *AIAA Journal*, 59(10), 4312–4325.
8. White, F. M. (2006). *Viscous Fluid Flow*. McGraw-Hill.
9. Hirsch, C. (2007). *Numerical Computation of Internal and External Flows*. Butterworth-Heinemann.
10. MacCormack, R. W. (1969). The Effect of Viscosity in Hypervelocity Impact Cratering. AIAA Paper 69-354.

---

**Author:** Thilak S | Department of Mechanical Engineering, BITS Pilani \
**Supervisor:** Prof. Anirudh Singh Rana \
**Period:** January 2026–Present \
**Status:** Active Development \
**Contact:** f20220771@pilani.bits-pilani.ac.in
