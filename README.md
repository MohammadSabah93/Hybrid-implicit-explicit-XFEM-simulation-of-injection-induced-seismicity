# Hybrid IMEX XFEM for Injection-Induced Seismicity

[![MATLAB](https://img.shields.io/badge/MATLAB-research%20code-e16737)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Research software](https://img.shields.io/badge/status-research%20software-4c1)](#scope-and-limitations)
[![ORCID](https://img.shields.io/badge/ORCID-0000--0002--7260--6138-A6CE39)](https://orcid.org/0000-0002-7260-6138)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Mohammad%20Sabah-0A66C2)](https://www.linkedin.com/in/mohammad-sabah)

A MATLAB research implementation developed by **Mohammad Sabah** for a hybrid implicit–explicit (IMEX) strategy in two-dimensional coupled hydromechanical XFEM simulation of **injection-induced fault slip and dynamic rupture**.

The formulation targets the strong timescale separation between slow reservoir pressurization and rapid fault acceleration. The mechanics are integrated implicitly during the quasi-static phase and switched to an explicit central-difference-type update once the maximum fault slip rate exceeds a prescribed threshold. During the explicit mechanical phase, pore pressure is advanced separately with a backward-Euler solve.

The project sits at the intersection of computational geomechanics, induced seismicity, reservoir geomechanics, XFEM, rate-and-state friction, dynamic rupture, and coupled multiphysics simulation.

## Scientific scope

The code combines:

- coupled poroelastic deformation and fluid flow;
- XFEM representation of an embedded fault;
- matrix–fracture hydraulic coupling;
- nonlinear fault contact;
- rate-and-state friction;
- inertia and dynamic boundary damping;
- adaptive time stepping in the implicit regime;
- slip-rate-based implicit/explicit switching; and
- post-processing of displacement, pressure, stress, slip and seismicity metrics.

The current implementation is intended for numerical-method development and research reproducibility rather than operational seismic-hazard forecasting.

## Time-integration strategy

The default driver uses two regimes:

1. **Implicit regime** — Newmark-type mechanics with Newton–Raphson iterations for the coupled system.
2. **Dynamic regime** — explicit mechanical update triggered when the peak fault slip rate exceeds `Vth`; pressure is then advanced with a backward-Euler step.

The switching threshold and dynamic time step are defined in `X_FEM_PoroElastic_Im_Ex.m` and should be treated as numerical parameters requiring sensitivity analysis.

## Quick start

### Requirements

- MATLAB;
- all repository `.m` files available on the MATLAB path.

A minimum MATLAB release and toolbox compatibility matrix have not yet been formally established.

### Run the default case

```bash
git clone https://github.com/MohammadSabah93/Hybrid-implicit-explicit-XFEM-simulation-of-injection-induced-seismicity.git
cd Hybrid-implicit-explicit-XFEM-simulation-of-injection-induced-seismicity
```

Open the repository in MATLAB and run:

```matlab
X_FEM_PoroElastic_Im_Ex
```

Before running, review:

- geometry, mesh, initial conditions, time stepping and switching parameters in `X_FEM_PoroElastic_Im_Ex.m`;
- material, hydraulic, dynamic and frictional parameters in `defineModelParameters.m`;
- mechanical and hydraulic boundary conditions in `defineBoundaryConditions.m`.

The default example uses a 100 m × 100 m domain with an embedded diagonal fault. Runtime depends strongly on mesh resolution, dynamic time step and the duration of the explicit phase.

## Code map

| Component | Main files |
|---|---|
| Driver and IMEX integration | `X_FEM_PoroElastic_Im_Ex.m` |
| Parameters and boundary conditions | `defineModelParameters.m`, `defineBoundaryConditions.m` |
| Mesh and XFEM enrichment | `MeshGeneration_2D.m`, `levelSet.m`, `enrElem.m` |
| Matrix assembly | `StiffnessMatrix.m`, `MassMatrix.m`, `StorageMatrix.m`, `ConductanceMatrix.m`, `CouplingMatrix.m` |
| Fault contact and friction | `Interface.m`, `Lagrange.m`, `StabLagrange.m`, `updateInterface.m`, `updateInterface_expilict.m` |
| Matrix–fracture flow coupling | `interface_flow.m`, `CouplingInterface.m` |
| Dynamic boundaries and seismicity | `dampingBoundary.m`, `SeismicityParameters.m` |
| Visualization | `postProcess.m` |

## Reproducibility notes

For quantitative studies, users should document and test at least:

- mesh resolution;
- implicit and explicit time-step limits;
- slip-rate switching threshold `Vth`;
- nonlinear convergence tolerance;
- rate-and-state friction parameters;
- damping parameters; and
- boundary-condition sensitivity.

A robust scientific use of the code should include mesh- and time-step-convergence studies and comparison against an appropriate reference solution or benchmark.

## Scope and limitations

Current limitations include:

- a two-dimensional idealized configuration;
- model parameters defined directly in MATLAB source files;
- no packaged automated verification suite;
- no guaranteed compatibility across MATLAB releases;
- no semantic-versioned software release yet; and
- explicit mechanics coupled to a separately solved implicit pressure step rather than a fully explicit multiphysics update.

## Citation

Please use [`CITATION.cff`](CITATION.cff) when citing the software. When the associated IMEX-method article is published, its final bibliographic information should be added here and to the citation file.

Related published formulation:

> Sabah, M., Hofmann, H., Cacace, M., Jalali, M. R., & Kivi, I. R. (2026). Modeling injection-induced seismicity using a fully coupled poroviscoelasto-dynamic extended finite element approach with stabilized contact and rate-and-state friction. *Computers and Geotechnics, 191*, 107803. https://doi.org/10.1016/j.compgeo.2025.107803

## Author

**Mohammad Sabah, PhD**  
Computational geomechanics · induced seismicity · coupled multiphysics · XFEM · rate-and-state friction  
Technische Universität Berlin  
[ORCID](https://orcid.org/0000-0002-7260-6138) · [LinkedIn](https://www.linkedin.com/in/mohammad-sabah) · [GitHub](https://github.com/MohammadSabah93)

## Contributing

Bug reports, reproducibility feedback and focused improvements are welcome. Please read [`CONTRIBUTING.md`](CONTRIBUTING.md) before opening an issue or pull request.

## License

Released under the [MIT License](LICENSE).
