# Hybrid IMEX XFEM for Injection-Induced Seismicity

[![MATLAB](https://img.shields.io/badge/MATLAB-source%20code-e16737)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Research software](https://img.shields.io/badge/status-research%20software-4c1)](#scope-and-limitations)

A MATLAB implementation of a hybrid implicit–explicit (IMEX) time-integration strategy for fully coupled hydromechanical XFEM simulation of injection-induced fault slip.

The framework is designed for the severe timescale separation between slow reservoir pressurization and rapid rupture. It combines implicit integration during the quasi-static response with explicit integration when the fault accelerates.

## Capabilities

- two-dimensional coupled poroelastic deformation and fluid flow;
- XFEM representation of an embedded fault;
- rate-and-state friction and nonlinear fault contact;
- implicit–explicit switching based on the evolving fault response;
- inertia and absorbing/damping boundary contributions;
- matrix–fracture hydraulic coupling; and
- post-processing of displacement, pressure, stress, slip, and seismicity measures.

## Quick start

### Requirements

- MATLAB;
- all repository `.m` files available on the MATLAB path.

The code is a research prototype. A specific minimum MATLAB release and toolbox compatibility matrix have not yet been established.

### Run the example

1. Clone or download this repository.
2. Open the repository folder in MATLAB.
3. Review the simulation setup:
   - geometry, mesh, time stepping, and initial conditions in `X_FEM_PoroElastic_Im_Ex.m`;
   - material, hydraulic, dynamic, and frictional properties in `defineModelParameters.m`;
   - mechanical and hydraulic boundary conditions in `defineBoundaryConditions.m`.
4. Run:

```matlab
X_FEM_PoroElastic_Im_Ex
```

The default case uses a 100 m × 100 m domain with an embedded diagonal fault. Because the simulation can transition to very small dynamic time steps, runtime and memory demand depend strongly on the chosen mesh and switching parameters.

## Code map

| Component | Main files |
|---|---|
| Driver and time integration | `X_FEM_PoroElastic_Im_Ex.m` |
| Parameters and boundary conditions | `defineModelParameters.m`, `defineBoundaryConditions.m` |
| Mesh and XFEM enrichment | `MeshGeneration_2D.m`, `levelSet.m`, `enrElem.m` |
| Coupled matrix assembly | `StiffnessMatrix.m`, `MassMatrix.m`, `StorageMatrix.m`, `ConductanceMatrix.m`, `CouplingMatrix.m` |
| Fault contact and friction | `Interface.m`, `Lagrange.m`, `StabLagrange.m`, `updateInterface.m`, `updateInterface_expilict.m` |
| Dynamic boundaries and seismicity | `dampingBoundary.m`, `SeismicityParameters.m` |
| Visualization | `postProcess.m` |

## Scope and limitations

This repository is intended for scientific development and reproducible numerical experimentation, not operational seismic-hazard forecasting.

Current limitations include:

- a two-dimensional idealized configuration;
- input parameters defined directly in MATLAB source files;
- no packaged automated verification suite;
- no guaranteed compatibility across MATLAB releases; and
- no formal software release or semantic version tag yet.

Users should perform independent verification, mesh and time-step convergence studies, and site-specific calibration before interpreting results physically.

## Citation

Please use the repository’s [`CITATION.cff`](CITATION.cff) metadata when citing the software. When the associated IMEX-method article is published, its bibliographic details should be added to both this section and the citation file.

Related published formulation:

> Sabah, M., Hofmann, H., Cacace, M., Jalali, M. R., & Kivi, I. R. (2026). Modeling injection-induced seismicity using a fully coupled poroviscoelasto-dynamic extended finite element approach with stabilized contact and rate-and-state friction. *Computers and Geotechnics, 191*, 107803. https://doi.org/10.1016/j.compgeo.2025.107803

## Contributing

Bug reports, reproducibility feedback, and focused improvements are welcome. Please read [`CONTRIBUTING.md`](CONTRIBUTING.md) before opening an issue or pull request.

## License

Released under the [MIT License](LICENSE).
