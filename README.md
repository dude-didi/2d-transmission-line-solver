# 2D Transmission Line Solver

This project provides a compact finite-element solver for a 2D microstrip transmission line cross-section using FEniCS.

The goal is to compute a rough but useful approximation of the line capacitance and characteristic impedance from a cross-sectional electrostatic field solution.

## Features

- 2D microstrip cross-section geometry
- FEniCS-based finite-element solve
- Substrate and air dielectric regions
- Dirichlet boundary conditions for conductor and ground plane
- Capacitance and approximate characteristic impedance estimation
- Plot of voltage field and electric field

## Problem definition

The solver treats the cross-section of a microstrip line as a 2D electrostatic problem:

- ground plane at the bottom boundary: V = 0
- signal conductor strip at the top boundary: V = 1
- dielectric regions: air + substrate
- solve Laplace's equation in the cross-section:

  ∇ · (ε ∇V) = 0

Then it estimates capacitance from the stored electric energy:

- W = 1/2 ∫ ε|∇V|² dA
- C = 2W / V²
- Z0 ≈ 1 / (c C) for a first-order estimate

## Installation

Install FEniCS and the plotting dependencies:

```bash
pip install -r requirements.txt
```

If your environment is a special FEniCS distribution (for example, conda-based), use the installation method recommended for that system instead.

## Run

```bash
python microstrip_solver.py
```

This script will:

1. build the 2D mesh
2. solve the electrostatic boundary-value problem
3. compute a capacitance estimate
4. print the characteristic impedance estimate
5. plot the voltage field

## Notes

This is an educational and engineering-oriented starting point, not a full wave solver. It is appropriate for:

- first-pass microstrip design
- teaching and prototyping
- benchmarking a finite-element setup
- building a foundation for more advanced full-wave and time-domain solvers

## Next steps

Possible extensions include:

- Gmsh-based geometry and mesh generation
- stripline and coaxial cross-sections
- loss and conductivity terms
- frequency-dependent material models
- extraction of per-unit-length parameters for an actual transmission line model

