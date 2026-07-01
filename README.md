# Numerical Poisson Solver: 2D Lid-Driven Cavity Flow Simulation

Final Project for the **Physics Hishuvit 4** Course  
**Authors:** Alexey Frolov & Piotr Smirnov



https://github.com/user-attachments/assets/edb58e8f-d9d9-4ef4-b89e-e509d93ec3b4



---

## Project Overview
This project implements a numerical simulator for incompressible viscous fluid dynamics inside a closed square domain. It models the classic **Lid-Driven Cavity** hydrodynamic benchmark: the top fluid layer is driven by a constant shear velocity $U$, dragging the fluid beneath it due to viscous forces, while the other three boundaries remain rigid and stationary.

The simulation resolves the 2D Navier-Stokes equations formulated in the **Stream Function ($\psi$) — Vorticity ($\omega$)** variables.

---

## Mathematical Formulation & Physics

The system is governed by two coupled differential equations:

1. **Poisson Equation for the Stream Function:**
   $$\nabla^2 \psi = -\omega$$
   Where velocity components are coupled to the stream function via: $u = \frac{\partial \psi}{\partial y}$ and $v = -\frac{\partial \psi}{\partial x}$.

2. **Vorticity Transport Equation:**
   $$\frac{\partial \omega}{\partial t} + u\frac{\partial \omega}{\partial x} + v\frac{\partial \omega}{\partial y} = \nu \nabla^2 \omega$$
   Where $\nu$ represents the kinematic viscosity of the fluid.

### Numerical Methodology
* **Discretization:** Finite Difference Method (FDM) on a uniform Cartesian grid.
* **Linear Solver:** The elliptical Poisson system and the vorticity evolution step are solved numerically using the iterative **Gauss-Seidel** relaxation scheme until a convergence tolerance of $\epsilon < 10^{-5}$ is reached.

---

## Libraries used
The project is built entirely in **Python 3** using a clean modular architecture:
* `NumPy` — Matrix operations and grid processing.
* `Matplotlib` — Renders streamline contours (`contourf`) and velocity vector fields (`quiver`).
* `Imageio` — Compiles serialized frame sequences into the final `.mp4` video.

---

## Project Structure

* `main.py` — Core simulation orchestration script. Handles parameter definition, time-stepping loops, plotting pipelines, and `.mp4` compilation.
* `solver.py` — High-performance computational kernel containing:
  * `solve_poisson_gs`: Iterative Gauss-Seidel Poisson equation solver.
  * `solve_vorticity_gs`: Vorticity transport processing matrix.
  * `apply_boundary_conditions_psi` / `_omega`: Implementation of Dirichlet boundary conditions (no-slip/no-penetration constraints) and Thom's formula for boundary vorticity updates under a moving boundary plate.

---

## Installation & Usage

1. Install required dependencies:
   ```bash
   pip install numpy matplotlib imageio imageio-ffmpeg
