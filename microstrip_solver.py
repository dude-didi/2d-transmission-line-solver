import numpy as np
from fenics import *
import matplotlib.pyplot as plt


def solve_microstrip(
    width=0.02,
    height=0.012,
    strip_width=0.006,
    strip_thickness=0.0004,
    substrate_thickness=0.0016,
    eps_sub=4.0,
    nx=120,
    ny=120,
):
    """Solve a 2D microstrip electrostatic problem with FEniCS.

    Parameters
    ----------
    width : float
        Total domain width in meters.
    height : float
        Total domain height in meters.
    strip_width : float
        Signal strip width in meters.
    strip_thickness : float
        Signal strip thickness in meters.
    substrate_thickness : float
        Height of dielectric substrate in meters.
    eps_sub : float
        Relative permittivity of the substrate.
    nx, ny : int
        Number of mesh cells along x and y.
    """

    # Create mesh for the full cross-section.
    mesh = RectangleMesh(Point(0.0, 0.0), Point(width, height), nx, ny)
    Vh = FunctionSpace(mesh, "Lagrange", 1)

    # Define unknown and test function.
    u = TrialFunction(Vh)
    v = TestFunction(Vh)
    U = Function(Vh)

    # Relative permittivity distribution:
    # - substrate under the signal plane -> eps_sub
    # - otherwise air -> 1.0
    # Use a simple coordinate-based piecewise expression.
    eps_r = conditional(lt(x[1], substrate_thickness), Constant(eps_sub), Constant(1.0))
    eps0 = 8.8541878128e-12

    # Dirichlet boundaries.
    # Ground plane: bottom edge.
    ground = CompiledSubDomain("near(x[1], 0.0)")

    # Strip conductor: top boundary segment near the center.
    x_left = 0.5 * (width - strip_width)
    x_right = 0.5 * (width + strip_width)
    strip_top = height - strip_thickness / 2.0
    strip = CompiledSubDomain(
        "x[1] > z - tol && x[1] < z + tol && x[0] > x_left && x[0] < x_right",
        z=strip_top,
        x_left=x_left,
        x_right=x_right,
        tol=1e-8,
    )

    bc_ground = DirichletBC(Vh, Constant(0.0), ground)
    bc_strip = DirichletBC(Vh, Constant(1.0), strip)
    bcs = [bc_ground, bc_strip]

    # Weak form for electrostatic DC problem.
    a = eps0 * eps_r * dot(grad(v), grad(u)) * dx
    L = Constant(0.0) * v * dx

    solve(a == L, U, bcs)

    # Estimate electrostatic energy and capacitance per length.
    gradU = grad(U)
    energy_density = 0.5 * eps0 * eps_r * dot(gradU, gradU)
    W = assemble(energy_density * dx)
    V_diff = 1.0  # V = 1 V across the strip to ground
    C = 2.0 * W / (V_diff ** 2)

    # First-order characteristic impedance estimate.
    c0 = 299792458.0
    Z0 = 1.0 / (c0 * C)

    return mesh, U, C, Z0


def plot_solution(mesh, U):
    """Plot the solution and the electric field."""
    fig, ax = plt.subplots(figsize=(8, 5))
    p = plot(U, mode="contourf", cmap="viridis")
    fig.colorbar(p, ax=ax, label="Voltage [V]")
    ax.set_title("2D microstrip voltage distribution")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    plt.tight_layout()
    plt.show()


def main():
    width = 0.02
    height = 0.012
    strip_width = 0.006
    strip_thickness = 4e-4
    substrate_thickness = 1.6e-3
    eps_sub = 4.0

    mesh, U, C, Z0 = solve_microstrip(
        width=width,
        height=height,
        strip_width=strip_width,
        strip_thickness=strip_thickness,
        substrate_thickness=substrate_thickness,
        eps_sub=eps_sub,
        nx=140,
        ny=140,
    )

    print(f"Estimated capacitance per unit length: {C:.6e} F/m")
    print(f"Approximate characteristic impedance: {Z0:.2f} ohms")

    plot_solution(mesh, U)


if __name__ == "__main__":
    main()

