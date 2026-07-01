# Physics Hishuvit 4 final project
# Alexey Frolov and Piotr Smirnov

import math
import matplotlib.pyplot as plt
import imageio
import numpy as np
import os
folder = "frames"

import solver

# //1. Initial variables
N = 15 # Number of squares
T = 7 # Time in seconds
dt = 0.01 
viscosity = 0.001
U = 1.0 # Speed of the top plate

N_steps = round(T/dt)
h = 1/(N-1)
dx = h
dy = h

# //2. Initialization
psi = np.zeros((N, N)) # Stream
omega = np.zeros((N, N)) # Vorticity
f = np.zeros_like(omega) # omega (t+1) = omega(t) + f
u = np.zeros_like(psi) # Velocity x
v = np.zeros_like(psi) # Velocity y

x = np.linspace(0, 1, N)
y = np.linspace(0, 1, N)
X, Y = np.meshgrid(x, y) # For using numpy visualization functions

# //3. Iteration over time
for t in range(N_steps):
    # //Solving Poisson eq and calculating the speeds
    solver.solve_poisson_gs(psi, omega, h) # Imported function for solving Poisson numerically
    u[1:-1, 1:-1] = (psi[2:, 1:-1] - psi[:-2, 1:-1]) / (2 * dy) # u = d(psi)/dy
    v[1:-1, 1:-1] = -(psi[1:-1, 2:] - psi[1:-1, :-2]) / (2 * dx) # v = -d(psi)/dx

    # //Solving vorticity evolution
    solver.solve_vorticity_gs(f, psi, omega, viscosity)
    omega += (dt/h**2)*f
    solver.apply_boundary_conditions_omega(omega, psi, U, h)

    # //Plotting
    plt.figure(figsize=(6, 6))
    plt.contourf(X, Y, psi, levels=50, cmap='viridis')
    plt.quiver(X, Y, u, v, scale=2)
    plt.title(f"t = {round(t * dt, 4)} s")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")

    # Saving the frame
    filename = os.path.join(folder, f"frame_{round(t * dt, 4)}.png")
    plt.savefig(filename, dpi=120)
    plt.close()

    # Printing progress
    print(round(t/N_steps, 3)*100, "%")

# //4. Collecting the frames into a video
os.makedirs(folder, exist_ok=True)
output = "FlowVideo.mp4"
writer = imageio.get_writer(output, fps=round(1/dt))
files = sorted(f for f in os.listdir(folder) if f.endswith(".png"))
for f in files:
    image = imageio.imread(os.path.join(folder, f))
    writer.append_data(image)
writer.close()