def apply_boundary_conditions_psi(psi):
    psi[0, :] = 0
    psi[-1, :] = 0
    psi[:, 0] = 0
    psi[:, -1] = 0

def solve_poisson_gs(psi, omega, h, max_iter=500, tol=1e-5):
    N = psi.shape[0]
    for it in range(max_iter):
        max_diff = 0.0
        for i in range(1, N - 1):
            for j in range(1, N - 1):
                psi_new = 0.25 * (
                    psi[i+1, j] + psi[i-1, j] +
                    psi[i, j+1] + psi[i, j-1] +
                    h**2 * omega[i, j]
                )
                diff = abs(psi_new - psi[i, j])
                max_diff = max(max_diff, diff)
                psi[i, j] = psi_new
        apply_boundary_conditions_psi(psi)
        if max_diff < tol:
            #print(f"Converged in {it} iterations")
            break

def solve_vorticity_gs(f, psi, omega, viscosity, max_iter=500, tol=1e-5):
    N = f.shape[0]
    for it in range(max_iter):
        max_diff = 0.0
        for i in range(1, N - 1):
            for j in range(1, N - 1):
                f_new = 0.25*(
                    (psi[i+1,j] - psi[i-1,j])*(omega[i,j+1] - omega[i,j-1]) -
                    (psi[i,j+1] - psi[i,j-1])*(omega[i+1,j] - omega[i-1,j])
                ) + viscosity*(
                    omega[i-1,j] + omega[i,j-1] + omega[i,j+1] + omega[i+1,j] - 4*omega[i,j]
                )
                diff = abs(f_new - f[i, j])
                max_diff = max(max_diff, diff)
                f[i, j] = f_new
        if max_diff < tol:
            #print(f"Converged in {it} iterations")
            break

def apply_boundary_conditions_omega(omega, psi, speed, h):
    N = omega.shape[0]
    for i in range(1,N-1):
        omega[0,i] = -(2/h**2)*psi[N-2,i]
        omega[i,0] = -(2/h**2)*psi[i,1]
        omega[i,N-1] = -(2/h**2)*psi[i,N-2]
        omega[N-1,i]= -(2/h**2)*psi[1,i] - (speed/h)