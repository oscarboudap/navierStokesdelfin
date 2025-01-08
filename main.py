import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Parámetros del dominio
nx, ny = 200, 100
lx, ly = 10.0, 5.0
dx, dy = lx / nx, ly / ny
dt = 0.001  # Paso temporal reducido
nu = 1e-3  # Viscosidad cinemática

# Campos iniciales
u = np.ones((nx, ny)) * 0.1  # Flujo inicial en x
v = np.zeros((nx, ny))  # Flujo inicial en y
p = np.zeros((nx, ny))  # Presión inicial
u_new, v_new, p_new = np.copy(u), np.copy(v), np.copy(p)

# Delfín como obstáculo
delfin_mask = np.zeros((nx, ny))
for i in range(nx):
    for j in range(ny):
        if (i - nx // 4) ** 2 + (j - ny // 2) ** 2 < (ny // 6) ** 2:
            delfin_mask[i, j] = 1

# Función para resolver Poisson
def solve_poisson(p, b, tol=1e-4, max_iter=500):
    for _ in range(max_iter):
        p_new = np.copy(p)
        p_new[1:-1, 1:-1] = (
            (b[1:-1, 1:-1] +
             (p[2:, 1:-1] + p[:-2, 1:-1]) / dx**2 +
             (p[1:-1, 2:] + p[1:-1, :-2]) / dy**2)
            / (2 / dx**2 + 2 / dy**2)
        )
        if np.linalg.norm(p_new - p) < tol:
            break
        p = np.copy(p_new)
    return p

# Función para actualizar campos
def update_fields():
    global u, v, p

    # Cálculo de términos convectivos
    u_conv = u[1:-1, 1:-1] * (u[2:, 1:-1] - u[:-2, 1:-1]) / (2 * dx) + \
             v[1:-1, 1:-1] * (u[1:-1, 2:] - u[1:-1, :-2]) / (2 * dy)
    v_conv = u[1:-1, 1:-1] * (v[2:, 1:-1] - v[:-2, 1:-1]) / (2 * dx) + \
             v[1:-1, 1:-1] * (v[1:-1, 2:] - v[1:-1, :-2]) / (2 * dy)

    # Cálculo de términos viscosos
    u_visc = nu * ((u[2:, 1:-1] - 2 * u[1:-1, 1:-1] + u[:-2, 1:-1]) / dx**2 +
                   (u[1:-1, 2:] - 2 * u[1:-1, 1:-1] + u[1:-1, :-2]) / dy**2)
    v_visc = nu * ((v[2:, 1:-1] - 2 * v[1:-1, 1:-1] + v[:-2, 1:-1]) / dx**2 +
                   (v[1:-1, 2:] - 2 * v[1:-1, 1:-1] + v[1:-1, :-2]) / dy**2)

    # Cálculo de presión
    b = np.zeros_like(p)
    b[1:-1, 1:-1] = (1 / dt) * ((u[2:, 1:-1] - u[:-2, 1:-1]) / (2 * dx) +
                                (v[1:-1, 2:] - v[1:-1, :-2]) / (2 * dy))
    p = solve_poisson(p, b)

    # Actualización de velocidades
    u_new[1:-1, 1:-1] = u[1:-1, 1:-1] - dt * u_conv + dt * u_visc - \
                        dt * (p[2:, 1:-1] - p[:-2, 1:-1]) / (2 * dx)
    v_new[1:-1, 1:-1] = v[1:-1, 1:-1] - dt * v_conv + dt * v_visc - \
                        dt * (p[1:-1, 2:] - p[1:-1, :-2]) / (2 * dy)

    # Condiciones de frontera
    u_new[0, :], u_new[-1, :], u_new[:, 0], u_new[:, -1] = 0, 0, 1, 1
    v_new[0, :], v_new[-1, :], v_new[:, 0], v_new[:, -1] = 0, 0, 0, 0

    # Aplicar máscara del delfín
    u_new[delfin_mask == 1] = 0
    v_new[delfin_mask == 1] = 0

    # Actualizar campos
    u[:], v[:] = u_new[:], v_new[:]

# Crear animación
fig, ax = plt.subplots(figsize=(8, 4))
ax.set_xlim(0, lx)
ax.set_ylim(0, ly)
cbar = ax.quiver(np.linspace(0, lx, nx), np.linspace(0, ly, ny),
                 u.T, v.T, scale=1, scale_units='xy')

def animate(frame):
    update_fields()
    cbar.set_UVC(u.T, v.T)
    return cbar,

ani = animation.FuncAnimation(fig, animate, frames=200, interval=50)
ani.save("flujo_delfin_corregidoAA.gif", writer="imagemagick")
plt.show()
