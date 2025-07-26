import numpy as np
import matplotlib.pyplot as plt

# Carregar trajetória salva pelo Fortran
traj = np.loadtxt("trajetoria.txt")
x_traj, y_traj = traj[:,0], traj[:,1]

# Definir a função f(x, y)
def f(x, y):
    return 4*x + 2*y + x**2 - 2*x**4 + 2*x*y - 3*y**2

# Gerar malha para curvas de nível
x = np.linspace(-2, 2, 400)
y = np.linspace(-2, 2, 400)
X, Y = np.meshgrid(x, y)
Z = f(X, Y)

# Plotar
plt.figure(figsize=(8,6))
cs = plt.contour(X, Y, Z, levels=40, cmap='viridis')
plt.clabel(cs, inline=True, fontsize=8)
plt.plot(x_traj, y_traj, 'r-o', label="Trajetória")
plt.title("Curvas de Nível e Trajetória de busca")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid(True)
plt.show()
