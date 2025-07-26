import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('mapa_temperatura.dat')  # temperaturas em K
data_celsius = data - 273.15                # converter para °C

x = np.linspace(0, 0.01, 50)  # eixo espacial
t_steps = data.shape[0]
time = np.linspace(60, 0, t_steps)  # eixo temporal

plt.figure(figsize=(8,6))
plt.imshow(data_celsius, aspect='auto', extent=[0, 0.01,0, 60, 0], origin='lower' ,cmap='inferno')
plt.colorbar(label='Temperatura (°C)')
plt.xlabel('Posição (m)')
plt.ylabel('Tempo (s)')
plt.title('Distribuição da Temperatura no Tempo e Espaço (°C)')
plt.show()
