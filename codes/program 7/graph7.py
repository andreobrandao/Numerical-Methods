import matplotlib.pyplot as plt
import numpy as np

# Substitua pelo nome do seu arquivo
arquivo = "resultados_temp.txt"

# Leitura dos dados
n_vals = []
trap_vals = []
simp13_vals = []
simp38_vals = []

with open(arquivo, 'r') as f:
    next(f)  # pula cabeçalho
    for linha in f:
        dados = linha.strip().split()
        if len(dados) >= 4:
            n = int(dados[0])
            trap = float(dados[1])
            s13 = float(dados[2])
            s38 = float(dados[3])
            n_vals.append(n)
            trap_vals.append(trap)
            simp13_vals.append(s13)
            simp38_vals.append(s38)

# Converte listas em arrays
n_vals = np.array(n_vals)
trap_vals = np.array(trap_vals)
simp13_vals = np.array(simp13_vals)
simp38_vals = np.array(simp38_vals)

# Plot
plt.figure(figsize=(10, 6))
plt.plot(n_vals, trap_vals, 'o-', label='Trapézio')
plt.plot(n_vals, simp13_vals, 's-', label='Simpson 1/3')
plt.plot(n_vals, simp38_vals, 'd-', label='Simpson 3/8')

plt.xlabel('n (nós em cada direção)')
plt.ylabel('Temperatura Média (°C)')
plt.title('Convergência da Temperatura Média vs. n')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
