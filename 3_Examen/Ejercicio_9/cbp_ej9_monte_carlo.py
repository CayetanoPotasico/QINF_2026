import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import os



# === Configuración general ===
plt.rcParams.update({
    "text.usetex": True,        # Usa LaTeX real
    "font.size": 12,
    # "font.weight": "bold",       # Hace el texto más grueso
    "font.family": "serif",
    # "text.latex.preamble": r"\usepackage{lmodern}\usepackage{bm}",
    # "axes.grid": True,
    "grid.linestyle": "--",
    "grid.alpha": 0.6
})


input_dir = "out_monte_carlo"
output_dir = "plots_monte_carlo"
os.makedirs(output_dir, exist_ok=True)

gammas = [0.1, 1.0, 10.0]
gamma_labels = [r"$\gamma = 0.1$ (Débil)", r"$\gamma = 1.0$ (Crítico)", r"$\gamma = 10.0$ (Fuerte)"]
states = [
    r"Polo Norte $|0\rangle$ $[0,0,1]$",
    r"Estado $|-\rangle$ $[-1,0,0]$",
    r"Estado $|+\rangle$ $[1,0,0]$"
]

colors = {'x': '#E66101', 'y': '#5E3C99', 'z': '#FDB863'}
state_colors = ['#1b9e77', '#d95f02', '#7570b3']

# ==============================================================================
# 1: PANEL GENERAl 3x3 (Evolución de x, y, z vs Tiempo)
# ==============================================================================
fig, axes = plt.subplots(3, 3, figsize=(15, 11), sharex='col')

for i in range(3):     # Bucle para filas (Estados Iniciales)
    for j in range(3): # Bucle para columnas (Valores de Gamma)

        g_idx = i + 1
        s_idx = j + 1
        file_id = g_idx*10 + s_idx # Formato 11, 12, 13...
        file_path = os.path.join(input_dir, f"cbp_ej9_monte_carlo_{file_id}.dat")
        
        if os.path.exists(file_path):
            
            df = pd.read_csv(file_path, sep=r'\s+', comment='#', names=['t', 'x', 'y', 'z'])

            ax = axes[j, i] # Fila = Estado, Columna = Gamma
            ax.plot(df['t'], df['x'], label='$x(t)$', color=colors['x'], lw=2)
            ax.plot(df['t'], df['y'], label='$y(t)$', color=colors['y'], lw=2)
            ax.plot(df['t'], df['z'], label='$z(t)$', color=colors['z'], lw=1.5, linestyle='--')
            
            ax.grid(True, linestyle=':', alpha=0.6)
            ax.set_ylim(-1.05, 1.05)
            
            # Títulos de las columnas (solo arriba)
            if j == 0:
                ax.set_title(gamma_labels[i], fontsize=12, fontweight='bold', pad=10)
            # Etiquetas de las filas (solo a la izquierda)
            if i == 0:
                ax.set_ylabel(states[j], fontsize=11, fontweight='bold', labelpad=10)
            # Etiqueta del eje X (solo abajo)
            if j == 2:
                ax.set_xlabel("Tiempo ($t$)", fontsize=10)
                
            if i == 2 and j == 0:
                ax.legend(loc='upper right', frameon=True, shadow=False)
        else:
            print(f"[-] No se pudo encontrar el archivo: {file_path}")

plt.tight_layout()
grid_out = os.path.join(output_dir, "grid_evolucion_temporal_monte_carlo.pdf")
plt.savefig(grid_out, format='pdf', bbox_inches='tight')
plt.close()
print(f"[+] Guardado panel temporal en: {grid_out}")


# ==============================================================================
# 2: TRAYECTORIAS 3D EN LA ESFERA DE BLOCH (Un gráfico por cada Gamma)
# ==============================================================================
# Esfera unitaria de fondo
u = np.linspace(0, 2 * np.pi, 30)
v = np.linspace(0, np.pi, 30)
x_sphere = np.outer(np.cos(u), np.sin(v))
y_sphere = np.outer(np.sin(u), np.sin(v))
z_sphere = np.outer(np.ones(np.size(u)), np.cos(v))

for i in range(3):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Dibujamos la esfera de Bloch translúcida de fondo
    ax.plot_wireframe(x_sphere, y_sphere, z_sphere, color='#999999', alpha=0.15, linewidth=0.5)
    
    # Dibujamos los ejes cardinales principales internos
    ax.plot([-1, 1], [0, 0], [0, 0], color='k', linestyle=':', alpha=0.4)
    ax.plot([0, 0], [-1, 1], [0, 0], color='k', linestyle=':', alpha=0.4)
    ax.plot([0, 0], [0, 0], [-1, 1], color='k', linestyle=':', alpha=0.4)
    
    # Marcar explícitamente el punto estacionario final
    x_ss = -gammas[i]**2 / (gammas[i]**2 + 8)
    y_ss = 4 / gammas[i] * x_ss
    z_ss = 0
    ax.scatter(x_ss, y_ss, z_ss, color='black', s=30, label='Estado Estacionario', zorder=10)

    g_idx = i + 1
    for j in [2,1,0]: # Metemos los 3 estados iniciales en el mismo espacio 3D
        s_idx = j + 1
        file_id = g_idx*10 + s_idx
        file_path = os.path.join(input_dir, f"cbp_ej9_monte_carlo_{file_id}.dat")
        
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, sep=r'\s+', comment='#', names=['t', 'x', 'y', 'z'])

            # Trayectoria continua
            ax.plot(df['x'], df['y'], df['z'], color=state_colors[j], lw=2.5, label=f"Inicio: {states[j].split()[1]}")
            # Punto de inicio resaltado
            ax.scatter([df['x'].iloc[0]], [df['y'].iloc[0]], [df['z'].iloc[0]], color=state_colors[j], s=30)

    ax.set_xlim([-1.05, 1.05])
    ax.set_ylim([-1.05, 1.05])
    ax.set_zlim([-1.05, 1.05])
    ax.set_xlabel("X", fontsize=10)
    ax.set_ylabel("Y", fontsize=10)
    ax.set_zlabel("Z", fontsize=10)
    ax.set_title(f"Trayectorias en la Esfera de Bloch\n{gamma_labels[i]}", fontsize=12, fontweight='bold', pad=15)
    ax.legend(loc='lower left', bbox_to_anchor=(0.0, 0.75), frameon=True, fontsize=9)
    
    # Ajuste de perspectiva
    ax.view_init(elev=20, azim=45)

    fig_3d_out = os.path.join(output_dir, f"bloch_3d_gamma_monte_carlo_{g_idx}.pdf")
    plt.savefig(fig_3d_out, format='pdf', dpi=300)
    plt.show()
    plt.close()
    print(f"[+] Guardado gráfico 3D en: {fig_3d_out}")