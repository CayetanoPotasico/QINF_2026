import pandas as pd
import matplotlib.pyplot as plt

# === Configuración general ===
plt.rcParams.update({
    "text.usetex": True,
    "font.size": 12,
    "font.family": "serif",
    "grid.linestyle": "--",
    "grid.alpha": 0.5
})

# Cargar datos
df = pd.read_csv('cbp_ej10_entropias_S15.dat', sep=r'\s+', comment='#', 
                 names=['Gamma', 'Simulacion', 'S15'])

# Seleccionar 9 valores representativos de Gamma para el panel 3x3
gammas_to_plot = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0]

fig, axes = plt.subplots(3, 3, figsize=(12, 10), sharex=True, sharey=True)
axes = axes.flatten()

max_entropy = 1

color_hist = "#E66101"
color_line = "black"

for i, g in enumerate(gammas_to_plot):
    ax = axes[i]
    
    # Filtramos los datos solo para el Gamma actual
    data_g = df[df['Gamma'] == g]['S15']
    
    # Calculamos media y std para este Gamma
    mean_val = data_g.mean()
    std_val = data_g.std()
    
    # Dibujamos el histograma
    ax.hist(data_g, bins=60, range=(0, max_entropy), density=True, 
            color=color_hist, alpha=0.7, edgecolor='white', linewidth=0.5)
    
    # Línea vertical para marcar la media
    ax.axvline(mean_val, color=color_line, linestyle='dashed', linewidth=1.5, 
               label=f'$\\mu = {mean_val:.2f}$\n$\\sigma = {std_val:.2f}$')
    
    ax.set_title(f'$\\Gamma = {g}$', fontsize=13, fontweight='bold')
    ax.grid(True)
    ax.legend(loc='upper right', fontsize=10, handlelength=1)
    
    # Ajustar ejes
    ax.set_xlim(0, max_entropy)
    
    # Solo poner etiquetas en los bordes
    if i >= 6:
        ax.set_xlabel(r'Entropía $\langle S_{15} \rangle$', fontsize=12)
    if i % 3 == 0:
        ax.set_ylabel(r'Frecuencia (Densidad)', fontsize=12)

plt.suptitle('Distribución Estadística de la Entropía de Entrelazamiento', fontsize=16, y=0.98)
plt.tight_layout()

plt.savefig('panel_histogramas.pdf', dpi=300, bbox_inches='tight')
print("Done!")