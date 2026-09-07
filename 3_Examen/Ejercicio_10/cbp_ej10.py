import pandas as pd
import matplotlib.pyplot as plt

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

# Cargar datos
df = pd.read_csv('cbp_ej10_entropias_S15.dat', sep=r'\s+', comment='#', 
                 names=['Gamma', 'Simulacion', 'S15'])

# Agrupar por Gamma para sacar la media y la desviación estándar
stats = df.groupby('Gamma')['S15'].agg(['mean', 'std']).reset_index()

print(stats)

plt.figure(figsize=(10, 6))
plt.plot(stats['Gamma'], stats['mean'], '.-', color="#E63301", label=r'$\langle S_{15} \rangle \pm 1\sigma$')
plt.fill_between(stats['Gamma'], stats['mean'] - stats['std'], stats['mean'] + stats['std'], color="#E66101", alpha=0.2)

plt.title('Entropía Media de Entrelazamiento de los 4 primeros qubits', fontsize=14, fontweight='bold')
plt.xlabel(r'$\Gamma$', fontsize=12)
plt.ylabel(r'$\langle S_{15} \rangle$', fontsize=12)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()
plt.tight_layout()

# Guardar y mostrar
plt.tight_layout()
plt.savefig('entropia_media.pdf', dpi=300)
print("Done!")