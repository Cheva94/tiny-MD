import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from cycler import cycler

# Paleta de colores (deben ser hex)
V = '#00a189'; N = '#fa6200'; R = '#ed3b3b'; Az = '#5ca2f7'; Am = '#ebe842'
Az_claro = '#97C4FA'; G = '#2F4F4F'; S = '#C0C0C0'; S_osc = '#626262'

# Parámetro de figuras
plt.rc('font', weight = 'bold', size = 30)
plt.rc('axes', labelweight = "bold", lw = 3)
plt.rc('lines', lw = 1, ls = '-')
plt.rc('xtick.major', size = 10, width = 3)
plt.rc('ytick.major', size = 10, width = 3)
# plt.rcParams['axes.prop_cycle'] = cycler(color=[V, N, Az, R, Am, G])

Secuencial = pd.read_csv(f'/home/cheva/tiny_md/outs/v.1-seq/part_D.csv', comment = '#').groupby(['N']).describe()
ISPC = pd.read_csv(f'/home/cheva/tiny_md/outs/v.2-vec/v.2.F/particulas_256-6912.csv', comment = '#').groupby(['N']).describe()
OMP = pd.read_csv(f'/home/cheva/tiny_md/outs/v.3-omp/rapidos_hilos-{21}.csv', comment = '#').groupby(['BLOCK', 'N']).describe().loc[128]
CUDA = pd.read_csv(f'/home/cheva/tiny_md/outs/v.4-cuda/NvsB.csv', comment = '#').groupby(['BLOCK', 'N']).describe().loc[32]

fig, ax = plt.subplots(figsize=(12.5,10))

x = Secuencial.index.to_numpy()
y = Secuencial.to_numpy()[:,1]
yerr = Secuencial.to_numpy()[:,2]

ax.errorbar(x, y, yerr=yerr, capsize=8, elinewidth=5, capthick = 2,
            linestyle='--', fmt='o', markersize=10, label=f'Secuencial')

x = ISPC.index.to_numpy()
y = ISPC.to_numpy()[:,1]
yerr = ISPC.to_numpy()[:,2]

ax.errorbar(x, y, yerr=yerr, capsize=8, elinewidth=5, capthick = 2,
            linestyle='--', fmt='o', markersize=10, label=f'ISPC')

x = OMP.index.to_numpy()
y = OMP.to_numpy()[:,1]
yerr = OMP.to_numpy()[:,2]

ax.errorbar(x, y, yerr=yerr, capsize=8, elinewidth=5, capthick = 2,
            linestyle='--', fmt='o', markersize=10, label=f'OMP')

x = CUDA.index.to_numpy()
y = CUDA.to_numpy()[:,1]
yerr = CUDA.to_numpy()[:,2]

ax.errorbar(x, y, yerr=yerr, capsize=8, elinewidth=5, capthick = 2,
            linestyle='--', fmt='o', markersize=10, label=f'CUDA')


ax.set_ylabel('GFLOPS')
ax.set_xlabel('Cantidad de partículas')
ax.legend(loc='upper left', fontsize=16)
ax.set_xscale('log')
ax.set_xticks([256, 500, 864, 1372, 2048, 2916, 4000, 5324, 6912, 8788, 10976, 13500])#, 16384, 19652, 23328, 27436])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.xticks(rotation=60)

plt.tight_layout()

plt.savefig(f'/home/cheva/tiny_md/outs/v.4-cuda/interlab')

plt.show()
