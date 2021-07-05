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

F = pd.read_csv(f'/home/cheva/tiny_md/outs/v.4-cuda/NvsB.csv', comment = '#').groupby(['BLOCK', 'N'])

F = F.describe()

fig, ax = plt.subplots(figsize=(12.5,10))

for B in [32,64,128,256]:

    x = F.loc[B].index.to_numpy()
    y = F.loc[B].to_numpy()[:,1]
    yerr = F.loc[B].to_numpy()[:,2]

    ax.errorbar(x, y, yerr=yerr, capsize=8, elinewidth=5, capthick = 2,
                linestyle='--', fmt='o', markersize=10, label=f'Tamaño de bloque {B}')


    ax.set_ylabel('GFLOPS')
    ax.set_xlabel('Cantidad de partículas')
    ax.legend(loc='upper left', fontsize=16)
    ax.set_xscale('log')
    ax.set_xticks([256, 500, 864, 1372, 2048, 2916])#, 4000, 5324, 6912, 8788, 10976, 13500])#, 16384, 19652, 23328, 27436])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xticks(rotation=60)

    plt.tight_layout()

    plt.savefig(f'/home/cheva/tiny_md/outs/v.4-cuda/NvsB')

plt.show()
