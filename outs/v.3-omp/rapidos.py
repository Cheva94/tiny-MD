import numpy as np
import pandas as pd
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

################################################################################
############################# GFLOPS
################################################################################

#############################

for h in [7, 14, 21, 28]:

    fig, ax = plt.subplots(figsize=(12.5,10))

    ax.set_title(f'# Hilos = {h}')

    h_h = pd.read_csv(f'/home/cheva/CP/tiny_md/outs/v.3-omp/rapidos_hilos-{h}.csv', comment = '#')
    h_hG = h_h.groupby(['N','BLOCK'])

    for part in h_h.iloc[:, 0].drop_duplicates().to_numpy():
        x = h_hG.describe().loc[part].index.to_numpy()
        y = h_hG.describe().loc[part].to_numpy()[:,1]
        yerr = h_hG.describe().loc[part].to_numpy()[:,2]

        ax.errorbar(x, y, yerr=yerr, capsize=8, elinewidth=5, capthick = 2,
                    linestyle='--', fmt='o', markersize=10, label=f'# Partículas={part}')

    ax.set_ylabel('GFLOPS')
    ax.set_xlabel('Tamaño de bloques')
    ax.legend(loc='upper left', fontsize=10)
    ax.set_xticks([64, 128, 256])

    plt.tight_layout()

    plt.savefig(f'/home/cheva/CP/tiny_md/outs/v.3-omp/rapidos_hilos-{h}.png')

plt.show()
