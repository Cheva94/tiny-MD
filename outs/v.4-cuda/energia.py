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
plt.rc('lines', lw = 6, ls = '-')
plt.rc('xtick.major', size = 10, width = 3)
plt.rc('ytick.major', size = 10, width = 3)
plt.rcParams['axes.prop_cycle'] = cycler(color=[V, N, Az, R, Am, G])

for N in [256,500,864]:
    fig, ax = plt.subplots(figsize=(12.5,10))

    ax.set_title(f'# partículas = {N}')

    E = pd.read_csv(f'/home/cheva/tiny_md/outs/v.4-cuda/energia_{N}.csv', comment = '#')

    d = E.columns.to_numpy()[1:].astype(float)
    E = E.to_numpy()

    for lab in [0,1,2,3]:
        y = E[lab,1:]
        ax.plot(d, y, label=f'{E[lab,0]}')

        ax.set_ylabel('Epot')
        ax.set_xlabel('Densidad')
        ax.legend(loc='upper right', fontsize=15)

        plt.tight_layout()

        plt.savefig(f'/home/cheva/tiny_md/outs/v.4-cuda/energia_{N}')

plt.show()
