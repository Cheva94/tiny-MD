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

################################################################################
############################# GFLOPS
################################################################################

N = pd.read_csv('/home/cheva/CP/tiny_md/outs/v.3-omp/tiempos.csv').iloc[:, 0].to_numpy()
Ts = pd.read_csv('/home/cheva/CP/tiny_md/outs/v.3-omp/tiempos.csv').iloc[:, 1].to_numpy()
Tp1 = pd.read_csv('/home/cheva/CP/tiny_md/outs/v.3-omp/tiempos.csv').iloc[:, 2].to_numpy()
Tp7 = pd.read_csv('/home/cheva/CP/tiny_md/outs/v.3-omp/tiempos.csv').iloc[:, 3].to_numpy()
Tp14 = pd.read_csv('/home/cheva/CP/tiny_md/outs/v.3-omp/tiempos.csv').iloc[:, 4].to_numpy()
Tp21 = pd.read_csv('/home/cheva/CP/tiny_md/outs/v.3-omp/tiempos.csv').iloc[:, 5].to_numpy()
Tp28 = pd.read_csv('/home/cheva/CP/tiny_md/outs/v.3-omp/tiempos.csv').iloc[:, 6].to_numpy()

Ac1 = Ts / Tp1
Ac7 = Ts / Tp7
Ac14 = Ts / Tp14
Ac21 = Ts / Tp21
Ac28 = Ts / Tp28

E1 = Ac1 / 1
E7 = Ac7 / 7
E14 = Ac14 / 14
E21 = Ac21 / 21
E28 = Ac28 / 28

print(E21)
# fig, ax = plt.subplots(figsize=(12.5,10))
#
# ax.set_title('Aceleración')
#
# # ax.plot(N, Ac1, linestyle='--', marker='o', markersize=10, label='# Hilos = 1')
# ax.plot(N, Ac7, linestyle='--', marker='o', markersize=10, label='# Hilos = 7')
# ax.plot(N, Ac14, linestyle='--', marker='o', markersize=10, label='# Hilos = 14')
# ax.plot(N, Ac21, linestyle='--', marker='o', markersize=10, label='# Hilos = 21')
# ax.plot(N, Ac28, linestyle='--', marker='o', markersize=10, label='# Hilos = 28')
#
# ax.set_ylabel('Sp = Ts/Tp')
# ax.set_xlabel('Cantidad de partículas')
# ax.legend(loc='upper left', fontsize=16)
# ax.set_xscale('log')
# ax.set_xticks(N)
# ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
#
# plt.tight_layout()
#
# plt.savefig('/home/cheva/CP/tiny_md/outs/v.3-omp/aceleracion.png')
#
#
#
#
# fig, ax = plt.subplots(figsize=(12.5,10))
#
# ax.set_title('Eficiencia')
#
# # ax.plot(N, E1, linestyle='--', marker='o', markersize=10, label='# Hilos = 1')
# ax.plot(N, E7, linestyle='--', marker='o', markersize=10, label='# Hilos = 7')
# ax.plot(N, E14, linestyle='--', marker='o', markersize=10, label='# Hilos = 14')
# ax.plot(N, E21, linestyle='--', marker='o', markersize=10, label='# Hilos = 21')
# ax.plot(N, E28, linestyle='--', marker='o', markersize=10, label='# Hilos = 28')
# ax.axhline(0.6, color='black', lw = 3)
#
# ax.set_ylabel('E = Sp/P')
# ax.set_xlabel('Cantidad de partículas')
# ax.legend(loc='upper left', fontsize=16)
# ax.set_xscale('log')
# ax.set_xticks(N)
# ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
#
# plt.tight_layout()
#
# plt.savefig('/home/cheva/CP/tiny_md/outs/v.3-omp/eficiencia.png')
#
# plt.show()
