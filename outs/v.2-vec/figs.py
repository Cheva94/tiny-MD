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
plt.rc('lines', lw = 8, ls = '-')
plt.rc('xtick.major', size = 10, width = 3)
plt.rc('ytick.major', size = 10, width = 3)
plt.rcParams['axes.prop_cycle'] = cycler(color=[V, N, Az, R, Am, G])

################################################################################
############################# Barrido de N
################################################################################

# Variables

barrido_D = pd.read_csv('/home/cheva/Documents/CP/tiny_md/outs/v.1-seq/part_D.csv', comment = '#').groupby('N')
barrido_J = pd.read_csv('/home/cheva/Documents/CP/tiny_md/outs/v.2-vec/v.2.N/part_J.csv', comment = '#').groupby('N')
barrido_S = pd.read_csv('/home/cheva/Documents/CP/tiny_md/outs/v.2-vec/v.2.F/particulas_SOA.csv', comment = '#').groupby('N')
barrido_SI = pd.read_csv('/home/cheva/Documents/CP/tiny_md/outs/v.2-vec/v.2.F/particulas_256-6912.csv', comment = '#').groupby('N')

prom_D = np.concatenate((barrido_D.mean().to_numpy()[:,0], [np.nan, np.nan, np.nan]))
SD_D = np.concatenate((barrido_D.std().to_numpy()[:,0], [np.nan, np.nan, np.nan]))

prom_J = np.concatenate((barrido_J.mean().to_numpy()[:,0], [np.nan]))
SD_J = np.concatenate((barrido_J.std().to_numpy()[:,0], [np.nan]))

prom_S = np.concatenate((barrido_S.mean().to_numpy()[:,0], [np.nan, np.nan]))
SD_S = np.concatenate((barrido_S.std().to_numpy()[:,0], [np.nan, np.nan]))

prom_SI = barrido_SI.mean().to_numpy()[:,0]
SD_SI = barrido_SI.std().to_numpy()[:,0]

ene = np.array([256, 500, 864, 1372, 2048, 2916, 4000, 5324, 6912])

# print(len(SD_D), len(SD_J), len(SD_S), len(SD_SI))

# Desktop - AoS

fig, ax = plt.subplots(figsize=(12.5,10))

ax.errorbar(ene, prom_D, yerr=SD_D, color=V, capsize=8, elinewidth=5, capthick = 2,
            linestyle='none', fmt='o', markersize=20, label='Desktop - AoS')

ax.set_ylabel('GFLOPS')
ax.set_xlabel('Cantidad de partículas')
ax.set_xscale('log')
ax.set_xlim(10**2,10**4)
ax.set_ylim(3.5,6.5)
ax.legend(loc='lower right')

plt.tight_layout()

plt.savefig('/home/cheva/Documents/CP/tiny_md/outs/Desktop_AoS')

# Jupiterace - AoS

fig, ax = plt.subplots(figsize=(12.5,10))

ax.errorbar(ene, prom_J, yerr=SD_J, color=N, capsize=8, elinewidth=5, capthick = 2,
            linestyle='none', fmt='o', markersize=20, label='Jupiterace - AoS')

ax.set_ylabel('GFLOPS')
ax.set_xlabel('Cantidad de partículas')
ax.set_xscale('log')
ax.set_xlim(10**2,10**4)
ax.set_ylim(17,21.5)
ax.legend(loc='lower right')

plt.tight_layout()

plt.savefig('/home/cheva/Documents/CP/tiny_md/outs/Jupiterace_AoS')

# Jupiterace - AoS vs SoA

fig, ax = plt.subplots(figsize=(12.5,10))

ax.errorbar(ene, prom_J, yerr=SD_J, color=N, capsize=8, elinewidth=5, capthick = 2,
            linestyle='none', fmt='o', markersize=20, label='Jupiterace - AoS')

ax.errorbar(ene, prom_S, yerr=SD_S, color=Az, capsize=8, elinewidth=5, capthick = 2,
            linestyle='none', fmt='o', markersize=20, label='Jupiterace - SoA')


ax.set_ylabel('GFLOPS')
ax.set_xlabel('Cantidad de partículas')
ax.set_xscale('log')
ax.set_xlim(10**2,10**4)
ax.set_ylim(17,22)
ax.legend(loc='upper left', fontsize=20)

plt.tight_layout()

plt.savefig('/home/cheva/Documents/CP/tiny_md/outs/AoSvsSoa')

# Jupiterace - AoS vs SoA con SoA+ISPC

fig, ax = plt.subplots(figsize=(12.5,10))

ax.errorbar(ene, prom_J, yerr=SD_J, color=N, capsize=8, elinewidth=5, capthick = 2,
            linestyle='none', fmt='o', markersize=20, label='Jupiterace - AoS')

ax.errorbar(ene, prom_S, yerr=SD_S, color=Az, capsize=8, elinewidth=5, capthick = 2,
            linestyle='none', fmt='o', markersize=20, label='Jupiterace - SoA')

ax.errorbar(ene, prom_SI, yerr=SD_SI, color=S_osc, capsize=8, elinewidth=5, capthick = 2,
            linestyle='none', fmt='o', markersize=20, label='Jupiterace - SoA+ISPC')

ax.set_ylabel('GFLOPS')
ax.set_xlabel('Cantidad de partículas')
ax.set_xscale('log')
ax.set_xlim(10**2,10**4)
ax.set_ylim(17,40)
ax.legend(loc='upper left', fontsize=20)

plt.tight_layout()

plt.savefig('/home/cheva/Documents/CP/tiny_md/outs/SoA+ISPC')

################################################################################
############################# Energía potencial
################################################################################

fig, ax = plt.subplots(figsize=(12.5,10))

densidad = [1.2,1.1,1,0.9,0.8,0.7,0.6,0.5,0.4]

#,N=256
Original = [-1406.820989,-1517.925785,-1333.553662,-1286.186051,-1216.948914,-1103.896974,-957.025008,-808.991494,-640.473456]
Lab2 = [-1406.703522,-1518.2026,-1293.801747,-1294.946608,-1216.611033,-1100.691575,-959.598585,-807.563515,-645.203302]
ax.plot(densidad, Original, label='N=256 (Orig)')
ax.plot(densidad, Lab2, label='N=256 (Lab2)', ls=':')

#,N=500
Original = [-2746.298803,-2961.935492,-2574.871299,-2518.14288,-2379.308324,-2151.166261,-1879.128806,-1573.684952,-1274.079605]
Lab2 = [-2746.993649,-2972.078716,-2521.888269,-2521.636934,-2377.84898,-2153.604974,-1876.028795,-1575.063486,-1265.668041]
ax.plot(densidad, Original, label='N=500 (Orig)')
ax.plot(densidad, Lab2, label='N=500 (Lab2)', ls=':')

#,N=864
Original = [-4748.551972,-5126.220498,-4367.91629,-4344.931651,-4115.46179,-3719.810983,-3238.702725,-2719.926199,-2195.569164]
Lab2 = [-4751.367838,-5128.308214,-4547.031085,-4358.822304,-4106.70241,-3718.763525,-3250.484801,-2720.060973,-2189.192933]
ax.plot(densidad, Original, label='N=864 (Orig)')
ax.plot(densidad, Lab2, label='N=864 (Lab2)', ls=':')

ax.set_ylabel('Epot')
ax.set_xlabel('Densidad')
# ax.set_ylim(3.5,6.5)
ax.legend(loc='lower left', fontsize=20)

plt.savefig('/home/cheva/Documents/CP/tiny_md/outs/epot')
