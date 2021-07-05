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
############################# Energía potencial
################################################################################

fig, ax = plt.subplots(figsize=(12.5,10))

densidad = [1.2,1.1,1,0.9,0.8,0.7,0.6,0.5,0.4]

#,N=256
Lab2 = [-1406.703522,-1518.2026,-1293.801747,-1294.946608,-1216.611033,-1100.691575,-959.598585,-807.563515,-645.203302]
Lab3 = [-1411.993196, -1524.492037, -1294.790960, -1298.679460, -1216.773334, -1112.284644, -969.776740, -809.027684, -652.913989]
ax.plot(densidad, Lab2, label='N=256 (Lab2)')
ax.plot(densidad, Lab3, label='N=256 (Lab3)', ls=':')

#,N=500
Lab2 = [-2746.993649,-2972.078716,-2521.888269,-2521.636934,-2377.84898,-2153.604974,-1876.028795,-1575.063486,-1265.668041]
Lab3 = [-2768.567416, -2990.782757, -2588.875079, -2541.078084, -2392.264389, -2166.329516, -1893.392950, -1584.471826, -1275.337347]
ax.plot(densidad, Lab2, label='N=500 (Lab2)')
ax.plot(densidad, Lab3, label='N=500 (Lab3)', ls=':')

#,N=864
Lab2 = [-4751.367838,-5128.308214,-4547.031085,-4358.822304,-4106.70241,-3718.763525,-3250.484801,-2720.060973,-2189.192933]
Lab3 = [-4805.591865, -5185.046842, -4427.650950, -4406.172567, -4150.879872, -3757.022419, -3270.309410, -2754.602653, -2207.630661]
ax.plot(densidad, Lab2, label='N=864 (Lab2)')
ax.plot(densidad, Lab3, label='N=864 (Lab3)', ls=':')

ax.set_ylabel('Epot')
ax.set_xlabel('Densidad')
# ax.set_ylim(3.5,6.5)
ax.legend(loc='lower left', fontsize=15)

plt.tight_layout()

plt.savefig('/home/cheva/CP/tiny_md/outs/v.3-omp/potencial')

plt.show()
