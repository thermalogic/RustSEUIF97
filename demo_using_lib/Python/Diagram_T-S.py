"""
The Python example to call the seuif97

T-s Diagram

1 Calculating isoenthalpic lines isoh(200, 3600)kJ/kg
2 Calculating isobar lines  isop(611.657e-6,100)MPa
3 Calculating saturation lines x=0,x=1
4 Calculating isoquality lines x(0.1,0.9)

Works with both local build and PyPI install:

Local build:
  1. cargo build -r --features python
  2. rename the shared library (*.dll / *.so) to seuif97.pyd
  3. add the path of seuif97.pyd to sys.path before importing

PyPI install:
  pip install seuif97

Author: Cheng Maohua
Email: cmh@seu.edu.cn

"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# add the path of seuif97.pyd to sys.path
sys.path.insert(0, os.path.abspath(r'D:\sim_rankine_workbench\iapws-if97\RustSEUIF97\target\release'))

from seuif97 import ph, tx, px

OT = 1
OS = 5

xAxis = "s"
yAxis = "T"
title = {"T": "T, °C", "s": "s, kJ/kgK"}

plt.title("%s-%s Diagram" % (yAxis, xAxis))
plt.xlabel(title[xAxis])
plt.ylabel(title[yAxis])
plt.grid()

Pt = 611.657e-6
isoh = np.linspace(100, 3600, 27)
isop = np.array([Pt, 0.001, 0.002, 0.004, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5,
                 1.0, 2.0, 5.0, 10.0, 15, 20.0, 22, 25, 28, 30, 35, 40, 45, 50.0, 60, 80, 100.0])

# Calculating isoenthalpic lines isoh(200, 3600)kJ/kg
for h in isoh:
    T = np.array([ph(p, h, OT) for p in isop])
    S = np.array([ph(p, h, OS) for p in isop])
    plt.plot(S, T, 'g', lw=0.5)

# Calculating isobar lines  isop(611.657e-6,100)MPa
for p in isop:
    T = np.array([ph(p, h, OT) for h in isoh])
    S = np.array([ph(p, h, OS) for h in isoh])
    plt.plot(S, T, 'b', lw=0.5)

# saturate water to wet steam
for p in [Pt, 0.001, 0.002]:
    s = ph(p, 100, OS)
    t = ph(p, 100, OT)
    s_sat_water = px(p, 0.0, OS)
    t_sat_water = px(p, 0.0, OT)
    plt.plot([s_sat_water, s], [t_sat_water, t], 'b', lw=0.5)

tc = 647.096 - 273.15
T = np.linspace(0.11, tc, 100)

for x in np.array([0, 1.0]):
    S = np.array([tx(t, x, OS) for t in T])
    plt.plot(S, T, 'r', lw=1.0)

for x in np.linspace(0.1, 0.9, 11):
    S = np.array([tx(t, x, OS) for t in T])
    plt.plot(S, T, 'r--', lw=0.5)

plt.show()
