"""
The Python example to call the seuif97

p-h (Mollier) Diagram

1. Saturated liquid line (x=0) in blue
2. Saturated vapor line (x=1) in red
3. Isoquality lines (10%~90%) in green dashed
4. Isotherms in black with temperature labels at top

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

from seuif97 import tx, px, pt, ph

OP = 0
OT = 1
OH = 4

# ============================================================
# Setup figure
# ============================================================
fig, ax = plt.subplots(figsize=(12, 7))
ax.set_title("p-h Diagram", fontsize=14)
ax.set_xlabel("h, kJ/kg")
ax.set_ylabel("p, MPa")
ax.set_yscale('log')
ax.grid(True, alpha=0.3)

# ============================================================
# Constants
# ============================================================
Pcrit = 22.064   # MPa
Tcrit = 647.096  # K
tc = Tcrit - 273.15  # °C
Pmin = 611.657e-6  # MPa

# ============================================================
# 1. Saturated liquid line (x=0) - blue
# ============================================================
T_sat = np.linspace(0.01, tc, 300)
h_liq = np.array([tx(t, 0.0, OH) for t in T_sat])
p_liq = np.array([tx(t, 0.0, OP) for t in T_sat])
ax.plot(h_liq, p_liq, 'b-', lw=2.0)

# ============================================================
# 2. Saturated vapor line (x=1) - red
# ============================================================
h_vap = np.array([tx(t, 1.0, OH) for t in T_sat])
p_vap = np.array([tx(t, 1.0, OP) for t in T_sat])
ax.plot(h_vap, p_vap, 'r-', lw=2.0)

# ============================================================
# 3. Isoquality lines inside the dome - green dashed
# ============================================================
T_dome = np.linspace(0.01, tc, 200)
for x in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    h_x = np.array([tx(t, x, OH) for t in T_dome])
    p_x = np.array([tx(t, x, OP) for t in T_dome])
    ax.plot(h_x, p_x, 'g--', lw=0.5, alpha=0.7)

# ============================================================
# 4. Isotherms - black
# ============================================================
# Temperature list matching the reference diagram
T_list = [1, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300,
          325, 350, 375, 400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800]

# Pressure range for isotherms
p_iso = np.logspace(np.log10(Pmin), np.log10(100.0), 200)

for t_c in T_list:
    h_vals = np.array([pt(p, t_c, OH) for p in p_iso])
    ax.plot(h_vals, p_iso, 'k-', lw=0.4, alpha=0.5)
    # Add temperature label at top of each isotherm
    h_top = h_vals[-1]
    p_top = p_iso[-1]
    ax.annotate(f'{t_c}°C', xy=(h_top, p_top),
                xytext=(h_top, p_top * 1.15),
                fontsize=6, ha='center', color='black', rotation=90)

# ============================================================
# Region 2/3 boundary: B23 curve
# ============================================================
# IAPWS-IF97 Eq.(6): T = n4 + sqrt((p - n5) / n3)
# Range: 623.15K (16.5292 MPa) to 863.15K (100 MPa)
# Coefficients from IAPWS-IF97 Table 1
n3_b23 = 0.10192970039326e-2   # 0.0010192970039326
n4_b23 = 0.57254459862746e3    # 572.54459862746
n5_b23 = 0.13918839778870e2    # 13.918839778870
Ps_623 = 16.5291642526045  # MPa, saturation pressure at 623.15K

# Start from saturated vapor point, then B23 curve to 100 MPa
T_start = 623.15  # K
h_start = tx(T_start - 273.15, 1.0, OH)

p_b23 = np.linspace(Ps_623 + 0.01, 100.0, 150)
T_b23 = np.array([n4_b23 + ((p - n5_b23) / n3_b23) ** 0.5 for p in p_b23])
h_b23 = np.array([pt(p, T - 273.15, OH) for p, T in zip(p_b23, T_b23)])

# Prepend start point from saturated vapor line
h_b23 = np.concatenate(([h_start], h_b23))
p_b23 = np.concatenate(([Ps_623], p_b23))

ax.plot(h_b23, p_b23, 'm-', lw=1.5, label='Region 2/3 (B23)')

# ============================================================
# Critical point
# ============================================================
h_crit = pt(Pcrit, tc, OH)
ax.plot(h_crit, Pcrit, 'ko', markersize=5)

# ============================================================
# Display
# ============================================================
ax.set_xlim(0, 4200)
ax.set_ylim(1e-3, 100)

ax.legend(loc='upper right', fontsize=9)
plt.tight_layout()
plt.show()
