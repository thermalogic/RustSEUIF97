"""
The simplest example of SEUIF97 Python binding.

Works with both local build and PyPI install:

Local build:
  1. cargo build -r --features python
  2. rename the shared library (*.dll / *.so) to seuif97.pyd
  3. add the path of seuif97.pyd to sys.path before importing

PyPI install:
  pip install seuif97
"""
import sys
import os

# For local build: add the directory containing seuif97.pyd
# Comment out the following line if using the PyPI package
sys.path.insert(0, os.path.abspath(r'D:\sim_rankine_workbench\iapws-if97\RustSEUIF97\target\release'))

from seuif97 import pt, pt2s

OH = 4
OS = 5

p = 16.0
t = 535.1

# universal functions with o_id parameter
h = pt(p, t, OH)
# direct property functions
s = pt2s(p,t)
print(f"p={p}, t={t} h={h:.3f} s={s:.3f}")

h_2000_1 =pt(0.000611212677444, 2000.0, OH)
h_2000_2 =pt(50.0, 2000.0, OH)
print(f"h_2000_1={h_2000_1:.3f} h_2000_2={h_2000_2:.3f}")

s_2000 =pt(0.000611212677444, 2000.0, OS)
print(f"s_2000={s_2000:.16f}")
