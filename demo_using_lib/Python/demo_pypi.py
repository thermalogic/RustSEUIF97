"""
The Python example to call seuif97

Author:   Cheng Maohua
Email:    cmh@seu.edu.cn

"""
from seuif97 import pt, pt2s

OH = 4
OS = 5

p = 16.0
t = 535.1

# Universal Functions (with o_id parameter))
h = pt(p, t, OH)
# Direct Property Functions
s = pt2s(p,t)
print(f"p={p}, t={t} h={h:.3f} s={s:.3f}")
