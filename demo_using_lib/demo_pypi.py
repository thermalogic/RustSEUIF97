"""
The Python example to call seuif97

Author:   Cheng Maohua
Email:    cmh@seu.edu.cn

"""
from seuif97 import pt, pt2h,pt2s

OH = 4
OS = 5

p = 16.0
t = 535.1

# ??(in1,in2,o_id)
h = pt(p, t, OH)
s = pt(p, t, OS)
print(f"??(in1,in2,o_id): p={p}, t={t} h={h:.3f} s={s:.3f}")

# ??2?(in1,in2)
h = pt2h(p,t)
s = pt2s(p, t)
print(f"   ??2?(in1,in2): p={p}, t={t} h={h:.3f} s={s:.3f}")
