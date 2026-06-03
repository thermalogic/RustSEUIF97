"""
The Python example to call the shared library

Author:   Cheng Maohua
Email:    cmh@seu.edu.cn

"""
from platform import system
from ctypes import CFUNCTYPE, cdll, c_double, c_int

prototype = CFUNCTYPE(c_double, c_double, c_double, c_int)
prototype_c = CFUNCTYPE(c_double, c_double, c_double)
cdll_names = {'Linux': '../../target/release/libseuif97.so',
              'Windows': '../../target/release/seuif97.dll',
              'Darwin': '../../target/release/libseuif97.dylib'}

osplat = system()
if osplat in cdll_names:
    flib = cdll.LoadLibrary(cdll_names[osplat])
else:
    raise OSError(f"Unsupported platform: {osplat}")


def pt(p, t, pid):
    f = prototype(("pt", flib),)
    result = f(p, t, pid)
    return result

def pt2s(p, t):
    f = prototype_c(("pt2s", flib),)
    result = f(p, t)
    return result

p = 16
t = 535.1
h = pt(p, t, 4)
s = pt2s(p, t)

print("h=", h)  
print("s=", s)
