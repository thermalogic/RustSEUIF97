"""
The Python example to call the shared library

Author:   Cheng Maohua
Email:    cmh@seu.edu.cn

"""
from platform import *
from ctypes import *

prototype = CFUNCTYPE(c_double, c_double, c_double, c_int)
cdll_names = {'Linux': '../target/release/libif97.so',
              'Windows': '../target/release/if97.dll'}

osplat = system()
if (osplat == 'Linux'):
    flib = cdll.LoadLibrary(cdll_names[osplat])
elif (osplat == 'Windows'):
    flib = cdll.LoadLibrary(cdll_names[osplat])

def pt(p, t, pid):
    f = prototype(("pt", flib),)
    result = f(p, t, pid)
    return result

p=16
t=535.1
h=pt(p,t,4)
print("h=",h)