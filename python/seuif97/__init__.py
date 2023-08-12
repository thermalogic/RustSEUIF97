from .seuif97 import *

from seuif97 import pt, ph, ps, hs, pv, tv, th, ts, px, tx, hx, sx

OP = 0
OT = 1
OV = 3
OH = 4
OS = 5
OX = 15


def pt2h(p, t):
    return pt(p, t, OH)


def pt2s(p, t):
    return pt(p, t, OS)


def pt2v(p, t):
    return pt(p, t, OV)


def pt2x(p, t):
    return pt(p, t, OX)

# p,h


def ph2t(p, h):
    return ph(p, h, OT)


def ph2s(p, h):
    return ph(p, h, OS)


def ph2v(p, h):
    return ph(p, h, OV)


def ph2x(p, h):
    return ph(p, h, OX)

# p,s


def ps2h(p, s):
    return ps(p, s, OP)


def ps2t(p, s):
    return ps(p, s, OT)


def ps2v(p, s):
    return ps(p, s, OV)


def ps2x(p, s):
    return ps(p, s, OX)

# h,s


def hs2p(h, s):
    return hs(h, s, OP)


def hs2t(h, s):
    return hs(h, s, OT)


def hs2v(h, s):
    return hs(h, s, OV)


def hs2x(h, s):
    return hs(h, s, OX)

# the extended input pairs
# p,v


def pv2t(p, v):
    return pv(p, v, OT)


def pv2h(p, v):
    return pv(p, v, OH)


def pv2s(p, v):
    return pv(p, v, OS)


def pv2x(p, v):
    return pv(p, v, OX)

# t,v


def tv2p(t, v):
    return tv(t, v, OP)


def tv2h(t, v):
    return tv(t, v, OH)


def tv2s(t, v):
    return tv(t, v, OS)


def tv2x(t, v):
    return tv(t, v, OX)

# t,h


def th2p(t, h):
    return th(t, h, OP)


def th2s(t, h):
    return th(t, h, OS)


def th2v(t, h):
    return th(t, h, OV)


def th2x(t, h):
    return th(t, h, OX)

# t,s


def ts2p(t, s):
    return ts(t, s, OP)


def ts2h(t, s):
    return ts(t, s, OH)


def ts2v(t, s):
    return ts(t, s, OV)


def ts2x(t, s):
    return ts(t, s, OX)

# p,x


def px2t(p, x):
    return px(p, x, OT)


def px2h(p, x):
    return px(p, x, OH)


def px2s(p, x):
    return px(p, x, OS)


def px2v(p, x):
    return px(p, x, OV)

# t,x


def tx2p(t, x):
    return tx(t, x, OP)


def tx2h(t, x):
    return tx(t, x, OH)


def tx2s(t, x):
    return tx(t, x, OS)


def tx2v(t, x):
    return tx(t, x, OV)

# h,x


def hx2p(h, x):
    return hx(h, x, OP)


def hx2t(h, x):
    return hx(h, x, OT)


def hx2s(h, x):
    return hx(h, x, OS)


def hx2v(h, x):
    return hx(h, x, OV)
# s,x


def sx2p(s, x):
    return sx(s, x, OP)


def sx2t(s, x):
    return sx(s, x, OT)


def sx2h(s, x):
    return sx(s, x, OH)


def sx2v(s, x):
    return sx(s, x, OV)
