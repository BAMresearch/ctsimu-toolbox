# -*- coding: UTF-8 -*-
from ctsimu.primitives import *
from math import pi

M = Matrix(cols=3, rows=3)
M.make_identity()
M.scale(pi)
print(M)

N = Matrix(cols=3, rows=3)
N.make_identity()
N.scale(4)
#print(N)

O = M + 3
print(O)

U = Vector(x=pi, y=2*pi, z=3*pi)
V = Vector(x=5, y=6, z=7)

#print(U)
#print(V)
#print(U.to(V))