#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from pyfem2 import *

V = Plane2DModel()
V.GenesisMesh('../meshes/PlateWithHoleQuad4.g')

V.Material('Material-1')
V.materials['Material-1'].Elastic(E=100, Nu=.2)
V.AssignProperties('ElementBlock1', PlaneStrainQuad4, 'Material-1', t=1)
V.PrescribedBC('LeftHandSide', X, 0.)
V.FixNodes('PinNode')
V.SurfaceLoad(IHI, [1,0])
V.Solve()
V.WriteResults('PlateWithHoleQuad4.exo')

F = File('PlateQuad.exo')
max_p = [0., None]
max_u = [0., None]
for step in F.steps:
    u = step.field_outputs['U']
    for value in u.values:
        u1 = value.magnitude
        if max_u[0] < u1:
            max_u = [u1, value]

    s = step.field_outputs['S']
    for value in s.values:
        s1 = value.max_principal
        if max_p[0] < s1:
            max_p = [s1, value]

# External and internal element numbers
xel = max_p[1].label
x = F.get_elem_coord(xel)
