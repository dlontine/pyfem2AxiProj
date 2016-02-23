#!/usr/bin/env python
import sys
sys.path.insert(0, '../')
from pyfem2 import *

V = Plane2DModel()
V.GenesisMesh('../meshes/PlateWithHoleQuad4QuarterSym.g')
V.Material('Material-1')
V.materials['Material-1'].Elastic(E=100, Nu=.2)
V.AssignProperties('', PlaneStrainQuad4, 'Material-1', t=1)
V.PrescribedBC('SymYZ', X, 0.)
V.PrescribedBC('SymXZ', Y, 0.)
V.SurfaceLoad('RightHandSide', [1,0])
V.Solve()
V.WriteResults('PlateWithHoleQuad4QuarterSym.exo')

F = File('QuarterPlate.exo')
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
print max_p[0]
