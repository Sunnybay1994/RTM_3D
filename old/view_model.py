#!/usr/bin/env python
from numpy import *
from matplotlib.pyplot import *
from par_RTM import *
from mayavi import mlab

num_processor = 8
order = 2
eps = zeros([nx, ny, nz])

ix = 0
for i in range(num_processor):
    print('./Input/eps.in_' + str(i).zfill(3)) 
    dum_eps = loadtxt('./Input/eps.in_' + str(i).zfill(3))
    # modified by mbw at 20190417
    dum_nx = int(len(dum_eps) / ny)
    # dum_nx = len(dum_eps) / ny
    dum_nx2 = dum_nx - order * 2
    # dum_nx2 = dum_nx - 4
    print(dum_nx2)

    dum_eps = reshape(dum_eps, (dum_nx, ny, nz))
    eps[ix:ix + dum_nx2, :, :] = dum_eps[order:-order, :, :]
    ix += dum_nx2

field = mlab.pipeline.scalar_field(eps)
#mlab.pipeline.volume(field,vmin=0, vmax=8)
cut = mlab.pipeline.scalar_cut_plane(field, plane_orientation="x_axes")
#cut.implicit_plane.widget.enabled = False
cut = mlab.pipeline.scalar_cut_plane(field, plane_orientation="y_axes")
#cut.implicit_plane.widget.enabled = False
cut = mlab.pipeline.scalar_cut_plane(field, plane_orientation="z_axes")
#cut.implicit_plane.origin = (30-2, 0, 0)
#cut.implicit_plane.widget.enabled = False
mlab.outline()