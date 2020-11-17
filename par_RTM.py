#!/usr/bin/env python
import os
def read_slice(fname):
    with open(fname) as fslice:
        slice_nx,slice_ny,slice_nz = fslice.readline().split(',')
        slice_x = fslice.readline().split(',')
        slice_y = fslice.readline().split(',')
        slice_z = fslice.readline().split(',')
        slice_nx = int(slice_nx);slice_ny = int(slice_ny);slice_nz = int(slice_nz)
        return slice_nx,slice_ny,slice_nz

def read_par(workdir='.'):
    global nx,ny,nz,slice_nx,slice_ny,slice_nz,nt,dx,dy,dz,dt,step_t_wavefield,step_x_wavefield
    with open(os.path.join(workdir,'Input','par.in')) as fpar:
        fpar.readline()
        dx,dy,dz,dt = fpar.readline().split(',')
        print('dx dy dz dt: ',dx,dy,dz,dt) 
        fpar.readline()
        nx,ny,nz,nt = fpar.readline().split(',')
        nx = int(nx);ny = int(ny);nz = int(nz);nt=int(nt)
        print('nx ny nz nt: ',nx,ny,nz,nt) 
        fpar.readline()
        nt_src = fpar.readline()
        print('nt of src: ',nt_src) 
        fpar.readline()
        step_t_wavefield,step_x_wavefield = fpar.readline().split(',')
        step_t_wavefield = int(step_t_wavefield)
        step_x_wavefield = int(step_x_wavefield)
        print('output time step and space step of wavefidld: ',step_t_wavefield,step_x_wavefield) 
        fpar.readline()
        step_slice = fpar.readline()
        print('output step of slice: ',step_slice) 
        fpar.readline()
        npml_x,npml_y,npml_z= fpar.readline().split(',')
        print('npml x y z: ',npml_x,npml_y,npml_z) 
        fpar.readline()
        fpar.readline() #pml m kapxmax kapymax kapzmax alpha
        fpar.readline()
        fsrc= fpar.readline().strip('\n')
        print('src.in: ',fsrc) 
        fpar.readline()
        frec= fpar.readline().strip('\n')
        print('rec.in: ',frec) 
        fpar.readline()
        feps = fpar.readline().strip('\n')
        fpar.readline()
        fmu = fpar.readline().strip('\n')
        fpar.readline()
        fsig= fpar.readline().strip('\n')
        fpar.readline()
        fslice= fpar.readline().strip('\n')
        slice_nx,slice_ny,slice_nz = read_slice(os.path.join(workdir,'Input',fslice))

print('current path: '+os.getcwd())
try:
    read_par()
except Exception as e:
    # raise
    print(e)
else:
    pass
finally:
    pass
#if __name__ == "__main__":
