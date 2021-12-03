#!/usr/bin/env python
import os,sys,functools
def writeit(out=False):
    def dec(func):
        @functools.wraps(func)
        def wrapper(*args,**kwargs):
            if out==True:
                wrapper_result = func(*args,**kwargs)
            elif isinstance(out,str):
                with open(out,'a+') as fo:
                    wrapper_result = fo.write(*args,**kwargs)
            else:
                wrapper_result = None
            return wrapper_result
        return wrapper
    return dec

@writeit()
def myprint(*args,**kwargs):
    print(*args,**kwargs)


def read_slice(fname):
    lines = [line.strip() for line in open(fname).readlines()]
    lines.reverse()
    slice_nx,slice_ny,slice_nz = [int(n) for n in lines.pop().split(',')]
    slice_x = [int(lines.pop().split(',')[0]) for i in range(slice_nx)]
    slice_y = [int(lines.pop().split(',')[0]) for i in range(slice_ny)]
    slice_z = [int(lines.pop().split(',')[0]) for i in range(slice_nz)] 
    return slice_nx,slice_ny,slice_nz,slice_x,slice_y,slice_z

def read_par(workdir='.'):
    global nx,ny,nz,slice_nx,slice_ny,slice_nz,slice_x,slice_y,slice_z,nt,dx,dy,dz,dt,step_t_wavefield,step_x_wavefield
    with open(os.path.join(workdir,'Input','par.in')) as fpar:
        myprint('***** Reading Parameters *****')
        fpar.readline()
        dx,dy,dz,dt = [float(x) for x in fpar.readline().split(',')]
        myprint('dx=%gm, dy=%gm, dz=%gm, dt=%gns'%(dx,dy,dz,dt/1e9)) 
        fpar.readline()
        nx,ny,nz,nt = fpar.readline().split(',')
        nx = int(nx);ny = int(ny);nz = int(nz);nt=int(nt)
        myprint('nx=%d, ny=%d, nz=%d, nt=%d'%(nx,ny,nz,nt))
        fpar.readline()
        nt_src = int(fpar.readline())
        myprint('nt_src=%d'%nt_src) 
        fpar.readline()
        step_t_wavefield,step_x_wavefield = fpar.readline().split(',')
        step_t_wavefield = int(step_t_wavefield)
        step_x_wavefield = int(step_x_wavefield)
        myprint('output time and space step of wavefield: %d %d'%(step_t_wavefield,step_x_wavefield)) 
        fpar.readline()
        step_slice = int(fpar.readline())
        myprint('output time step of slice: %d'%step_slice) 
        fpar.readline()
        npml_x,npml_y,npml_z= [int(x) for x in fpar.readline().split(',')]
        myprint('npmlx=%d, npmly=%d, npmlz=%d'%(npml_x,npml_y,npml_z))
        fpar.readline()
        fpar.readline() #pml m kapxmax kapymax kapzmax alpha
        fpar.readline()
        fsrc= fpar.readline().strip('\n')
        myprint('source info: %s'%fsrc) 
        fpar.readline()
        frec= fpar.readline().strip('\n')
        myprint('receiver info: %s'%frec) 
        fpar.readline()
        feps = fpar.readline().strip('\n')
        fpar.readline()
        fmu = fpar.readline().strip('\n')
        fpar.readline()
        fsig= fpar.readline().strip('\n')
        fpar.readline()
        fslice= fpar.readline().strip('\n')
        slice_nx,slice_ny,slice_nz,slice_x,slice_y,slice_z = read_slice(os.path.join(workdir,'Input',fslice))
        myprint('slice_nx,slice_ny,slice_nz,slice_x,slice_y,slice_z: ',slice_nx,slice_ny,slice_nz,slice_x,slice_y,slice_z)
        myprint('***** End of Parameters Reading *****')

# if __name__ == "__main__":
try:
    print('%s Call "par_RTM.py" current path: %s'%(sys.argv[0],os.getcwd()))
    read_par()
except Exception as e:
    # raise
    print('%s: %s'%(sys.argv[0],e))
else:
    pass
finally:
    pass
