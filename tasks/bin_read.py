import sys,os,re,struct
import numpy as np
import matplotlib.pyplot as plt

mu0 = 4*np.pi*1e-7
epsl0 = 8.85*1e-12
dx = 0.02
dy = 0.02
dz = 0.02

nx = 128
ny=128
nz=64
nrec=25
nt=2006
xstep = 2
tstep = 5

xlist1 = os.popen('ls '+os.path.join('pstdtest_800MHz_2.0m_0.5m_fdtd_16','RTM','Output', 'slx_Ey_'+str(0).zfill(4)+'*.bin')).readlines()
print(xlist1[1].strip())
# name = os.path.join('output',xlist1[-1].strip())
name = xlist1[-1].strip()
with open(name,'rb') as fo:
    data_raw = struct.unpack('f'*nz*ny,fo.read())
data = np.reshape(data_raw,(nz,ny))
np.savetxt('data.txt',data)
plt.imshow(data)
plt.colorbar()
plt.show()
plt.savefig('1.png')