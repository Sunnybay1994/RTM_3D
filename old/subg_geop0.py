#!/usr/bin/env python
import os

nsrc = 121
job_cap = 8

proc_num = 8
list_src = range(0,nsrc)
cwd = os.getcwd()
stdpath = os.path.join(cwd,'STD')
rtmpath = os.path.join(cwd,'RTM')

fname_sub = 'log/sub.sh'
fname_rtm = 'sub_0offset.sh'

with open(fname_sub, 'w') as fp:
        fp.write('''#!/bin/sh
var1="sub_"; var3=".sh"; for((i=0;i<''' + str(job_cap) + ''';i++)); do var2=`echo $i | awk '{printf("%04d\\n",$0)}'`; qsub ${var1}${var2}${var3}; done''')

with open(fname_rtm, 'w') as fp:
    fp.write('''#!/bin/sh

# Lines begin with "#$" are parameters of qsub.
# Lines begin with "#" except "#!" and "#$" are comments.

# Useage: qsub openmpi.sh
# Output: <JOB_NAME>.o<JOB_ID>

#$ -cwd
#$ -m beas
#$ -j y
#$ -S /bin/sh
#$ -v PATH,LD_LIBRARY_PATH

######### set the name of this job
#$ -N rtm_0offset

######### set Parallel Environment and CORE numbers
module unload mpi/mpich-x86_64
module load mpi/openmpi-x86_64
#$ -pe openmpi 8

echo "Got $NSLOTS slots."
echo "PATH = $PATH"
echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
WORKPATH="''' + cwd + '''"
WORKPATHSTD="''' + stdpath + '''"
WORKPATHRTM="''' + rtmpath + '''"
echo "Current Directory = $WORKPATH"
echo 

#########  execute PROGRAM_NAME
echo "Computing is started at $(date)."

python pre_RTM_sub.py -m 0 
echo "Current Directory = $WORKPATHRTM"
~/software/openmpi-4.0.1/bin/mpiexec -np $NSLOTS -wdir $WORKPATHRTM $WORKPATHRTM/FDTD_MPI_geop 0

echo "Computing is stopped at $(date)."

exit 0
''')


for isrc in list_src:
    fname = 'log/sub_' + str(isrc).zfill(4) + '.sh'
    fname_next = 'sub_' + str(isrc + job_cap).zfill(4) + '.sh'

    text = '''#!/bin/sh

# Lines begin with "#$" are parameters of qsub.
# Lines begin with "#" except "#!" and "#$" are comments.

# Useage: qsub openmpi.sh
# Output: <JOB_NAME>.o<JOB_ID>

#$ -cwd
#$ -m beas
#$ -j y
#$ -S /bin/sh
#$ -v PATH,LD_LIBRARY_PATH

######### set the name of this job
#$ -N rtm'''+ str(isrc) +'''

######### set Parallel Environment and CORE numbers
module unload mpi/mpich-x86_64
module load mpi/openmpi-x86_64
#$ -pe openmpi '''+ str(proc_num) +'''

echo "Got $NSLOTS slots."
echo "PATH = $PATH"
echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
WORKPATH="''' + cwd + '''"
WORKPATHSTD="''' + stdpath + '''"
WORKPATHRTM="''' + rtmpath + '''"
echo "Current Directory = $WORKPATH"
echo 

#########  execute PROGRAM_NAME
echo "Computing is started at $(date)."

~/software/openmpi-4.0.1/bin/mpiexec -np $NSLOTS -wdir $WORKPATH $WORKPATH/FDTD_MPI_geop '''+ str(isrc) +'''
echo "Current Directory = $WORKPATHSTD"
~/software/openmpi-4.0.1/bin/mpiexec -np $NSLOTS -wdir $WORKPATHSTD $WORKPATHSTD/FDTD_MPI_geop '''+ str(isrc) +'''
echo "Current Directory = $WORKPATH"
~/software/openmpi-4.0.1/bin/mpiexec -np 1 -wdir $WORKPATH $WORKPATH/clean.py -f '''+ str(isrc) +'''
# ~/software/openmpi-4.0.1/bin/mpiexec -np 1 -wdir $WORKPATH $WORKPATH/pre_RTM_sub.py -m 0 '''+ str(isrc) +'''

echo "Computing is stopped at $(date)."
qsub ''' + fname_next + '''
'''
    if isrc == list_src[-1]:
        text += 'cd ..\nsh sub_0offset.sh\n'
    text += '\nexit 0\n'

    with open(fname, 'w') as fp:
        fp.write(text)
