#!/usr/bin/env python
import os,sys,getopt
from model_em import cleanfiles

def subg(dirname,nsrc,job_cap=8,proc_num=8):
    list_src = range(0,nsrc)
    cwd = os.getcwd()
    pypath = cwd
    workdir = os.path.join('tasks',dirname)
    workpath = os.path.join(cwd,workdir)
    stdpath = os.path.join(workpath,'STD')
    rtmpath = os.path.join(workpath,'RTM')
    rtm0path = os.path.join(workpath,'RTM0')

    fname_sub = os.path.join(workpath,'log','sub.sh')
    fname_rtm0 = os.path.join(workpath,'sub_0offset.sh')

    with open(fname_sub, 'w') as fp:
            fp.write('''#!/bin/sh
var1="sub_"; var3=".sh"; for((i=0;i<''' + str(job_cap) + ''';i++)); do var2=`echo $i | awk '{printf("%04d\\n",$0)}'`; qsub ${var1}${var2}${var3}; done''')

    if is_zRTM:
        with open(fname_rtm0, 'w') as fp:
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
PYPATH="''' + pypath + '''"
DIRNAME="''' + dirname + '''"
WORKPATH="''' + workpath + '''"
WORKPATHSTD="''' + stdpath + '''"
WORKPATHRTM="''' + rtm0path + '''"
echo "Current Directory = $WORKPATH"
echo 

#########  execute PROGRAM_NAME
echo "Computing is started at $(date)."

python $PYPATH/pre_RTM_sub.py -m 0 
echo "Current Directory = $WORKPATHRTM"
~/software/openmpi-4.0.3/bin/mpiexec -np $NSLOTS -wdir $WORKPATHRTM $WORKPATHRTM/FDTD_MPI 0

echo "Computing is stopped at $(date)."

exit 0
''')


    for isrc in list_src:
        fname = os.path.join(workpath,'log','sub_' + str(isrc).zfill(4) + '.sh')
        fname_next = 'sub_' + str(isrc + job_cap).zfill(4) + '.sh'
        fname_post = 'sub_' + str(isrc).zfill(4) + '_post.sh'

        text_head = '''#!/bin/sh

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
PYPATH="''' + pypath + '''"
DIRNAME="''' + dirname + '''"
WORKPATH="''' + workpath + '''"
WORKPATHSTD="''' + stdpath + '''"
WORKPATHRTM="''' + rtmpath + '''"
WORKPATHRTM0="''' + rtm0path + '''"
echo "Current Directory = $WORKPATH"
echo 

#########  execute PROGRAM_NAME
echo "Computing is started at $(date)."
'''
        if is_zRTM==1:
            text = '''
~/software/openmpi-4.0.3/bin/mpiexec -np $NSLOTS -wdir $WORKPATH $WORKPATH/FDTD_MPI '''+ str(isrc) +'''
echo "Current Directory = $WORKPATHSTD"
~/software/openmpi-4.0.3/bin/mpiexec -np $NSLOTS -wdir $WORKPATHSTD $WORKPATHSTD/FDTD_MPI '''+ str(isrc) +'''
echo "Current Directory = $WORKPATH"
#~/software/openmpi-4.0.3/bin/mpiexec -np 1 -wdir $WORKPATH $PYPATH/clean.py -f '''+ str(isrc) +'''
'''
        else:
            text = '''
~/software/openmpi-4.0.3/bin/mpiexec -np $NSLOTS -wdir $WORKPATH $WORKPATH/FDTD_MPI '''+ str(isrc) +'''
echo "Current Directory = $WORKPATHSTD"
~/software/openmpi-4.0.3/bin/mpiexec -np $NSLOTS -wdir $WORKPATHSTD $WORKPATHSTD/FDTD_MPI '''+ str(isrc) +'''
echo "Current Directory = $WORKPATH"
~/software/openmpi-4.0.3/bin/mpiexec -np 1 -wdir $WORKPATH $PYPATH/pre_RTM_sub.py '''+ str(isrc) +'''
echo "Current Directory = $WORKPATHRTM"
~/software/openmpi-4.0.3/bin/mpiexec -np $NSLOTS -wdir $WORKPATHRTM $WORKPATHRTM/FDTD_MPI '''+ str(isrc) +'''
#echo "Current Directory = $WORKPATH"
#~/software/openmpi-4.0.3/bin/mpiexec -np 1 -wdir $WORKPATH $PYPATH/corr_RTM_wavefield_sub.py '''+ str(isrc) +'''
#~/software/openmpi-4.0.3/bin/mpiexec -np 1 -wdir $WORKPATH $PYPATH/corr_RTM_slice_sub.py '''+ str(isrc) +'''
#~/software/openmpi-4.0.3/bin/mpiexec -np 1 -wdir $WORKPATH $PYPATH/clean.py '''+ str(isrc) +'''
'''

        text_tail = '''
echo "Computing is stopped at $(date)."

qsub ''' + fname_next + '''
'''

        if is_zRTM==1 or is_zRTM==2:
            if isrc == list_src[-1]:
                text_tail += '\ncd ..\nsh sub_0offset.sh\n'

        if is_zRTM == 1:
            text_tail += '''
cd $PYPATH
python post_put.py -z -t ''' + dirname + ' ' + str(isrc) +'''
'''
        else:
            text_tail += '''
cd $PYPATH
python post_put.py -t ''' + dirname + ' ' + str(isrc) +'''
'''
        text_tail += '\nexit 0\n'

        with open(fname, 'w') as fp:
            fp.write(text_head+text+text_tail)

########################### rtm_post ######################
#
#  depicted because of high I/O cost.
#  replaced by post_master.py & post_put.py
#
###########################################################


#         if not is_zRTM:
#             with open(os.path.join(workpath,'log',fname_post), 'w') as fp:
#                 fp.write('''#!/bin/sh

# # Lines begin with "#$" are parameters of qsub.
# # Lines begin with "#" except "#!" and "#$" are comments.

# # Useage: qsub openmpi.sh
# # Output: <JOB_NAME>.o<JOB_ID>

# #$ -cwd
# #$ -m beas
# #$ -j y
# #$ -S /bin/sh
# #$ -v PATH,LD_LIBRARY_PATH

# ######### set the name of this job
# #$ -N rtm'''+ str(isrc) +'''_post

# echo "PATH = $PATH"
# echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
# PYPATH="''' + pypath + '''"
# DIRNAME="''' + dirname + '''"
# WORKPATH="''' + workpath + '''"
# echo "Current Directory = $WORKPATH"
# echo 

# #########  execute PROGRAM_NAME
# echo "Computing is started at $(date)."

# cd $WORKPATH
# python $PYPATH/corr_RTM_wavefield_sub.py '''+ str(isrc) +'''
# python $PYPATH/corr_RTM_slice_sub.py '''+ str(isrc) +'''
# python $PYPATH/clean.py '''+ str(isrc) +'''

# echo "Computing is stopped at $(date)."


# exit 0
# ''')
###########################################################


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "z:s:d:c:p", ["zero-offset=","src_num=","workdir=","job_capacity=","proc_num="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        # usage()
        sys.exit(2)

    is_zRTM = 0
    dirname = 'default'
    job_capacity = 8
    proc_num = 8
    for o, a in opts:
        if o in ('-z','--zero-offset'):
            is_zRTM = int(a)
            print('Zero-offset Mode: %d'%is_zRTM)
        elif o in ('-s','--src_num'):
            src_num = int(a)
        elif o in ('-d','--workdir'):
            dirname = a
        elif o in ('-c','--job_capacity'):
            job_capacity = int(a)
        elif o in ('-p','--proc_num'):
            proc_num = int(a)
        else:
            assert False, "unhandled option"

    logdir = os.path.join('tasks',dirname,'log')
    if not os.path.exists(logdir):
        os.mkdir(logdir)
    else:
        cleanfiles(logdir)

    subg(dirname,src_num,job_capacity,proc_num)