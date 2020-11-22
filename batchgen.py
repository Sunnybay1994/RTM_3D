#!/usr/bin/env python
import os,sys,getopt
import numpy as np
from model_em import cleanfiles

def subg(dirname,nsrc,job_cap=8,proc_num=12):
    list_src = range(0,nsrc)
    cwd = os.getcwd()
    exepath = cwd
    workdir = os.path.join('tasks',dirname)
    workpath = os.path.join(cwd,workdir)
    stdpath = os.path.join(workpath,'STD')
    rtmpath = os.path.join(workpath,'RTM')
    rtm0path = os.path.join(workpath,'RTM0')

    fname_sub = os.path.join(workpath,'log','sub.sh')
    fname_rtm0 = os.path.join(workpath,'sub_0offset.sh')
    
    with open(fname_sub, 'w') as fp:
            fp.write('''#!/bin/sh
var1="sub_"; var3=".sh"; for((i=0;i<''' + str(job_cap) + ''';i++)); do var2=`echo $i | awk '{printf("%04d\\n",$0)}'`; sbatch ${var1}${var2}${var3}; done''')

    if is_zRTM:
        isrc = 0
        if forward_method == 'fdtd':
            pstd_tag = ''
            execmd_rtm = 'mpiexec -np $NSLOTS -wdir $WORKPATHRTM $EXEPATH/FDTD_MPI.exe '+ str(isrc) +' > $WORKPATHRTM/Output/' + str(isrc) + '.out'
        elif forward_method == 'pstd':
            pstd_tag = ' --pstd'
            execmd_rtm ='''cd $WORKPATHRTM
./PSTD.exe $NSLOTS '''+ str(isrc) +'''
cd $WORKPATH'''

        with open(fname_rtm0, 'w') as fp:
            fp.write('''#!/bin/bash

NSLOTS=''' + str(proc_num) + '''
echo "Got $NSLOTS slots."
echo "PATH = $PATH"
echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
EXEPATH="''' + exepath + '''"
DIRNAME="''' + dirname + '''"
WORKPATH="''' + workpath + '''"
WORKPATHSTD="''' + stdpath + '''"
WORKPATHRTM="''' + rtm0path + '''"

cd $WORKPATH
echo "Current Directory = $WORKPATH"
echo 

#########  execute PROGRAM_NAME
echo "Computing is started at $(date)."

cd $WORKPATH
python $EXEPATH/pre_RTM_sub.py -m ''' + pstd_tag + ''' 0
echo "Current Directory = $WORKPATHRTM"
''' + execmd_rtm +'''

echo "Computing is stopped at $(date)."

exit 0
''')

    if server_name == 'local':
        post_tag = ""
        mpicmd = "mpiexec -np $NSLOTS "
    elif server_name == 'freeosc':
        mpicmd = "$MPI_HOME/bin/mpiexec -n $NSLOTS -iface ib0 -machinefile slurm.hosts "
        post_tag = '#'

    for isrc in list_src:
        if forward_method == 'fdtd':
            pstd_tag = ''
            execmd = mpicmd + '-wdir $WORKPATH $EXEPATH/FDTD_MPI.exe '+ str(isrc) +' > $WORKPATH/Output/' + str(isrc) + '.out'
            execmd_std = mpicmd + '-wdir $WORKPATHSTD $EXEPATH/FDTD_MPI.exe '+ str(isrc) +' > $WORKPATHSTD/Output/' + str(isrc) + '.out'
            execmd_rtm = mpicmd + '-wdir $WORKPATHRTM $EXEPATH/FDTD_MPI.exe '+ str(isrc) +' > $WORKPATHRTM/Output/' + str(isrc) + '.out'
        elif forward_method == 'pstd':
            pstd_tag = ' --pstd'
            execmd = '''cd $WORKPATH
./PSTD.exe $NSLOTS '''+ str(isrc) +' > '+ str(proc_num) +'_threads.out'
            execmd_std = '''cd $WORKPATHSTD
./PSTD.exe $NSLOTS '''+ str(isrc) +'''
cd $WORKPATH'''
            execmd_rtm ='''cd $WORKPATHRTM
./PSTD.exe $NSLOTS '''+ str(isrc) +'''
cd $WORKPATH'''

        fname = os.path.join(workpath,'log','sub_' + str(isrc).zfill(4) + '.sh')
        fname_next = 'sub_' + str(isrc + job_cap).zfill(4) + '.sh'
        fname_post = 'sub_' + str(isrc).zfill(4) + '_post.sh'

        if server_name == 'freeosc':
            ntasks_per_node = 24
            text_sub_next = "sbatch " + fname_next
            text_head = '''#!/bin/sh

# Lines begin with "#SBATCH" set slurm parameters.
# Lines begin with "#" except "#!" and "#SBATCH" are comments.
# Slurm parameters must appear before shell command. 

# Useage: sbatch intel.sh
# Output: slurm-<JOB_ID>.out

#SBATCH --get-user-env
#SBATCH --mail-type=end

######### set job's name
#SBATCH -J ''' + os.path.basename(workpath) + '(' + str(isrc) + ''')
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

######### set NODE and TASK values(CORES = nodes * ntasks-per-node)
#SBATCH --nodes=''' + str(int(np.ceil(proc_num/ntasks_per_node))) + '''
#SBATCH --ntasks-per-node=''' + str(ntasks_per_node) + '''
NSLOTS=''' + str(proc_num) + '''
echo "Got $NSLOTS slots."

######### set Parallel Environment
## load environment before submitting this job
##     module load intel/2019.1.144
export FI_PROVIDER=sockets
export I_MPI_FABRICS=ofi
export FI_SOCKETS_IFACE=ib0

echo "PATH = $PATH"
echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
EXEPATH="''' + exepath + '''"
DIRNAME="''' + dirname + '''"
WORKPATH="''' + workpath + '''"
WORKPATHSTD="''' + stdpath + '''"
WORKPATHRTM="''' + rtmpath + '''"
WORKPATHRTM0="''' + rtm0path + '''"
echo "Current Directory = $WORKPATH"

#env|sort|grep "SLURM"

#########  execute PROGRAM_NAME
echo  "Computing is started at $(date)."

srun hostname | sort -n > slurm.hosts
#sed -i -e 's|compute|fast|g' -e 's|.local||g' slurm.hosts

# -n CORES
'''
        elif server_name == 'local':
            text_sub_next = '''nohup sh "log/''' + fname_next + '''" > "log/''' + str(isrc + job_cap) + '''.out" 2>&1 &'''
            text_head = '''#!/bin/bash

NSLOTS=''' + str(proc_num) + '''
echo "Got $NSLOTS slots."
echo "PATH = $PATH"
echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
EXEPATH="''' + exepath + '''"
DIRNAME="''' + dirname + '''"
WORKPATH="''' + workpath + '''"
WORKPATHSTD="''' + stdpath + '''"
WORKPATHRTM="''' + rtmpath + '''"
WORKPATHRTM0="''' + rtm0path + '''"

cd $WORKPATH
echo "Current Directory = $WORKPATH"
echo 

#########  execute PROGRAM_NAME
echo "Computing is started at $(date)."
'''
        if is_zRTM==1:
            text = '''
''' + execmd +'''
''' + execmd_cmd + '''

''' + text_sub_next + '''

python $EXEPATH/clean.py -f '''+ str(isrc) +'''
'''
        else:
            text = '''
''' + execmd +'''
''' + execmd_std + '''
python $EXEPATH/pre_RTM_sub.py ''' + pstd_tag + ' ' + str(isrc) + '''
''' + execmd_rtm +'''

''' + text_sub_next + '''

''' + post_tag + '''python $EXEPATH/corr_RTM_wavefield_sub.py '''+ str(isrc) +'''
''' + post_tag + '''python $EXEPATH/corr_RTM_slice_sub.py '''+ str(isrc) +'''
''' + post_tag + '''python $EXEPATH/clean.py '''+ str(isrc) +'''
'''

        text_tail = '''
echo "Computing is stopped at $(date)."
'''

        if is_zRTM==1 or is_zRTM==2:
            if isrc == list_src[-1]:
                text_tail += '\ncd $WORKPATH\nsh sub_0offset.sh\n'

        if server_name != 'local':
            if is_zRTM == 1:
                text_tail += '''
python $EXEPATH/post_put.py -z -t ''' + dirname + ' ' + str(isrc) +'''
    '''
            else:
                text_tail += '''
python $EXEPATH/post_put.py -t ''' + dirname + ' ' + str(isrc) +'''
    '''

        text_tail += '\nexit $exit_code\n'

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
# EXEPATH="''' + exepath + '''"
# DIRNAME="''' + dirname + '''"
# WORKPATH="''' + workpath + '''"
# echo "Current Directory = $WORKPATH"
# echo 

# #########  execute PROGRAM_NAME
# echo "Computing is started at $(date)."

# cd $WORKPATH
# python $EXEPATH/corr_RTM_wavefield_sub.py '''+ str(isrc) +'''
# python $EXEPATH/corr_RTM_slice_sub.py '''+ str(isrc) +'''
# python $EXEPATH/clean.py '''+ str(isrc) +'''

# echo "Computing is stopped at $(date)."


# exit 0
# ''')
###########################################################


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "z:s:d:c:p:", ["zero-offset=","src_num=","workdir=","job_capacity=","proc_num=","pstd","server="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        # usage()
        sys.exit(2)

    is_zRTM = 0
    dirname = 'default'
    job_capacity = 8
    proc_num = 8
    forward_method = 'fdtd'
    server_name = 'local'
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
        elif o in ('--pstd',):
            forward_method = 'pstd'
        elif o in ('--server',):
            server_name = a
            print('servername:"' + server_name + '"')
        else:
            assert False, "unhandled option"

    logdir = os.path.join('tasks',dirname,'log')
    if not os.path.exists(logdir):
        os.mkdir(logdir)
    else:
        cleanfiles(logdir)
    print(forward_method)
    subg(dirname,src_num,job_capacity,proc_num)