#!/usr/bin/env python
import os,sys,argparse
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

    if mode == 'z':
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
python $EXEPATH/pre_RTM_sub.py -m z''' + pstd_tag + '''
echo "Current Directory = $WORKPATHRTM"
''' + execmd_rtm +'''

echo "Computing is stopped at $(date)."

exit 0
''')

    for isrc in list_src:
        if isrc >= nsrc:
            next_task_tag="# "
        else:
            next_task_tag=""
        if server_name == 'local':
            post_tag = ""
            mpicmd = "mpiexec -np $NSLOTS "
        elif server_name == 'freeosc':
            hostfile = "slurm" + str(isrc).zfill(4) + ".hosts"
            mpicmd = "$MPI_HOME/bin/mpiexec -n $NSLOTS -iface ib0 -machinefile " + hostfile + " "
            post_tag = '' if 'i' in steps else '# '
        if forward_method == 'fdtd':
            pstd_tag = ''
            exit_txt = 'exit_code=$?;echo $exit_code;'
            execmd = ('' if 'g' in steps else '# ') + 'echo "($(date))Generate data...";' + mpicmd + '-wdir $WORKPATH $EXEPATH/FDTD_MPI.exe '+ str(isrc) +' > $WORKPATH/Output/' + str(isrc) + '.out;' + exit_txt
            execmd_std = ('' if 'f' in steps else '# ') + 'echo "($(date))Calculate forward wavafield...";' + mpicmd + '-wdir $WORKPATHSTD $EXEPATH/FDTD_MPI.exe '+ str(isrc) +' > $WORKPATHSTD/Output/' + str(isrc) + '.out;' + exit_txt
            execmd_rtm = ('' if 'b' in steps else '# ') + 'echo "($(date))Calculate backward wavafield...";' + mpicmd + '-wdir $WORKPATHRTM $EXEPATH/FDTD_MPI.exe '+ str(isrc) +' > $WORKPATHRTM/Output/' + str(isrc) + '.out;' + exit_txt
            execmd_rtm0 = ('' if 'b' in steps else '# ') + 'echo "($(date))Calculate backward wavafield...";' + mpicmd + '-wdir $WORKPATHRTM0 $EXEPATH/FDTD_MPI.exe '+ str(isrc) +' > $WORKPATHRTM0/Output/' + str(isrc) + '.out;' + exit_txt
        elif forward_method == 'pstd':
            pstd_tag = ' --pstd'
            execmd = '''cd $WORKPATH;echo $PWD
./PSTD.exe $NSLOTS '''+ str(isrc) +' > '+ str(proc_num) +'_threads.out'
            execmd_std = '''cd $WORKPATHSTD;echo $PWD
./PSTD.exe $NSLOTS '''+ str(isrc) +'''
cd $WORKPATH'''
            execmd_rtm ='''cd $WORKPATHRTM;echo $PWD
./PSTD.exe $NSLOTS '''+ str(isrc) +'''
cd $WORKPATH'''

        fname = os.path.join(workpath,'log','sub_' + str(isrc).zfill(4) + '.sh')
        fname_next = 'sub_' + str(isrc + job_cap).zfill(4) + '.sh'
        fname_post = 'sub_' + str(isrc).zfill(4) + '_post.sh'

        if server_name == 'freeosc':
            ntasks_per_node = 24
            if ntasks_per_node > proc_num:
                ntasks_per_node = proc_num
                nodes = 1
            else:
                nodes = int(np.ceil(proc_num/ntasks_per_node))
            text_sub_next = "cd log;sbatch " + fname_next + ";cd .."
            text_head = '''#!/bin/sh

# Lines begin with "#SBATCH" set slurm parameters.
# Lines begin with "#" except "#!" and "#SBATCH" are comments.
# Slurm parameters must appear before shell command. 

# Useage: sbatch intel.sh
# Output: slurm-<JOB_ID>.out

#SBATCH --get-user-env
#SBATCH --mail-type=end

######### set job's name
#SBATCH -J (''' + str(isrc) + ')' + os.path.basename(workpath) + '''
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

######### set NODE and TASK values(CORES = nodes * ntasks-per-node)
#SBATCH --nodes=''' + str(nodes) + '''
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

cd $WORKPATH
echo "Current Directory = $WORKPATH"
echo "Current Source NO.: ''' + str(isrc) + '''"

#env|sort|grep "SLURM"

# sleep ''' + str(isrc*5) + '''
#########  execute PROGRAM_NAME
echo  "Computing is started at $(date)."

srun hostname | sort -n > ''' + hostfile + '''
#sed -i -e 's|compute|fast|g' -e 's|.local||g' ''' + hostfile + '''

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
        if mode == 'z':
            text = '''
''' + execmd +'''
''' + execmd_std + '''
# echo "($(date))Prepare for RTM..."
# python $EXEPATH/pre_RTM_sub.py ''' + pstd_tag + ' ' + str(isrc) + ''' -m z
# ''' + execmd_rtm0 +'''

''' + next_task_tag + text_sub_next + '''

python $EXEPATH/clean.py -f '''+ str(isrc) +'''
'''
        else:
            text = '''
''' + execmd +'''
''' + execmd_std + '''
''' + ('' if 'b' in steps else '# ') + '''echo "($(date))Prepare for RTM..."
''' + ('' if 'b' in steps else '# ') + '''python $EXEPATH/pre_RTM_sub.py ''' + pstd_tag + ' ' + str(isrc) + '''
''' + execmd_rtm +'''

echo "($(date))Submitting next task..."
''' + next_task_tag + text_sub_next + '''

''' + post_tag + '''echo "($(date))Applying image condition..."
''' + post_tag + '''python $EXEPATH/corr_RTM_wavefield_sub.py '''+ str(isrc) +''' &
''' + post_tag + '''python $EXEPATH/corr_RTM_slice_sub.py '''+ str(isrc) +''' &
''' + post_tag + '''wait
''' + post_tag + '''echo "($(date))Cleaning..."
''' + post_tag + '''python $EXEPATH/clean.py '''+ str(isrc) +'''
'''

        text_tail = '''
echo "Computing is stopped at $(date)."
'''

        if 'z' in mode:
            if isrc == list_src[-1]:
                text_tail += '\ncd $WORKPATH\nsh sub_0offset.sh\n'

        if server_name == 'x3850':
            if mode == 'z':
                text_tail += '''
python $EXEPATH/post_put.py -z -t ''' + dirname + ' ' + str(isrc) +'''
    '''
            else:
                text_tail += '''
python $EXEPATH/post_put.py -t ''' + dirname + ' ' + str(isrc) +'''
    '''
        if server_name == 'freeosc':
            text_tail += '''
/bin/rm -f ''' + hostfile + '''
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


#         if mode == 'm':
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
    parser = argparse.ArgumentParser(description='Generate batch jobs for server to run.',conflict_handler='resolve')
    parser.add_argument('-d','--workdir',default='default',dest='dirname',help="work directory of the task.")
    parser.add_argument('-s','--src_num',type=int,required=True,help="The total number of shot gathers (sources).")
    parser.add_argument('-m','--mode',choices=['m','z','mz'],default='m',help="Mode: 'm' for multi-offset, 'z' for zero-offset, 'mz' for both multi- and zero-offset.")
    parser.add_argument('--forward_method',choices=['fdtd','pstd'],default='fdtd',help="Forward method used in RTM.")
    parser.add_argument('--fdtd',action='store_const',const='fdtd',dest='forward_method',help='Use finite difference time domain as the forward method.')
    parser.add_argument('--pstd',action='store_const',const='pstd',dest='forward_method',help='Use pseudo spectral time domain as the forward method.')
    parser.add_argument('--server',choices=['local','freeosc','x3850'],default='local',help="Where to run the code.")
    parser.add_argument('-p','--np',type=int,default=-1,help='Number of processers/threads used in parallel FDTD/PSTD. Default: 6/12 for local, 8/16 for x3850 and freeosc.')
    parser.add_argument('-c','--job_cap',type=int,default=-1,help='How may shot gathers (sources) should the server handle at the same time. Default: 1 for local, 8 for x3850 and 32 for freeosc.')
    parser.add_argument('-y',action='store_const',const='y',dest='noprompt',default=False,help="Input 'y' in all input prompts with no disturbing.")
    parser.add_argument('-n',action='store_const',const='n',dest='noprompt',default=False,help="Input 'n' in all input prompts with no disturbing.")
    parser.add_argument('--steps',type=str,default='gfbi',help="Select which steps are involved. 'g' for generate data; 'f' for source(forward) wavefield; 'b' for receiver(backward) wavefield; 'i' for imaging.")
    args = parser.parse_args()
    
    forward_method = args.forward_method
    server = args.server
    noprompt = args.noprompt
    steps = args.steps
    if server == 'local':
        pnum = 6
        job_cap = 1
    else:
        pnum = 8
        if server == 'x3850':
            job_cap = 8
        elif server == 'freeosc':
            job_cap = 32
    if forward_method == 'pstd':
        pnum *= 2
    if args.np > 0:
        pnum = args.np
    if args.job_cap > 0:
        job_cap = args.job_cap
    print('pnum=%d, job_capacity=%d'%(pnum,job_cap))

    mode = args.mode
    server_name = args.server
    print('servername:"%s", mode:"%s"'%(server_name,mode))

    logdir = os.path.join('tasks',args.dirname,'log')
    if not os.path.exists(logdir):
        os.mkdir(logdir)
    else:
        cleanfiles(logdir,noprompt)
    print('forward_method:%s'%forward_method)
    if job_cap > args.src_num:
        job_cap = args.src_num
    subg(args.dirname,args.src_num,job_cap,pnum)