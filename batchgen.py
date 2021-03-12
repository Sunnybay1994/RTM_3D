#!/usr/bin/env python
import os,sys,argparse
import numpy as np
from model_em import cleanfiles

class rtm_workflow:
    def __init__(self,taskname='default',nsrc=1,job_cap=1,proc_num=4,method='fdtd',steps='gfbizc',mpipath='mpiexec',subcmd='nohup bash {script_name} {logcmd} 2>&1 &',script_name='script{s_isrc}.sh'):
        self._taskname = taskname
        self._nsrc = nsrc
        self._job_cap = job_cap
        self._proc_num = proc_num
        self._method = method
        self._mpipath = mpipath
        self._subcmd = subcmd
        self._steps = steps
        self._name_script = script_name
        self._rootpath = os.getcwd()
        self._exepath = self._rootpath
        self._workpath = os.path.join(self._rootpath,'tasks',self._taskname)

    @property
    def rootpath(self):
        return self._rootpath
    @rootpath.setter
    def rootpath(self,path):
        self._rootpath = path
    @property
    def exepath(self):
        return self._exepath
    @exepath.setter
    def exepath(self,path):
        self._exepath = path
    @property
    def workpath(self):
        return self._workpath
    @workpath.setter
    def workpath(self,path):
        self._workpath = path
    @property
    def stdpath(self):
        return os.path.join(self.workpath,'STD')
    @property
    def rtmpath(self):
        return os.path.join(self.workpath,'RTM')
    @property
    def rtm0path(self):
        return os.path.join(self.workpath,'RTM0')
    @property
    def logpath(self):
        return os.path.join(self.workpath,'log')
    @property
    def name_script(self):
        return self._name_script
    @name_script.setter
    def name_script(self,prefix):
        self._name_script = '%s{s_isrc}.sh'%prefix
    def get_script_name(self,isrc):
        if isinstance(isrc,int):
            s_isrc = s_isrc='_%04d'%isrc
        else:
            s_isrc = isrc
        return self.name_script.format(s_isrc=s_isrc)
    @property
    def name_script_sub(self):
        return 'sub_{0}'.format(self.get_script_name(''))

    @property
    def txt_head(self):
        txt = '''#!/bin/bash
{additional_head_txt}
######### parameters
ISRC={isrc}
STEPS=''' + self._steps + '''
METHOD_TAG=''' + ('--%s'%self._method if self._method != 'fdtd' else '') + '''
NSRC=''' + str(self._nsrc) + '''
echo "ISRC=$ISRC, STEPS=$STEPS, method_tag=$METHOD_TAG, nsrc=$NSRC"

######### environment
NP=''' + str(self._proc_num) + '''
echo "Using $NP processes."
echo "PATH = $PATH"
echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
EXEPATH="''' + self.exepath + '''"
DIRNAME="''' + self._taskname + '''"
WORKPATH="''' + self.workpath + '''"
WORKPATHSTD="''' + self.stdpath + '''"
WORKPATHRTM="''' + self.rtmpath + '''"
WORKPATHRTM0="''' + self.rtm0path + '''"
LOGPATH="''' + self.logpath + '''"
cd $WORKPATH
echo "Current Directory = $WORKPATH"
echo "Current Source NO.: $ISRC"

#########  execute PROGRAM_NAME
echo  "Computing is started at $(date)."
'''
        return txt
    def gen_txt_head(self,isrc,txt=''):
        return self.txt_head.format(isrc=isrc,additional_head_txt=txt)

    def gen_mpicmd_txt(self,wdir):
        return '{mpipath} -np $NP -wdir {wdir} $EXEPATH/FDTD_MPI.exe $ISRC > "{wdir}/Output/$ISRC($NP).out"\n'.format(mpipath=self._mpipath,wdir=wdir)
    def gen_forward_txt(self,step):
        head_txt = 'if [[ $STEPS =~ %s ]];then\n'%step
        exit_txt = '    exit_code=$?;echo exit_code=$exit_code\nfi\n'
        if step == 'g':
            head_txt += '    echo "($(date))Generate data..."\n'
            path = '$WORKPATH'
        elif step == 'f':
            head_txt += '    echo "($(date))Calculating forward wavafield..."\n'
            path = '$WORKPATHSTD'
        elif step == 'b':
            head_txt += '''    echo "($(date))Prepare for RTM..."
    python $EXEPATH/pre_RTM_sub.py $METHOD_TAG $ISRC
    echo "($(date))Calculating backward wavafield..."
'''
            path = '$WORKPATHRTM'
        elif step == 'z':
            head_txt += '''    echo "($(date))Prepare for zero-offset RTM..."
    python $EXEPATH/pre_RTM_sub.py -m z $METHOD_TAG $ISRC
    status=%?;echo status=$status
    if [ $status -eq 1 ];then
    echo "($(date))Calculating backward wavafield..."
'''
            exit_txt += 'fi\n'
            path = '$WORKPATHRTM0'
        else:
            raise ValueError('rtm_workflow.gen_forward_cmd_txt(step=%s):forward step must be in "gfbz"'%step)

        if self._method == 'fdtd':
            main_txt = '    ' + self.gen_mpicmd_txt(path)
        elif self._method == 'pstd':
            main_txt = "    cd %s;./PSTD.exe $NP $ISRC;cd $WORKPATH;\n"%path
        return head_txt + main_txt + exit_txt
    def gen_txt_main(self):
        return ''.join([self.gen_forward_txt(step) for step in 'gfbz'])

    def gen_txt_tail(self,txt=''):
        if not txt:
            txt = '''
if [[ $STEPS =~ r ]];then
    if [[ $STEPS =~ c ]];then
        cleanflag=""
    else
        cleanflag="--no_clean"
    fi
    if [[ $STEPS =~ i ]];then
        python $EXEPATH/post_put.py -t $DIRNAME $cleanflag
    elif [[ $STEPS =~ z ]];then
        python $EXEPATH/post_put.py -z -t $DIRNAME $cleanflag
    fi
else
    zoflag=""
    if [[ $STEPS =~ i ]];then
        echo "($(date))Applying image condition..."
        python $EXEPATH/corr_RTM_wavefield_sub.py $ISRC &
        python $EXEPATH/corr_RTM_slice_sub.py $ISRC &
        wait
    elif [[ $STEPS =~ z ]];then
        zoflag="-f"
    fi
    if [[ $STEPS =~ c ]];then
        echo "($(date))Cleaning..."
        python $EXEPATH/clean.py $zoflag $ISRC
    fi
fi
echo "Computing is stopped at $(date)."

exit $exit_code
'''
        return txt

    @property
    def sub_next_txt(self):
        return "echo submitting {script_name_next};\ncd $LOGPATH;{subcmd}\ncd $WORKPATH"
    def gen_txt_next(self,isrc):
        isrc_next = isrc + self._job_cap
        if isrc_next < self._nsrc:
            script_name_next = self.get_script_name(isrc_next)
            logcmd = '> {isrc_next}.out'.format(isrc_next=isrc_next)
            return self.sub_next_txt.format(script_name_next=script_name_next,subcmd=self._subcmd.format(script_name=script_name_next,logcmd=logcmd))
        else:
            return ''

    def gen_script(self,isrc):
        content = self.gen_txt_head(isrc)+self.gen_txt_main()+self.gen_txt_next(isrc)+self.gen_txt_tail()
        with open(os.path.join(self.logpath,self.get_script_name(isrc)),'w+') as fo:
            fo.write(content)

    def gen_script_submit(self):
        var = '''do var=`echo $i | awk '{printf("%04d\\n",$0)}'`'''
        isrc = '${var}'
        script_name=self.get_script_name('_'+isrc)
        logcmd = '> {isrc}.out'.format(isrc=isrc)
        subcmd = self._subcmd.format(script_name=script_name,logcmd=logcmd)
        content = '#!/bin/bash\nfor((i=0;i<{job_cap};i++)); {var}; `{subcmd}`; done'.format(job_cap=self._job_cap,var=var,subcmd=subcmd)
        with open(os.path.join(self.logpath,self.name_script_sub),'w+') as fo:
            fo.write(content)

    def gen(self):
        self.gen_script_submit()
        for i in range(self._nsrc):
            self.gen_script(i)





def subg(dirname,nsrc,job_cap=8,proc_num=12):
    list_src = range(0,nsrc)
    cwd = os.getcwd()
    exepath = cwd
    workdir = os.path.join('tasks',dirname)
    workpath = os.path.join(cwd,workdir)

    fname_sub = os.path.join(workpath,'log','sub.sh')
    fname_rtm0 = os.path.join(workpath,'sub_0offset.sh')

    with open(fname_sub, 'w') as fp:
            fp.write('''#!/bin/sh
var1="sub_"; var3=".sh"; for((i=0;i<''' + str(self._job_cap) + ''';i++)); do var2=`echo $i | awk '{printf("%04d\\n",$0)}'`; sbatch ${var1}${var2}${var3}; done''')

    def subtxt(isrc,np,dirname,exepath,workpath,steps):
        stdpath = os.path.join(workpath,'STD')
        rtmpath = os.path.join(workpath,'RTM')
        rtm0path = os.path.join(workpath,'RTM0')

        def gen_wp_txt(isrc,workdir,np,method,cmd='mpiexec',info="'Calculating wavefield propagating...'"):
            if method == 'pstd':
                return "echo %s;cd %s;echo $PWD;./PSTD.exe %d %d  >  %s/Output/%d(%d).out;cd $WORKPATH;exit_code=$?;echo exit_code=$exit_code;"%(info,workdir,np,isrc,wrkdir,isrc,np)
            elif method == 'fdtd':
                return "echo %s;%s -n %d -wdir %s $EXEPATH/FDTD_MPI.exe %d > %s/Output/%d(%d).out;exit_code=$?;echo exit_code=$exit_code;"%(info,cmd,np,workdir,isrc,workdir,isrc,np)
            else:
                txt = "echo %s;%s;exit_code=$?;echo exit_code=$exit_code;"%(info,cmd)
                # print("self-defined command: %s"%txt)
                return txt

        txt = '''
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

#########  execute PROGRAM_NAME
echo  "Computing is started at $(date)."

''' + gen_wp_txt(isrc,'$WORKPATH',np,method,cmd,info='($(date))Generate data...') +'''
''' + gen_wp_txt(isrc,'$WORKPATHSTD',np,method,cmd,info='($(date))Generate data...') + '''
''' + ('' if 'b' in steps else '# ') + '''echo "($(date))Prepare for RTM...;python $EXEPATH/pre_RTM_sub.py ''' + pstd_tag + ' %d'%isrc + '''
''' + gen_wp_txt(isrc,'$WORKPATHRTM',np,method,cmd,info='($(date))Calculate backward wavafield...') + '''

echo "($(date))Submitting next task..."
''' + next_task_tag + text_sub_next + '''

''' + post_tag + '''echo "($(date))Applying image condition..."
''' + post_tag + '''python $EXEPATH/corr_RTM_wavefield_sub.py '''+ str(isrc) +''' &
''' + post_tag + '''python $EXEPATH/corr_RTM_slice_sub.py '''+ str(isrc) +''' &
''' + post_tag + '''wait
''' + post_tag + '''echo "($(date))Cleaning..."
''' + post_tag + '''python $EXEPATH/clean.py '''+ str(isrc) +'''
'''



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
        if isrc + job_cap >= nsrc:
            next_task_tag="# "
        else:
            next_task_tag=""
        if server_name == 'local':
            post_tag = ""
            mpicmd = "mpiexec -np $NSLOTS "
        elif server_name == 'freeosc':
            # hostfile = "slurm" + str(isrc).zfill(4) + ".hosts"
            # mpicmd = "$MPI_HOME/bin/mpiexec -n $NSLOTS -iface ib0 -machinefile " + hostfile + " "
            mpicmd = "$MPI_HOME/bin/mpiexec "
            post_tag = '' if 'i' in steps else '# '
        if forward_method == 'fdtd':
            pstd_tag = ''
            exit_txt = 'exit_code=$?;echo exit_code=$exit_code;'
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
            cpus_per_task = proc_num if forward_method=='pstd' else 1
            ntasks_per_node = proc_num if forward_method=='fdtd' else 1
            nodes = 1
            text_sub_next = "cd log;sbatch " + fname_next + ";cd .."
            text_head = '''#!/bin/sh

# Lines begin with "#SBATCH" set sbatch parameters.
# Lines begin with "#" except "#!" and "#SBATCH" are comments.
# Sbatch parameters must appear before shell command.

# Useage: sbatch intel.sh
# Output: slurm-<JOB_ID>.out

#SBATCH --get-user-env
#SBATCH --mail-type=end
##SBATCH --mail-user=username@pku.edu.cn

######### set job's name
#SBATCH --job-name=(''' + str(isrc) + ')' + os.path.basename(workpath) + '''
##SBATCH --output=slurm-%j.out
##SBATCH --error=slurm-%j.err

######### set NODE and TASK values(CORES = nodes * ntasks-per-node * cpus-per-task)
#SBATCH --nodes=''' + str(nodes) + '''
#SBATCH --ntasks-per-node=''' + str(ntasks_per_node) + '''
#SBATCH --cpus-per-task=''' + str(cpus_per_task) + '''

######### set Parallel Environment
## load environment before submitting this job
##     module load intel/2020.4.304

echo "JOB_NODELIST: ${SLURM_JOB_NODELIST}"
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

#########  execute PROGRAM_NAME
echo  "Computing is started at $(date)."

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
        if 'z' in mode:
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
    parser.add_argument('-d','--workdir',default='default',dest='dirname',help="work directory (basename) of the task.")
    parser.add_argument('-s','--src_num',type=int,required=True,help="The total number of shot gathers (sources).")
    parser.add_argument('-m','--mode',choices=['m','z','mz'],default='m',help="Mode: 'm' for multi-offset, 'z' for zero-offset, 'mz' for both multi- and zero-offset.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--forward_method',choices=['fdtd','pstd'],default='fdtd',help="Forward method used in RTM.")
    group.add_argument('--fdtd',action='store_const',const='fdtd',dest='forward_method',help='Use finite difference time domain as the forward method.')
    group.add_argument('--pstd',action='store_const',const='pstd',dest='forward_method',help='Use pseudo spectral time domain as the forward method.')
    parser.add_argument('--server',choices=['local','freeosc','x3850'],default='local',help="Where to run the code.")
    parser.add_argument('-p','--np',type=int,default=-1,help='Number of processers/threads used in parallel FDTD/PSTD. Default: 6/12 for local, 8/16 for x3850 and freeosc.')
    parser.add_argument('-c','--job_cap',type=int,default=-1,help='How may shot gathers (sources) should the server handle at the same time. Default: 1 for local, 8 for x3850 and 32 for freeosc.')
    parser.add_argument('-y',action='store_const',const='y',dest='noprompt',default=False,help="Input 'y' in all input prompts with no disturbing.")
    parser.add_argument('-n',action='store_const',const='n',dest='noprompt',default=False,help="Input 'n' in all input prompts with no disturbing.")
    parser.add_argument('--steps',type=str,default='',help="Manually Select which steps are involved. 'g' for generate data; 'f' for source(forward) wavefield; 'b' for receiver(backward) wavefield; 'i' for cross-corrlation image condition; 'z' for cal zero-offset backward wvf(image condition); 'c' for clean middle result, add 'r' option if you want to run the 'i' and 'c' process on other servers for they do not involve parallel calculation(set up 'post_master/worker.py' first). Default: Use 'gfbic' for a whole multi-offset workflow, 'gfzc' for a whole zero-offset workflow, 'gfbizc' for both.")
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
    if args.steps:
        steps = args.steps
    else:
        if mode == 'z':
            steps = 'gfzc'
        elif mode == 'm':
            steps = 'gfbic'
        elif mode == 'mz':
            steps = 'gfbizc'

    logdir = os.path.join('tasks',args.dirname,'log')
    if not os.path.exists(logdir):
        os.mkdir(logdir)
    else:
        cleanfiles(logdir,noprompt)
    print('forward_method:%s'%forward_method)
    if job_cap > args.src_num:
        job_cap = args.src_num

    workflow = rtm_workflow(args.dirname,args.src_num,job_cap,pnum,forward_method,steps)
    workflow.gen()
    # subg(args.dirname,args.src_num,job_cap,pnum)