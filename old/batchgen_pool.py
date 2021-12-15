#!/usr/bin/env python
import os,sys,argparse
import numpy as np
from common import *

class rtm_workflow(object):
    def __init__(self,taskname='default',nsrc=1,job_cap=2,proc_num=2,method='fdtd',mode='z',steps='gfbizc',mpipath='mpiexec',subcmd='nohup bash {script_name} {logcmd} 2>&1 &',script_name='script{s_isrc}.sh',add_head='',add_tail=''):
        self._taskname = taskname
        self._nsrc = nsrc
        self._job_cap = job_cap
        self._proc_num = proc_num
        self._method = method
        self._mode = mode
        self._mpipath = mpipath
        self._subcmd = subcmd
        self._steps = steps
        self._name_script = script_name
        self._add_head = add_head
        self._add_tail = add_tail

        self._rootpath = rootdir
        self._srcpath = srcpath
        self._binpath = binpath
        self._workpath = os.path.join(taskpath,self._taskname)

    @property
    def rootpath(self):
        return self._rootpath
    @rootpath.setter
    def rootpath(self,path):
        self._rootpath = path
    @property
    def srcpath(self):
        return self._srcpath
    @property
    def binpath(self):
        return self._binpath
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
    def resultpath(self):
        return os.path.join(self.workpath,'Result')
    @property
    def logpath(self):
        path = os.path.join(self.workpath,'log')
        if not os.path.isdir(path):
            os.mkdir(path)
        return path
    @property
    def name_script(self):
        return self._name_script
    @name_script.setter
    def name_script(self,prefix):
        self._name_script = '%s{s_isrc}.sh'%prefix
    def get_script_name(self,isrc):
        if isinstance(isrc,int):
            s_isrc = '_%04d'%isrc
        else:
            s_isrc = isrc
        return self.name_script.format(s_isrc=s_isrc)
    def name_script_sub(self,pref=''):
        return 'sub{pref}_{name}'.format(pref=pref,name=self.get_script_name(''))

    @property
    def txt_head(self):
        txt = '''#!/bin/bash
{additional_head_txt}

######### parameters
DIRNAME="''' + self._taskname + '''"
STEPS=''' + self._steps + '''
METHOD_TAG=''' + ('--%s'%self._method if self._method != 'fdtd' else '') + '''
NSRC=''' + str(self._nsrc) + '''

######### environment
NP=''' + str(self._proc_num) + '''
NTASK=''' + str(self._job_cap) + '''
echo "PATH = $PATH"
echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
SRCPATH="''' + self.srcpath + '''"
BINPATH="''' + self.binpath + '''"
WORKPATH="''' + self.workpath + '''"
WORKPATHSTD="''' + self.stdpath + '''"
WORKPATHRTM="''' + self.rtmpath + '''"
WORKPATHRTM0="''' + self.rtm0path + '''"
WORKPATHRESULT="''' + self.resultpath + '''"
LOGPATH="''' + self.logpath + '''"
echo
echo "NP=$NP, STEPS=$STEPS, method_tag=$METHOD_TAG"
echo

myjob(){{
ISRC=$1
echo "Current Source NO.: $ISRC/$NSRC"
SISRC=`echo $ISRC | awk '{{printf("%04d",$0)}}'`;
echo isrc_str=$SISRC

#########  execute PROGRAM_NAME
echo  "($ISRC)=====Computing started at $(date)====="
cd $WORKPATH
echo "Current Directory = $WORKPATH"

if [[ -f $WORKPATHRESULT/result_wavefield_corr_${{SISRC}}.dat && -f $WORKPATHRESULT/result_xcorr_${{SISRC}}_00.dat \\
&& -f $WORKPATHRESULT/result_ycorr_${{SISRC}}_00.dat && -f $WORKPATHRESULT/result_zcorr_${{SISRC}}_00.dat ]]; then
    doneflag=true
else
    doneflag=false
fi
echo "($ISRC)doneflag=$doneflag"

'''
        return txt
    def gen_txt_head(self,isrc,txt=''):
        return self.txt_head.format(isrc=isrc,s_isrc="%04d"%isrc,additional_head_txt=txt.format(isrc=isrc))

    def gen_mpicmd_txt(self,wdir,execmd,isrctag='$ISRC'):
        return '{mpipath} -wdir {wdir} {execmd} {isrctag} > "{wdir}/Output/{isrctag}($NP).out"\n'.format(mpipath=self._mpipath,wdir=wdir,execmd=execmd,isrctag=isrctag)
    def gen_forward_txt(self,step):
        if step != 'z':
            head_txt = '''if [[ $STEPS =~ ''' + "%s"%step + ''' && $doneflag = "false" ]];then
    exit_code=1
    while (( $exit_code != 0 ));do
'''
            exit_txt = '''    exit_code=$?;echo "($ISRC)exit_code=$exit_code"
    if (( $exit_code != 0 ));then
        sleep 60
    fi
    done
fi\n'''
        else:
            head_txt = '''if [[ $STEPS =~ ''' + "%s"%step + ''' && $doneflag = "false" ]];then\n'''
            exit_txt = '''    exit_code=$?;echo "($ISRC)exit_code=$exit_code"\nfi\n'''
        isrctag = '$ISRC'
        if step == 'g':
            head_txt += '    echo "($ISRC)($(date +%F\ %T))Generate data..."\n'
            path = '$WORKPATH'
        elif step == 'f':
            head_txt += '    echo "($ISRC)($(date +%F\ %T))Calculating forward wavafield..."\n'
            path = '$WORKPATHSTD'
        elif step == 'b':
            head_txt += '''    echo "($ISRC)($(date +%F\ %T))Prepare for RTM..."
    python $SRCPATH/pre_RTM_sub.py $METHOD_TAG $ISRC
    echo "($ISRC)($(date +%F\ %T))Calculating backward wavafield..."
'''
            path = '$WORKPATHRTM'
        elif step == 'z':
            head_txt += '''    echo "($ISRC)($(date +%F\ %T))Prepare for zero-offset RTM..."
    python $SRCPATH/pre_RTM_sub.py -m z $METHOD_TAG $ISRC
    status=$?;echo "($ISRC)status=$status"
    if [ $status -eq 100 ];then
    echo "($ISRC)($(date +%F\ %T))Calculating backward wavafield..."
'''
            exit_txt += 'fi\n'
            path = '$WORKPATHRTM0'
            isrctag = '0'
        else:
            raise ValueError('rtm_workflow.gen_forward_cmd_txt(step=%s):forward step must be in "gfbz"'%step)

        if self._method == 'fdtd':
            main_txt = '    ' + self.gen_mpicmd_txt(path,'-n $NP $BINPATH/FDTD_MPI.exe',isrctag)
        elif self._method == 'pstd':
            main_txt = '    ' + self.gen_mpicmd_txt(path,'-n 1 $BINPATH/PSTD.exe $NP',isrctag)
        return head_txt + main_txt + exit_txt
    def gen_txt_main(self):
        return ''.join([self.gen_forward_txt(step) for step in 'gfbz'])
    
    @property
    def txt_tail(self):
        return '''

if [[ $STEPS =~ i && $doneflag = "false" ]];then
    echo "($ISRC)($(date +%F\ %T))Applying image condition..."
    corrflag=
else
    corrflag='--nocorr'
fi

if [[ $STEPS =~ c ]];then
    echo "($ISRC)($(date +%F\ %T))clean files..."
    cleanflag=
else
    cleanflag='--noclean'
fi
python $SRCPATH/apply_image_condition.py $ISRC $cleanflag $corrflag -m ''' + self._mode + '''


echo "($ISRC)=====Computing stopped at $(date)====="
}}

{additional_tail_txt}

exit $exit_code
'''
    def gen_txt_tail(self,txt=''):
        return self.txt_tail.format(additional_tail_txt=txt)

    @property
    def sub_next_txt(self):
        return "# echo submitting {script_name_next};\n# cd $LOGPATH;{subcmd}\n# cd $WORKPATH"
    def gen_txt_next(self,isrc):
        isrc_next = isrc + self._job_cap
        if isrc_next < self._nsrc:
            script_name_next = self.get_script_name(isrc_next)
            logcmd = '> {isrc_next}.out'.format(isrc_next=isrc_next)
            return self.sub_next_txt.format(script_name_next=script_name_next,subcmd=self._subcmd.format(script_name=script_name_next,logcmd=logcmd))
        else:
            return ''

    def gen_script(self,isrc,mode='standalone'):
        pool_tail = '''##### Process Pool Begin #####
# 设置并发数
TASK_NUM=$NTASK
echo "TASK_NUM=$TASK_NUM"

# 以主进程PID命名管道文件
FIFO_FILE="$$.fifo"
rm -f ${FIFO_FILE}

# 新建管道文件
mkfifo ${FIFO_FILE}
echo pipefile=${FIFO_FILE}

# 生成文件描述符9指向管道文件;"<"表示可读,">"表示可写
exec 9<>${FIFO_FILE}

# 向文件描述符中写入设置的并发数量的行数,模拟生成进程池;一行就是一个子进程
for process_num in $(seq ${TASK_NUM})
do
		  echo "${process_num}" >&9
done

availCores(){
    echo `sinfo -Nt idle,mixed  -o "%14C %4n %9P %10T %8O %8e %E" | grep cpu | sed 's/^\(...........\).*/\\1/' | tr "/" "\\n" | sed -n '2~4p' |  tr "\\n" "+" | sed -e 's/+$//'` | bc
}

# 按并发数设置执行任务
for i in $(seq $NSRC);do
    isrc=$(($NSRC-$i))
#     ac=`availCores`

    # while (( $ac < $NP ));do
    #     sleep 30
    #     ac=`availCores`
    # done

    # SISRC=`echo $isrc | awk '{printf("%04d",$0)}'`
    # echo submitting script_$SISRC.sh
    # sbatch script_$SISRC.sh
    # sleep 10

    # if (( $isrc > $PROC_NUM && $ac > $NP )); then
    #     SISRC=`echo $isrc | awk '{printf("%04d",$0)}'`
    #     echo submitting script_$SISRC.sh
    #     sbatch script_$SISRC.sh
    #     sleep 10
    # else
        #从文件描述符中读取一行,读到后生成子进程向下执行,读不到就等待;模拟进程锁
        read -u 9 P       
        {
            # 子进程执行的实际任务
            echo "(${P})$(date +%F\ %T) src${isrc} Start."
            myjob ${isrc} > ${isrc}.out
            echo "(${P})$(date +%F\ %T) src${isrc} End."
            # 将从文件描述符中读取的内容重新写回文件描述符;模拟释放进程锁
            echo ${P} >&9
        } &
    # fi
done


wait

echo "All Completed"

# 删除文件描述符
exec 9>&-

# 删除临时管道文件
rm -f ${FIFO_FILE}

# ref: https://blog.51cto.com/784687488/2501555
##### Process Pool END #####'''
        standalone_tail = "myjob %d"%isrc
        if mode == 'pool':
            tail = pool_tail
            fn = os.path.join(self.logpath,self.get_script_name('_pool'))
        elif mode == 'standalone':
            tail = standalone_tail
            fn = os.path.join(self.logpath,self.get_script_name(isrc))
        else:
            tail = self._add_tail
        content = self.gen_txt_head(isrc,self._add_head)+self.gen_txt_main()+self.gen_txt_next(isrc)+self.gen_txt_tail(tail)
        with open(fn,'w+') as fo:
            fo.write(content)

    def gen_script_submit(self):
        var = '''do var=`echo $i | awk '{printf("%04d\\n",$0)}'`'''
        isrc = '${var}'
        script_name=self.get_script_name('_'+isrc)
        logcmd = '> {isrc}.out'.format(isrc=isrc)
        subcmd = self._subcmd.format(script_name=script_name,logcmd=logcmd)
        content = '#!/bin/bash\nfor((i=0;i<{job_cap};i++)); {var}; `echo "{subcmd}"`; done'.format(job_cap=self._job_cap,var=var,subcmd=subcmd)
        with open(os.path.join(self.logpath,self.name_script_sub()),'w+') as fo:
            fo.write(content)

    def gen_script_submit_pool(self,pref='_pool'):
        content = '''#!/bin/bash
. ~/.shrc
cp {rawfn} script_pool.sh
NP={np}
# split=$1 # set to an uint if run multiples tasks simultaneously
# echo split=${{split:=1}} # default is 1
echo cores_per_job=$NP
COREAVAIL=`availCores`
echo COREAVAIL=$COREAVAIL
# COREMAX=$(($COREAVAIL/$split/$NP*$NP))
# COREMAX=1 # use as a master
COREMAX=$1
echo COREMAX=${{COREMAX:=24}}
sed -i "s/#SBATCH --cpus-per-task={np}/#SBATCH --cpus-per-task=$COREMAX/" script_pool.sh
sed -i "s/#MAXCORES#/$COREMAX/" script_pool.sh
# sbatch script_pool.sh
'''.format(rawfn=self.get_script_name('_pool_raw'),np=self._proc_num)
        with open(os.path.join(self.logpath,self.name_script_sub(pref)),'w+') as fo:
            fo.write(content)

    def gen(self):
        # self.gen_script_submit()
        # self.gen_script_submit_pool()
        for i in range(self._nsrc):
            self.gen_script(i)
        self.gen_script(-1,mode='pool')

class rtm_workflow_freeosc(rtm_workflow):
    additional_head_txt = '''
# Lines begin with "#SBATCH" set sbatch parameters.
# Lines begin with "#" except "#!" and "#SBATCH" are comments.
# Sbatch parameters must appear before shell command.

# Useage: sbatch intel.sh
# Output: slurm-<JOB_ID>.out

#SBATCH --get-user-env
#SBATCH --mail-type=end
#SBATCH --mail-user=mbw@pku.edu.cn

######### set job's name
#SBATCH --job-name={method}{isrc}_{taskname}
#SBATCH --output={scriptname}-%j.out
#SBATCH --error={scriptname}-%j.out

######### set NODE and TASK values(CORES = nodes * ntasks-per-node * cpus-per-task)
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH --ntasks={jobcap}
#SBATCH --cpus-per-task={np}

######### set Parallel Environment
## load environment before submitting this job
##     module load intel/2020.4.304

echo "JOB_NODELIST: ${{SLURM_JOB_NODELIST}}"
'''
    def __init__(self,taskname='default',nsrc=1,job_cap=60,proc_num=4,method='fdtd',mode='z',steps='gfbizc',mpipath='$MPI_HOME/bin/mpiexec',subcmd='sbatch {script_name}',script_name='script{s_isrc}.sh',add_head=additional_head_txt,add_tail=''):
        rtm_workflow.__init__(self,taskname,nsrc,job_cap,proc_num,method,mode,steps,mpipath,subcmd,script_name,add_head,add_tail)
    def gen_txt_head(self,isrc,txt=''):
        return self.txt_head.format(isrc=isrc,s_isrc="%04d"%isrc,additional_head_txt=txt.format(isrc=isrc,method=self._method[0],taskname=self._taskname,jobcap=self._job_cap if isrc==-1 else 1,np=self._proc_num,scriptname=self._name_script.format(s_isrc="_%04d"%isrc)))

def batchgen(args):
    if args.steps:
        steps = ''.join(args.steps)
    else:
        mode = args.mode
        if mode == 'z':
            steps = 'gfzc'
        elif mode == 'm':
            steps = 'gfbic'
        elif mode == 'mz':
            steps = 'gfbizc'
    if args.max_cpu:
        job_cap = args.max_cpu//args.np
    else:
        job_cap = args.max_job
    if job_cap > args.src_num:
        job_cap = args.src_num
    print('Job capacity:%d\nforward method: %s\nsteps: %s\n'%(job_cap,args.forward_method,steps))

    if args.server == 'freeosc':
        workflow = rtm_workflow_freeosc(args.dirname,args.src_num,job_cap,args.np,args.forward_method,args.mode,steps)
    else:
        workflow = rtm_workflow(args.dirname,args.src_num,job_cap,args.np,args.forward_method,args.mode,steps)
    workflow.gen()

if __name__ == '__main__':
    parser = parser_ini(2)
    args = parser.parse_args()
    batchgen(args)

    bashrc_bkp='''alias mrun="matlab -nodesktop -nosplash -logfile `date +%Y_%m_%d-%H_%M_%S`.log -r"
alias msq='squeue -o "%7i %1P %16j %5u %1T %11M %4l %3C %2D %R"'
msqm(){
    msq | grep $1 
} 
mscancel(){
    temp=`msq | grep $1 | sed 's/^\(.......\).*/\1/' | tr "\n" "," | sed -e 's/,$/\n/'`;
    echo scancel $temp
    scancel $temp
}
alias msi='sinfo -Nt idle,mixed  -o "%14C %4n %9P %10T %8O %8e %E"'
availCores(){
    echo `msi | grep cpu | sed 's/^\(...........\).*/\1/' | tr "/" "\n" | sed -n '2~4p' |  tr "\n" "+" | sed -e 's/+$//'` | bc
}
'''
