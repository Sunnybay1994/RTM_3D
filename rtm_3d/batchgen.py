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
    @property
    def name_script_sub(self):
        return 'sub_{0}'.format(self.get_script_name(''))

    @property
    def txt_head(self):
        txt = '''#!/bin/bash
{additional_head_txt}

######### parameters
ISRC={isrc}
DIRNAME="''' + self._taskname + '''"
STEPS=''' + self._steps + '''
METHOD_TAG=''' + ('--%s'%self._method if self._method != 'fdtd' else '') + '''
NSRC=''' + str(self._nsrc) + '''

######### environment
NP=''' + str(self._proc_num) + '''
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
echo "Current Source NO.: $ISRC/$NSRC"
SISRC=`echo $ISRC | awk '{{printf("%04d",$0)}}'`;
echo isrc_str=$SISRC

#########  execute PROGRAM_NAME
echo  "=====Computing started at $(date)====="
cd $WORKPATH
echo "Current Directory = $WORKPATH"

if [[ -f $WORKPATHRESULT/result_wavefield_corr_${{SISRC}}.bin && -f $WORKPATHRESULT/result_xcorr_${{SISRC}}_00.bin \\
&& -f $WORKPATHRESULT/result_ycorr_${{SISRC}}_00.bin && -f $WORKPATHRESULT/result_zcorr_${{SISRC}}_00.bin ]]; then 
    doneflag=true
else
    doneflag=false
fi
echo "($ISRC)doneflag=$doneflag"


'''
        return txt
    def gen_txt_head(self,isrc,txt=''):
        return self.txt_head.format(isrc=isrc,additional_head_txt=txt.format(isrc=isrc))

    def gen_mpicmd_txt(self,wdir,execmd,isrctag='$ISRC'):
        return '{mpipath} -wdir {wdir} {execmd} {isrctag} > "{wdir}/Output/{isrctag}($NP).out"\n'.format(mpipath=self._mpipath,wdir=wdir,execmd=execmd,isrctag=isrctag)
    def gen_forward_txt(self,step):
        head_txt = 'if [[ $STEPS =~ %s || ( $STEPS =~ %s && $doneflag = "false" ) ]];then\n'%(step.upper(),step)
        exit_txt = '    exit_code=$?;echo exit_code=$exit_code\nfi\n'
        isrctag = '$ISRC'
        if step == 'g':
            head_txt += '    echo "($(date))Generate data..."\n'
            path = '$WORKPATH'
        elif step == 'f':
            head_txt += '    echo "($(date))Calculating forward wavafield..."\n'
            path = '$WORKPATHSTD'
        elif step == 'b':
            head_txt = '''if [[ ${STEPS,,} =~ pb ]];then
    echo "($(date))Prepare for RTM..."
    python $SRCPATH/pre_RTM_sub.py $METHOD_TAG $ISRC
fi
''' + head_txt + '    echo "($(date))Calculating backward wavafield..."\n'
            path = '$WORKPATHRTM'
        elif step == 'z':
            head_txt += '''    if [[ ${STEPS,,} =~ pz ]];then
        echo "($(date))Prepare for zero-offset RTM..."
        python $SRCPATH/pre_RTM_sub.py -m z $METHOD_TAG $ISRC
        status=$?;echo status=$status
    else
        status=100
    fi
    if [ $status -eq 100 ];then
    echo "($(date))Calculating backward wavafield..."
'''
            exit_txt = '    fi\n' + exit_txt
            path = '$WORKPATHRTM0'
            isrctag = '0'
        else:
            raise ValueError('rtm_workflow.gen_forward_cmd_txt(step=%s):forward step must be in "gfbz"'%step)

        if self._method == 'fdtd':
            main_txt = '    ' + self.gen_mpicmd_txt(path,'-np $NP $BINPATH/FDTD_MPI.exe',isrctag)
        elif self._method == 'pstd':
            main_txt = '    ' + self.gen_mpicmd_txt(path,'-np 1 $BINPATH/PSTD.exe $NP',isrctag)
        return head_txt + main_txt + exit_txt
    def gen_txt_main(self):
        return ''.join([self.gen_forward_txt(step) for step in 'gfbz'])
    
    @property
    def txt_tail(self):
        return '''
if [[ $STEPS =~ I || ($STEPS =~ i && $doneflag = "false") ]];then
    echo "($(date))Applying image condition..."
    corrflag=
else
    corrflag='--nocorr'
fi
if [[ ${{STEPS,,}} =~ c ]];then
    echo "($(date))clean files..."
    cleanflag=
else
    cleanflag='--noclean'
fi
python $SRCPATH/apply_image_condition.py $ISRC $cleanflag $corrflag -m ''' + self._mode + '''


echo "=====Computing stopped at $(date)====="
{additional_tail_txt}
exit $exit_code
'''
    def gen_txt_tail(self,txt=''):
        return self.txt_tail.format(additional_tail_txt=txt)

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
        content = self.gen_txt_head(isrc,self._add_head)+self.gen_txt_main()+self.gen_txt_next(isrc)+self.gen_txt_tail(self._add_tail)
        with open(os.path.join(self.logpath,self.get_script_name(isrc)),'w+') as fo:
            fo.write(content)

    def gen_script_submit(self):
        var = '''do var=`echo $i | awk '{printf("%04d\\n",$0)}'`'''
        isrc = '${var}'
        script_name=self.get_script_name('_'+isrc)
        logcmd = '> {isrc}.out'.format(isrc=isrc)
        subcmd = self._subcmd.format(script_name=script_name,logcmd=logcmd)
        content = '#!/bin/bash\nfor((i=0;i<{job_cap};i++)); {var}; `echo "{subcmd}"`; done'.format(job_cap=self._job_cap,var=var,subcmd=subcmd)
        with open(os.path.join(self.logpath,self.name_script_sub),'w+') as fo:
            fo.write(content)

    def gen(self):
        self.gen_script_submit()
        for i in range(self._nsrc):
            self.gen_script(i)

class rtm_workflow_freeosc(rtm_workflow):
    additional_head_txt = '''
# Lines begin with "#SBATCH" set sbatch parameters.
# Lines begin with "#" except "#!" and "#SBATCH" are comments.
# Sbatch parameters must appear before shell command.

# Useage: sbatch intel.sh
# Output: slurm-<JOB_ID>.out

#SBATCH --get-user-env
#SBATCH --mail-type=end
##SBATCH --mail-user=username@pku.edu.cn

######### set job's name
#SBATCH --job-name={method}{isrc}_{taskname}
#SBATCH --output={scriptname}-%j.out
#SBATCH --error={scriptname}-%j.out

######### set NODE and TASK values(CORES = nodes * ntasks-per-node * cpus-per-task)
#SBATCH --ntasks={np}

######### set Parallel Environment
## load environment before submitting this job
##     module load intel/2020.4.304

echo "JOB_NODELIST: ${{SLURM_JOB_NODELIST}}"
'''
    def __init__(self,taskname='default',nsrc=1,job_cap=60,proc_num=4,method='fdtd',mode='z',steps='gfbizc',mpipath='$MPI_HOME/bin/mpiexec',subcmd='sbatch {script_name}',script_name='script{s_isrc}.sh',add_head=additional_head_txt,add_tail=''):
        # command changing exists file: sed -i "s/sbatch script/sbatch --exclude=cu17 script/g" *.sh 
        rtm_workflow.__init__(self,taskname,nsrc,job_cap,proc_num,method,mode,steps,mpipath,subcmd,script_name,add_head,add_tail)
    def gen_txt_head(self,isrc,txt=''):
        return self.txt_head.format(isrc=isrc,additional_head_txt=txt.format(isrc=isrc,method=self._method[0],taskname=self._taskname,np=self._proc_num,scriptname=self._name_script.format(s_isrc="_%04d"%isrc)))

def batchgen(args):
    global logger
    logger=addlogger(os.path.basename(sys.argv[0]),args.dirname,path=logpath)
    logger.info(' '.join(sys.argv))

    if args.steps:
        steps = ''.join(args.steps)
    else:
        mode = args.mode
        if mode == 'z':
            steps = 'gfpzc'
        elif mode == 'm':
            steps = 'gfpbic'
        elif mode == 'mz':
            steps = 'gfpbipzc'
    if args.max_cpu:
        job_cap = args.max_cpu//args.np
    else:
        job_cap = args.max_job
    if job_cap > args.src_num:
        job_cap = args.src_num
    logger.info('Job capacity:%d\nforward method: %s\nsteps: %s\n'%(job_cap,args.forward_method,steps))

    if args.server == 'freeosc':
        workflow = rtm_workflow_freeosc(args.dirname,args.src_num,job_cap,args.np,args.forward_method,args.mode,steps)
    else:
        workflow = rtm_workflow(args.dirname,args.src_num,job_cap,args.np,args.forward_method,args.mode,steps)
    workflow.gen()

if __name__ == '__main__':
    parser = parser_ini(2)
    args = parser.parse_args()
    batchgen(args)