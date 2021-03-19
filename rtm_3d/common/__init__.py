import os,re,argparse,glob,shutil,datetime,logging,struct
from .writesource import *
from .normal_moveout import *
from .par_RTM import *

rootdir = os.path.abspath('..')
srcpath = os.path.join(rootdir,'rtm_3d')
modelpath = os.path.join(srcpath,'make_model')
binpath = os.path.join(rootdir,'bin')
taskpath = os.path.join(rootdir,'tasks')
logpath = os.path.join(rootdir,'log')

def cp(f1,f2):
    for f in glob.glob(r'%s'%f1):
        # logger.debug('cp %s %s'%(f,f2))
        shutil.copy(f,f2)

def addlogger(name,taskname='',path=logpath,streamlevel=logging.INFO,loglevel=logging.INFO):
    today = datetime.date.today()
    fn = os.path.join(path,'{name}_{taskname}_{date}.log'.format(name=name,taskname=taskname,date=today.strftime('%Y%m%d')))
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
    fh = logging.FileHandler(fn)
    fh.setLevel(streamlevel)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('(%(process)5d)%(asctime)s-%(levelname)s: %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.info('log path:%s'%fn)
    return logger

def scandir_re_match(path='.',ftype='file',*patterns):
    # type in ['file','dir','all']
    results={pattern:[] for pattern in patterns}
    with os.scandir(path) as it:
        for entry in it:
            if (ftype == 'file' and entry.is_dir()) or (ftype == 'dir' and entry.is_file()):
                continue
            for pattern in patterns:
                if re.fullmatch(pattern,entry.name):
                    results[pattern].append(entry)
    return (results[pattern] for pattern in patterns)

def filename_re_match(path='.',ftype='file',*patterns):
    p_entrys = scandir_re_match(path,ftype,*patterns)
    return [[entry.name for entry in entrys] for entrys in p_entrys]

def read_bin_data(fn, dims, type='f'):
    with open(fn,'rb') as fo:
        data_raw = struct.unpack(type*dims[0]*dims[1]*dims[2],fo.read())
    dims.reverse()
    data = np.reshape(data_raw,dims)
    data = data.transpose(2,1,0)
    # print(np.shape(data))
    return data

def cleanfiles(paths,noprompt=False):
    def cleanChoice(path,noprompt):
        if noprompt == 'y':
            return 1
        elif noprompt == 'n':
            return 0
        if not paths:
            choice = input('Clear %s? Y/(N)'%path).lower()
            if choice in ['yes','y']:
                return 1
            else:
                return 0
        choice = input('Clear %s? A/Y/(N)'%path).lower()
        if choice in ['all', 'a']:
            confirm = input("Clear all files in '%s'? Y/(N)" % "' & '".join(paths+[path])).lower()
            if confirm in ['yes','y']:
                return 2
            else:
                return cleanChoice(path,noprompt)
        elif choice in ['yes','y']:
            return 1
        else:
            return 0
    if not isinstance(paths,list):
        paths = [paths]
    choice = 0
    while len(paths) > 0:
        path = paths.pop()
        if not choice == 2:
            choice = cleanChoice(path,noprompt)
        if not choice == 0:
            logger.info('cleaning %s...'%path)
            for f in os.listdir(path):
                fn = os.path.join(path,f)
                if os.path.isfile(fn):
                    os.remove(fn)

def parser_ini(mode=0):
    parser = argparse.ArgumentParser(description='Initialize task: Step1. Initialize folder and input base on the given model; Step2. Generate auto-submit&run workflow script.',conflict_handler='resolve',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-p','--np',type=int,default=2,help='Number of processers used in parallel FDTD/PSTD.')
    parser.add_argument('-m','--mode',choices=['m','z','mz'],default='m',help="Mode: 'm' for multi-offset, 'z' for zero-offset, 'mz' for both multi- and zero-offset.")
    group_method = parser.add_mutually_exclusive_group()
    group_method.add_argument('--forward_method',choices=['fdtd','pstd'],default='fdtd',help="Forward method used in RTM.")
    group_method.add_argument('--fdtd',action='store_const',const='fdtd',dest='forward_method',help='Use finite difference time domain as the forward method.')
    group_method.add_argument('--pstd',action='store_const',const='pstd',dest='forward_method',help='Use pseudo spectral time domain as the forward method.')
    group_prompt = parser.add_mutually_exclusive_group()
    group_prompt.add_argument('-y',action='store_const',const='y',dest='noprompt',help="Input 'y' in all input prompts with no disturbing.")
    group_prompt.add_argument('-n',action='store_const',const='n',dest='noprompt',help="Input 'n' in all input prompts with no disturbing.")

    group1 = parser.add_argument_group(title='model_em.py parameter', description='Parameters only used in step1.')
    group1.add_argument('--model',default=os.path.join(modelpath,'model.mat'),help='Path to the .mat file stores all the parameters for simulated data for RTM.')
    group1.add_argument('-f','--freq',type=float,default=-1.0,help='Appoint the main frequency (MHz) of the source to replace the one in the model file.')
    group1.add_argument('--dx_src',type=float,default=-1.0,help='Appoint the x space interval (m) of the source to replace the one in the model file.')
    group1.add_argument('--dy_src',type=float,default=-1.0,help='Appoint the y space interval (m) of the source to replace the one in the model file.')
    group1.add_argument('--dx_rec',type=float,default=-1.0,help='Appoint the x space interval (m) of the receiver to replace the one in the model file.')
    group1.add_argument('--dy_rec',type=float,default=-1.0,help='Appoint the y space interval (m) of the receiver to replace the one in the model file.')
    group1.add_argument('--no_gen_model',action='store_const',const=True,default=False,help="Just modifiy parameters without generating model as generating model needs a lot of time.")
    group1.add_argument('--half_span',type=int,default=-1,help='Span the point source to a (1+2*half_span)x(1+2*half_span) source while forwarding. It helps to reduce the sideslobe in FDTD forwarding.')
    group2 = parser.add_argument_group(title='genbatch.py parameter', description='Parameters only used in step2.')
    if mode == 2:
        group2.add_argument('-t','--taskname',required=True,default='default',dest='dirname',help="work directory (basename) of the task.")
        group2.add_argument('-s','--src_num',type=int,required=True,default=1,help="The total number of shot gathers (sources).")
    group2.add_argument('--server',choices=['local','freeosc','x3850'],default='local',help="Where to run the code.")
    group2.add_argument('--steps',type=str,choices='gfbizc',nargs='+',help="Manually Select which steps are involved. 'g' for generate data; 'f' for source(forward) wavefield; 'b' for receiver(backward) wavefield; 'i' for cross-corrlation image condition; 'z' for cal zero-offset backward wvf(image condition); 'c' for clean middle result. Default: Use 'g f b i c' for a whole multi-offset workflow, 'g f z c' for a whole zero-offset workflow, 'g f b i z c' for both.")
    group_jobcap = group2.add_mutually_exclusive_group()
    group_jobcap.add_argument('--max_job',type=int,default=2,help='Max number of jobs should the server handle at the same time.')
    group_jobcap.add_argument('--max_cpu',type=int,help='Max number of the cpu used by the server. max_job=floor(max_cpu/np).')

    return parser