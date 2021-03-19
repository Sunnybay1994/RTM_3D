import os,configparser

class rtm_base(object):
    def __init__(self):
        object.__init__(self)
        _rootpath = os.path.abspath('..')
        for d in [self.logpath, self.taskpath]:
            if not os.path.isdir(d):
                os.mkdir(d)
    @property
    def rootpath(self):
        return self._rootpath
    @property
    def srcpath(self):
        return os.path.join(self.rootpath,'rtm_3d')
    @property
    def modelpath(self):
        return os.path.join(self.srcpath,'make_model')
    @property
    def binpath(self):
        return sos.path.join(self.rootpath,'bin')
    @property
    def taskpath(self):
        return os.path.join(self.rootpath,'tasks')
    @property
    def logpath(self):
        return os.path.join(self.rootpath,'log')
    @property
    def config(self):
        return configparser.ConfigParser()

class task_path(rtm_base):
    _indir = 'Input'
    _outdir = 'Output'
    @property
    def indir(self):
        return self._indir
    @property
    def outdir(self):
        return self._outdir
    
    def __init__(self,taskname,mode):
        rtm_base.__init__(self)
        _taskname = taskname
        _mode = mode
        _workpath = ''
        _configfile = 'config.ini'
        for d in [self.workpath, self.rtmpath, self.stdpath, self.logpath, self.stdpath, self.rtmpath, self.inpath, self.outpath, self.stdinpath, self.stdoutpath, self.rtminpath, self.rtmoutpath]:
            if not os.path.isdir(d):
                os.mkdir(d)
    @property
    def taskname(self):
        return self._taskname
    @property
    def mode(self):
        return self._mode
    @property
    def absworkpath(self):
        return os.path.join(self.taskpath,self.taskname)
    @property
    def workpath(self):
        return self._workpath
    @property
    def configfile(self):
        return os.path.join(self.workpath,self._configfile)
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
    def tlogpath(self):
        return os.path.join(self.workpath,'log')
    @property
    def inpath(self):
        return os.path.join(self.workpath,self.indir)
    @property
    def outpath(self):
        return os.path.join(self.workpath,self.outdir)
    @property
    def stdinpath(self):
        return os.path.join(self.stdpath,self.indir)
    @property
    def stdoutpath(self):
        return os.path.join(self.stdpath,self.outdir)
    @property
    def rtminpath(self):
        return os.path.join(self.rtmpath,self.indir)
    @property
    def rtmoutpath(self):
        return os.path.join(self.rtmpath,self.outdir)
    @property
    def config(self):
        config = rtm_base.config
        config.add_section('task'):
        config.set('task','taskname',self.taskname)
        config.set('task','mode',self.mode)
