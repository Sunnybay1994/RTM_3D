#!/usr/bin/env python
import numpy as np
#import mayavi.mlab as mlab
import scipy.linalg as sllg
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import os,sys,shutil,logging,getopt,datetime
import glob

#logger
today = datetime.date.today()
logger = logging.getLogger('model_em_gly.py')
logger.setLevel(logging.INFO) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
fh = logging.FileHandler(os.path.join('log','model_em_gly-' + today.strftime('%Y%m%d') + '.log'))
fh.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('(%(process)5d)%(asctime)s-%(levelname)s: %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)

def cp(f1,f2):
    for f in glob.glob(r'%s'%f1):
        # logger.debug('cp %s %s'%(f,f2))
        shutil.copy(f,f2)

def cleanfiles(paths):

    def cleanChoice(path):
        choice = input('Clear %s? A/Y/(N)'%path).lower()
        if choice in ['all', 'a']:
            confirm = input("Clear all files in '%s'? Y/(N)" % "' & '".join(paths+[path])).lower()
            if confirm in ['yes','y']:
                return 2
            else:
                return cleanChoice(path)
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
            choice = cleanChoice(path)
        if not choice == 0:
            logger.info('cleaning %s...'%path)
            for f in os.listdir(path):
                fn = os.path.join(path,f)
                if os.path.isfile(fn):
                    os.remove(fn)

def eps_sig_mu(meps=1,meps_bg=False,msig=1e-11,msig_bg=False,mmiu=1,mmiu_bg=False):
    logger.info('Generating model...')

    for ii in range(NUM_OF_PROCESS):

        if not isinstance(meps_bg,bool):
            feps_STD = open(os.path.join(std_indir,'eps.in_' + str(ii).zfill(3)), 'w')
        if not isinstance(msig_bg,bool):
            fsig_STD = open(os.path.join(std_indir,'sig.in_' + str(ii).zfill(3)), 'w')
        if not isinstance(mmiu_bg,bool):
            fmiu_STD = open(os.path.join(std_indir,'mu.in_' + str(ii).zfill(3)), 'w')


        feps = open(os.path.join(indir,'eps.in_' + str(ii).zfill(3)), 'w')
        fsig = open(os.path.join(indir,'sig.in_' + str(ii).zfill(3)), 'w')
        fmiu = open(os.path.join(indir,'mu.in_' + str(ii).zfill(3)), 'w')


        logger.info('file: %s'%(feps))
        logger.info('rangex: %d~%d,%d'%(rangex[ii], rangex[ii + 1], rangex[ii + 1] - rangex[ii]))
        for dumx in np.arange(rangex[ii] - order, rangex[ii + 1] + order):

            # for dumx in x_axis:

            if dumx < 0:
                logger.info('begin: %d'%dumx) 
                dumx = 0

            if dumx >= rangex[NUM_OF_PROCESS]:
                logger.info('end: %d'%dumx) 
                dumx = rangex[NUM_OF_PROCESS] - 1

            # Modified by mbw at 20180607
            #print('dumx: ',dumx)
            dumx = int(float(dumx))
            #print(dumx)
            # End Modify

            # initialize model
            eps = np.ones((ny, nz))
            sig = np.ones((ny, nz))
            miu = np.ones((ny, nz))

            if not isinstance(meps,(int,float)):
                eps = meps[dumx,:,:]
            else:
                eps[:,:] = meps
            np.savetxt(feps, eps, fmt='%.3g')

            if not isinstance(msig,(int,float)):
                sig = msig[dumx,:,:]
            else:
                sig[:,:] = msig
            np.savetxt(fsig, sig, fmt='%.3g')

            if not isinstance(mmiu,(int,float)):
                miu = mmiu[dumx,:,:]
            else:
                miu[:,:] = mmiu
            np.savetxt(fmiu, miu, fmt='%.3g')

            if not isinstance(meps_bg,bool):
                eps_STD = np.ones((ny, nz))
                if not isinstance(meps_bg,(int,float)):
                    eps_STD = meps_bg[dumx,:,:]
                else:
                    eps_STD[:,:] = meps_bg
                np.savetxt(feps_STD, eps_STD, fmt='%.3g')

            if not isinstance(msig_bg,bool):
                sig_STD = np.ones((ny, nz))
                if not isinstance(msig_bg,(int,float)):
                    sig_STD = msig_bg[dumx,:,:]
                else:
                    sig_STD[:,:] = msig_bg
                np.savetxt(fsig_STD, sig_STD, fmt='%.3g')

            if not isinstance(mmiu_bg,bool):
                miu_STD = np.ones((ny, nz))
                if not isinstance(mmiu_bg,(int,float)):
                    miu_STD = mmiu_bg[dumx,:,:]
                else:
                    miu_STD[:,:] = mmiu_bg
                np.savetxt(fmiu_STD, miu_STD, fmt='%.3g')
            

        feps.close()
        fsig.close()
        fmiu.close()
        if not isinstance(meps_bg,bool):
            feps_STD.close()
        if not isinstance(msig_bg,bool):
            fsig_STD.close()
        if not isinstance(mmiu_bg,bool):
            fmiu_STD.close()

    
    # if isinstance(meps_bg,bool):
    #     cp(os.path.join(indir, 'eps.in*'), std_indir)
    # if isinstance(msig_bg,bool):
    #     cp(os.path.join(indir, 'sig.in*'), std_indir)
    # if isinstance(mmiu_bg,bool):
    #     cp(os.path.join(indir, 'mu.in*'), std_indir)

    # cp(os.path.join(std_indir, 'eps.in*'), rtm_indir)
    # cp(os.path.join(std_indir, 'sig.in*'), rtm_indir)
    # cp(os.path.join(std_indir, 'mu.in*'), rtm_indir)

    return 0


##############################################################################
def finddt(epmin, mumin, dx, dy, dz):
    epmin = epmin * ep0
    mumin = mumin * mu0
    dtmax = 6.0 / 7.0 * \
        np.sqrt(epmin * mumin / (1.0 / dx ** 2 + 1.0 / dy ** 2 + 1.0 / dz ** 2))
    logger.info("max dt = %fns"%(dtmax/1e-9)) 
    return dtmax


def finddx(epmax, mumax, fmax):
    epmax = epmax * ep0
    mumax = mumax * mu0
    wlmin = 1 / (fmax * np.sqrt(epmax * mumax))
    dxmax = wlmin
    logger.info("max dx = %fm"%dxmax) 
    return dxmax


def blackharrispulse(fmax, dt):
    a = [0.35322222, -0.488, 0.145, -0.010222222]
    T = 1.14 / fmax
    t = np.arange(0, T, dt)
    window = np.zeros(np.size(t))
    for n in range(4):
        window = window + a[n] * np.cos(2 * n * np.pi * t / T)

    window[t >= T] = 0
    p = window
    p = window[:] - np.append(0, window[:-1])
    p = p / np.max(np.abs(p))
    plt.plot(p)
    plt.savefig(os.path.join(workdir,'blackharrispulse.png'))
    return p


def gaussian(dt, t0, nt):
    global t
    T = 10 * t0
    t = linspace(0, T, nt)
    p = 10 ** 5 * exp(-((t - 3 * t0) / t0) ** 2)
    plot(t, p)
    savefig(os.path.join(workdir,'gaussian.png'))
    return p


def ricker(f, length, dt):
    t = np.linspace(-length / 2, (length - dt) / 2, length / dt)
    y = (1.0 - 2.0 * (np.pi ** 2) * (f ** 2) * (t ** 2)) * \
        np.exp(-(np.pi ** 2) * (f ** 2) * (t ** 2))
    try:
        plt.figure()
    except Exception as e:
        logger.error(e)
    else:
        plt.plot(t, y)
        plt.savefig(os.path.join(workdir,'ricker.png'))
    return y


def check_dx(srcpulse):
    n = round(2 ** np.ceil(np.log2(len(srcpulse))))
    freqs = np.linspace(0, 1 / dt / 2, n / 2 + 1)
    sp = np.fft.rfft(srcpulse, n) / n
    W = abs(sp)
    fmax2 = max(freqs[W > max(W) / 10.0])
    logger.info("Src's main frequency: %fMHz" % (freqs[np.argmax(W)]/1e6))
    # logger.info("!!check dx again (src_fmax(within 90%% of max amplitude)=%fMHz):"%(fmax2/1e6)) 
    dx_max = finddx(epmax, mumax, fmax2)
    try:
        plt.figure()
        plt.plot(freqs, W)
        plt.title('frequency spectrum of source')
        plt.savefig(os.path.join(workdir,'spectral_src.png'))
    except Exception as e:
        logger.error(e)
    return dx_max

def distance(x, y, z, x0, y0, z0):
    return np.sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2)


def par():
    content = []
    content.append("#dx dy dz dt\n")
    content.append("%e %e %e %e\n" % (dx, dy, dz, dt))
    content.append("#nx ny nz nt\n")
    content.append("%d %d %d %d\n" % (nx, ny, nz, nt))
    content.append("#nt of src\n")
    content.append("%d\n" % (nt_src))
    content.append("#output time step and space step of wavefield\n")
    content.append("%d %d\n" % (outstep_t_wavefield, outstep_x_wavefield))
    content.append("#output step of slice\n")
    content.append("%d\n" % (outstep_slice))
    content.append("#npml x y z\n")
    content.append("%d %d %d\n" % (npmlx, npmly, npmlz))
    content.append("#pml m kapxmax kapymax kapzmax alpha\n")
    content.append("4 5 5 5 0.00\n")
    content.append("#location of src\n")
    content.append("src.in\n")
    content.append("#location of rec\n")
    content.append("rec.in\n")
    content.append("#epsilon file\n")
    content.append("eps.in\n")
    content.append("#mu file\n")
    content.append("mu.in\n")
    content.append("#sigma file\n")
    content.append("sig.in\n")
    content.append("#slices file\n")
    content.append("slice.in\n")
    with open(os.path.join(indir, 'par.in'), 'w') as fpar:
        fpar.write(''.join(content))

    # cp(os.path.join(indir, 'par.in'), std_indir)
    # if is_zRTM:
    #     content[1] = "%e %e %e %e\n" % (dx, dy, dz, dt/2)
    #     with open(os.path.join(rtm_indir, 'par.in'), 'w') as fpar:
    #         fpar.write(''.join(content))
    # else:
    #     cp(os.path.join(indir, 'par.in'), rtm_indir)


def X_partition(nx, NUM_OF_PROCESS):
    # global
    rangex = np.zeros(NUM_OF_PROCESS + 1)
    nxSize = np.zeros(NUM_OF_PROCESS)
    rangex[0] = 0

    for i in range(1, NUM_OF_PROCESS):
        if i <= (nx % NUM_OF_PROCESS):
            rangex[i] = rangex[i - 1] + (nx // NUM_OF_PROCESS + 1) # change '//' to '/'
        else:
            rangex[i] = rangex[i - 1] + (nx // NUM_OF_PROCESS) # change '//' to '/'
        nxSize[i - 1] = rangex[i] - rangex[i - 1] + 2 * order
    rangex[NUM_OF_PROCESS] = nx
    nxSize[NUM_OF_PROCESS - 1] = rangex[NUM_OF_PROCESS] - \
        rangex[NUM_OF_PROCESS - 1] + 2 * order
    return rangex, nxSize


def islice(sxl,syl,szl):
    if not isinstance(sxl,list):
        sxl = [sxl]
    if not isinstance(syl,list):
        syl = [syl]
    if not isinstance(szl,list):
        szl = [szl]


    nslicex = len(sxl)
    nslicey = len(syl)
    nslicez = len(szl)

    slicex = [[int(nx),'Ey'] for nx in sxl]
    slicey = [[int(ny),'Ey'] for ny in syl]
    slicez = [[int(nz),'Ey'] for nz in szl]
    
    logger.info("slicex(%d):%s, slicey(%d):%s, slicez(%d):%s"%(nslicex,slicex,nslicey,slicey,nslicey,slicez)) 

    fn_slice = os.path.join(indir, 'slice.in')
    with open(fn_slice, 'w') as fslice:
        fslice.write("%d %d %d\n" % (nslicex, nslicey, nslicez))
        for i in range(len(slicex)):
            fslice.write("%d %s\n" % (slicex[i][0], slicex[i][1]))
        for i in range(len(slicey)):
            fslice.write("%d %s\n" % (slicey[i][0], slicey[i][1]))
        for i in range(len(slicez)):
            fslice.write("%d %s\n" % (slicez[i][0], slicez[i][1]))
    
    # cp(fn_slice, std_indir)
    # cp(fn_slice, rtm_indir)
    # if Dum_RTM == 1:
    #     # updated by mbw at 20190415
    #     shutil.copy('./Input/slice.in', './STD/Input/')
    #     shutil.copy('./Input/slice.in', './RTM/Input/')

##############################################################################
if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:", ["workdir=","no_gen_model"])
    except getopt.GetoptError as err:
        # print help information and exit:
        logger.error(err)  # will print something like "option -a not recognized"
        # usage()
        sys.exit(2)

    gen_model = True
    workdir = os.path.join('tasks','gly_400MHz_linear_smooth')
    fsrc = '2src0_400MHz.in_0000'
    for o, a in opts:
        if o in ('-d','--workdir'):
            if a.strip() == '1':
                mode_hz = 1
                workdirn = 'gly_100MHz'
                fsrc = '100MHz_8_0.5_src.in_0000'
            elif a.strip() == '1a':
                mode_hz = 1
                workdirn = 'gly_100MHz_agc'
                fsrc = '100MHz_8_0.5_src_agc.in_0000'
            elif a.strip() == '4':
                mode_hz = 4
                workdirn = 'gly_400MHz'
                fsrc = '400MHz_8_0.5_src.in_0000'
            elif a.strip() == '4a':
                mode_hz = 4
                workdirn = 'gly_400MHz_agc'
                fsrc = '400MHz_8_0.5_src_agc.in_0000'
            else:
                workdirn = a
            workdir = os.path.join('tasks',workdirn)
        elif o in ('--no_gen_model'):
            gen_model = False
        else:
            assert False, "unhandled option"
    




    #####################################################################
    mu0 = 1.2566370614e-6
    ep0 = 8.8541878176e-12

    epmax = 3.0
    mumax = 4.0
    epmin = epmax
    mumin = 1.0
    fmax = 400e6  # Hz

    dx_max = finddx(epmax, mumax, fmax)
    dx = 0.05
    dy = 0.1
    dz = 0.02
    logger.info("dx=%g, dy=%g, dz=%g"%(dx, dy, dz)) 
    assert np.max([dx,dy,dz]) < dx_max, 'dx,dy,dz too big!!!'

    nz_air = 10
    nx = 1200+10
    if mode_hz == 4:
        ny = 125+10 # 400MHz
        nz = round(6/dz)+ nz_air# 400MHz
        dt = 0.0587e-9# 400MHz
        nt = 952+1#400MHz
    elif mode_hz == 1:
        ny = 21+10 # 100MHz
        nz = round(12/dz)+ nz_air# 100MHz
        dt = 0.1466e-9# 100MHz
        nt = 924+1#100MHz
    logger.info('nx=%d, ny=%d, nz=%d'%(nx, ny, nz)) 

    npmlx = 8
    npmly = 8
    npmlz = 8
    nt_src = 50 # useless

    outstep_t_wavefield = 4
    outstep_x_wavefield = 2
    outstep_slice = 4

    dt_max = finddt(epmin, mumin, dx, dy, dz)
    # nt better be k*outstep_t_wavefield+1(k is integer), so that forward wavefield and 
    # backward wavefield will coincide perfectly when doing cross-correlation.
    # nt += outstep_t_wavefield + 1 - nt%outstep_t_wavefield #useless
    logger.info("dt=%g, nt=%d"%(dt, nt)) 
    assert dt < dt_max, 'dt too big!!! (%g>%g)'%(dt,dt_max)

    ###########################################################################
    workdir = workdir + '_' + str(epmax) + '_' + str(dx) + '_' + str(dy) + '_' + str(dz)
    logger.info('workdir="%s"'%(workdir))

    indir = os.path.join(workdir,'Input')
    outdir = os.path.join(workdir,'Output')
    dirlist = [workdir,indir,outdir]
    for idir in dirlist:
        if not os.path.exists(idir):
            os.mkdir(idir)
        else:
            cleanfiles([idir])


    # added by mbw at 20190415
    shutil.copy("FDTD_MPI",workdir)
    ###########################################################################

    NUM_OF_PROCESS = 4
    order = 2 # num of interchange layers of each process
    logger.info("NUM_OF_PROCESS: %d"%(NUM_OF_PROCESS)) 
    if gen_model:
        rangex, nxSize = X_partition(nx, NUM_OF_PROCESS)
        eps_sig_mu(epmax)

    shutil.copy(os.path.join('src_rec','rec.in'), os.path.join(workdir,'Input'))
    shutil.copy(os.path.join('src_rec',fsrc), os.path.join(workdir,'Input','src.in_0000'))
    par()
    islice(round(nx/2),round(ny/2),round(1/dz)+nz_air)

    cp('model_em_gly.py', workdir)

#######################################
    qsubfn = os.path.join(workdir,'sub_0000.sh')
    cwd = os.getcwd()
    workpath = os.path.join(cwd,workdir)
    with open(qsubfn,'w') as fo:
        fo.write('''#!/bin/sh

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
#$ -N rtm_gly

######### set Parallel Environment and CORE numbers
module unload mpi/mpich-x86_64
module load mpi/openmpi-x86_64
#$ -pe openmpi ''' + str(NUM_OF_PROCESS) +'''

echo "Got $NSLOTS slots."
echo "PATH = $PATH"
echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
PYPATH="''' + cwd + '''"
WORKPATH="$PYPATH/''' + workdir + '''"
echo "Current Directory = $WORKPATH"
echo 

#########  execute PROGRAM_NAME
echo "Computing is started at $(date)."

~/software/openmpi-4.0.3/bin/mpiexec -np $NSLOTS -wdir $WORKPATH $WORKPATH/FDTD_MPI 0

echo "Computing is stopped at $(date)."

exit 0
''')