#!/usr/bin/env python
import numpy as np
#import mayavi.mlab as mlab
import scipy.linalg as sllg
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import os,sys,shutil,logging,argparse,datetime
import glob
from writesource import *

#logger
today = datetime.date.today()
logger = logging.getLogger('model_em')
logger.setLevel(logging.INFO) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
fh = logging.FileHandler(os.path.join('log','model_em-' + today.strftime('%Y%m%d') + '.log'))
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

def src_rec(dnx_src,dny_src=False,dnx_rec=False,dny_rec=False,nzp_src=False,nzp_rec=False,nx_src=False,ny_src=False,nx_rec=False,ny_rec=False,nshift=0,marginx=0,marginy=False,marginx_rec=False,marginy_rec=False,half_span=0):
    logger.info('adding source and receiver...')
    logger.info('dnxs=%d,dnys=%d,dnxr=%d,mx=%d,my=%d,mxr=%d,myr=%d,half_span=%d'%(dnx_src,dny_src,dnx_rec,marginx,marginy,marginx_rec,marginy_rec,half_span))
    if not dny_src:
        dny_src = dnx_src
    if not nzp_src: # z position of src 
        nzp_src = nz_air - 1
    if not dnx_rec:
        dnx_rec = dnx_src
        if not dny_rec:
            dny_rec = dny_src
    else:
        if not dny_rec:
            dny_rec = dnx_rec
    if not nzp_rec: # z position of rec
        nzp_rec = nzp_src
    if not marginy:
        marginy = marginx
    if not marginx_rec:
        marginx_rec = marginx
    if not marginy_rec:
        marginy_rec = marginx_rec

    global src, rec
    def ant_pos(dnx_ant,dny_ant,nz_pos,nx_ant,ny_ant,nshift,marginx,marginy):
        ant = []
        dnx_ant = dnx_ant
        dny_ant = dny_ant
        nz_pos = nz_pos

        if not nx_ant:
            nx_ant = int((nx - 2*marginx) // dnx_ant)
            if nx_ant % 2 == 0:
                nx_ant -= 1
        if not ny_ant:
            ny_ant = int((ny - 2*marginy) // dny_ant)
            if ny_ant % 2 == 0:
                ny_ant -= 1
        
        for i in range(-(nx_ant-1)//2,(nx_ant+1)//2):
            dumx = nx0 + i*dnx_ant
            for j in range(-(ny_ant-1)//2,(ny_ant+1)//2):
                dumy = ny0 + j*dny_ant
                ant.append([dumx + nshift, dumy + nshift, nz_pos, 'Ey', 1])
        return ant

    src = ant_pos(dnx_src,dny_src,nzp_src,nx_src,ny_src,nshift,marginx,marginy)
    nsrc = len(src)
    logger.info("nsrc: %d"%nsrc) 

    rec = ant_pos(dnx_rec,dny_rec,nzp_rec,nx_rec,ny_rec,nshift,marginx_rec,marginy_rec)
    nrec = len(rec)
    logger.info("nrec: %d"%nrec) 

    for i in range(nsrc):
        fn_src = os.path.join(indir,'src.in_' + str(i).zfill(4))
        extend_and_write_one_source(fn_src, src[i], srcpulse,xhalfspan=half_span,yhalfspan=half_span)

    with open(os.path.join(indir,'rec.in'), 'w') as frec:
        frec.write("%d\n" % (nrec))
        for i in range(nrec):
            frec.write("%d,%d,%d,%s\n" %
                       (rec[i][0], rec[i][1], rec[i][2], rec[i][3]))

    cp(os.path.join(indir,'src.in*'),std_indir)
    cp(os.path.join(indir,'rec.in'),os.path.join(std_indir,'rec.in'))
    if 'z' in mode:
        with open(os.path.join(rtm0_indir,'rec.in'),'w+') as fo:
            fo.write('1\n')
            fo.write('%d,%d,%d,Ey\n'%(nx0,ny0,nz_air-2))
    if 'm' in mode:
        with open(os.path.join(rtm_indir,'rec.in'),'w+') as fo:
            fo.write('1\n')
            fo.write('%d,%d,%d,Ey\n'%(nx0,ny0,nz_air-2))
            
    return nsrc,nrec


def eps_sig_mu(meps=1,meps_bg=False,msig=1e-5,msig_bg=False,mmiu=1,mmiu_bg=False,pstd=False):
    logger.info('Generating model...')
    if pstd:
        pnum = 1
    else:
        pnum = NUM_OF_PROCESS

    for ii in range(pnum):

        if pstd:
            suffix = ''
        else:
            suffix = '_' + str(ii).zfill(3)

        if not isinstance(meps_bg,bool):
            feps_STD = open(os.path.join(std_indir,'eps.in' + suffix), 'w')
        if not isinstance(msig_bg,bool):
            fsig_STD = open(os.path.join(std_indir,'sig.in' + suffix), 'w')
        if not isinstance(mmiu_bg,bool):
            fmiu_STD = open(os.path.join(std_indir,'mu.in' + suffix), 'w')


        feps = open(os.path.join(indir,'eps.in' + suffix), 'w')
        fsig = open(os.path.join(indir,'sig.in' + suffix), 'w')
        fmiu = open(os.path.join(indir,'mu.in' + suffix), 'w')


        if pstd:
            dumx_range = range(nx)
        else:
            logger.info('file: %s'%(feps))
            logger.info('rangex: %d~%d,%d'%(rangex[ii], rangex[ii + 1], rangex[ii + 1] - rangex[ii]))
            dumx_range = np.arange(rangex[ii] - order, rangex[ii + 1] + order)

        for dumx in dumx_range:

            # for dumx in x_axis:

            if dumx < 0:
                continue
                # logger.info('begin: %d'%dumx) 
                # dumx = 0
            elif dumx >= nx:
                continue
                # logger.info('end: %d'%dumx) 
                # dumx = rangex[pnum] - 1

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

    
    if isinstance(meps_bg,bool):
        cp(os.path.join(indir, 'eps.in*'), std_indir)
    if isinstance(msig_bg,bool):
        cp(os.path.join(indir, 'sig.in*'), std_indir)
    if isinstance(mmiu_bg,bool):
        cp(os.path.join(indir, 'mu.in*'), std_indir)

    if 'z' in mode:
        cp(os.path.join(std_indir, 'eps.in*'), rtm0_indir)
        cp(os.path.join(std_indir, 'sig.in*'), rtm0_indir)
        cp(os.path.join(std_indir, 'mu.in*'), rtm0_indir)
    if 'm' in mode:
        cp(os.path.join(std_indir, 'eps.in*'), rtm_indir)
        cp(os.path.join(std_indir, 'sig.in*'), rtm_indir)
        cp(os.path.join(std_indir, 'mu.in*'), rtm_indir)

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
    t = np.linspace(0, T, nt)
    p = 10 ** 5 * np.exp(-((t - 3 * t0) / t0) ** 2)
    plt.plot(t, p)
    plt.savefig(os.path.join(workdir,'gaussian.png'))
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
    n = int(round(2 ** np.ceil(np.log2(len(srcpulse)))))
    freqs = np.linspace(0, int(1 / dt / 2), int(n / 2) + 1)
    sp = np.fft.rfft(srcpulse, n) / n
    W = abs(sp)
    fmax2 = max(freqs[W > max(W) / 10.0])
    logger.info("Src's main frequency: %fMHz" % (freqs[np.argmax(W)]/1e6))
    logger.info("!!check dx again (src_fmax(within 90%% of max amplitude)=%fMHz):"%(fmax2/1e6)) 
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
    content.append("%e,%e,%e,%e\n" % (dx, dy, dz, dt))
    content.append("#nx ny nz nt\n")
    content.append("%d,%d,%d,%d\n" % (nx, ny, nz, nt))
    content.append("#nt of src\n")
    content.append("%d\n" % (nt_src))
    content.append("#output time step and space step of wavefield\n")
    content.append("%d,%d\n" % (outstep_t_wavefield, outstep_x_wavefield))
    content.append("#output step of slice\n")
    content.append("%d\n" % (outstep_slice))
    content.append("#npml x y z\n")
    content.append("%d,%d,%d\n" % (npmlx, npmly, npmlz))
    content.append("#pml m kapxmax kapymax kapzmax alpha\n")
    content.append("4,5,5,5,0.00\n")
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
    # content.append("#slices(for pstd)\n")
    # content.append("%d,%d,%d\n"%(slx,sly,slz))
    # content.append("#nthreads(for pstd)\n")
    # content.append("%d\n"%nthreads)
    with open(os.path.join(indir, 'par.in'), 'w') as fpar:
        fpar.write(''.join(content))

    cp(os.path.join(indir, 'par.in'), std_indir)
    if 'z' in mode:
        content[1] = "%e,%e,%e,%e\n" % (dx, dy, dz, dt/2)
        with open(os.path.join(rtm0_indir, 'par.in'), 'w') as fpar:
            fpar.write(''.join(content))
    if 'm' in mode:
        cp(os.path.join(indir, 'par.in'), rtm_indir)


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
        fslice.write("%d,%d,%d\n" % (nslicex, nslicey, nslicez))
        for i in range(len(slicex)):
            fslice.write("%d,%s\n" % (slicex[i][0], slicex[i][1]))
        for i in range(len(slicey)):
            fslice.write("%d,%s\n" % (slicey[i][0], slicey[i][1]))
        for i in range(len(slicez)):
            fslice.write("%d,%s\n" % (slicez[i][0], slicez[i][1]))
    
    cp(fn_slice, std_indir)
    if 'z' in mode:
        cp(fn_slice, rtm0_indir)
    if 'm' in mode:
        cp(fn_slice, rtm_indir)

    return sxl[0],syl[0],szl[0]


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



##############################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate model for RTM',conflict_handler='resolve')
    parser.add_argument('--model',default=os.path.join('Model','model.mat'),help='An .mat file stores all the parameters for simulated data for RTM.')
    parser.add_argument('-f','--freq',type=float,default=-1.0,help='Appoint the main frequency (MHz) of the source to replace the one in the model file.')
    parser.add_argument('--dx_src',type=float,default=-1.0,help='Appoint the x space interval (m) of the source to replace the one in the model file.')
    parser.add_argument('--dy_src',type=float,default=-1.0,help='Appoint the y space interval (m) of the source to replace the one in the model file.')
    parser.add_argument('--dx_rec',type=float,default=-1.0,help='Appoint the x space interval (m) of the receiver to replace the one in the model file.')
    parser.add_argument('--dy_rec',type=float,default=-1.0,help='Appoint the y space interval (m) of the receiver to replace the one in the model file.')
    parser.add_argument('-m','--mode',choices=['m','z','mz'],default='m',help="Mode: 'm' for multi-offset, 'z' for zero-offset, 'mz' for both multi- and zero-offset.")
    # parser.add_argument('--forward_method',choices=['fdtd','pstd'],default='fdtd',help="Forward method used in RTM.")
    parser.add_argument('--fdtd',action='store_const',const='fdtd',dest='forward_method',default='fdtd',help='Use finite difference time domain as the forward method.')
    parser.add_argument('--pstd',action='store_const',const='pstd',dest='forward_method',default='fdtd',help='Use pseudo spectral time domain as the forward method.')
    parser.add_argument('--no_gen_model',action='store_const',const=True,default=False,help="Just modifiy parameters without generating model as generating model needs a lot of time.")
    parser.add_argument('--half_span',type=int,default=-1,help='Span the point source to a (1+2*half_span)x(1+2*half_span) source while forwarding. It helps to reduce the sideslobe in FDTD forwarding.')
    parser.add_argument('--server',choices=['local','freeosc','x3850'],default='local',help="Where to run the code.")
    parser.add_argument('-p','--np',type=int,default=-1,help='Number of processers/threads used in parallel FDTD/PSTD. Default: 6/12 for local, 8/16 for x3850 and freeosc.')
    parser.add_argument('-c','--job_cap',type=int,default=-1,help='How may shot gathers (sources) should the server handle at the same time. Default: 1 for local, 8 for x3850 and 32 for freeosc.')
    parser.add_argument('-y',action='store_const',const='y',dest='noprompt',default=False,help="Input 'y' in all input prompts with no disturbing.")
    parser.add_argument('-n',action='store_const',const='n',dest='noprompt',default=False,help="Input 'n' in all input prompts with no disturbing.")
    parser.add_argument('--steps',type=str,default='gfbi',help="Select which steps are involved. 'g' for generate data; 'f' for source(forward) wavefield; 'b' for receiver(backward) wavefield; 'i' for imaging.")
    args = parser.parse_args()

    model = args.model
    mode = args.mode
    forward_method = args.forward_method
    server = args.server
    gen_model = not args.no_gen_model
    noprompt = args.noprompt
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
    logger.info('pnum=%d, job_capacity=%d'%(pnum,job_cap))


    ### load model ###
    try:
        dic_model = sio.loadmat(model)
    except Exception as e:
        logger.error(e)
        raise e
    else:
        modelname = str(dic_model['modelname'][0])
    finally:
        pass
    ### load model end ###

    ### parameter ###
    if args.freq > 0:
        freq = args.freq * 1e6
    else:
        freq = float(dic_model['freq_src'])
    
    if args.half_span < 0:
        try:
            half_span = int(dic_model['src_span'])
        except Exception as e:
            logger.warning('%s, half_span set_to 0.'%e)
            half_span = 0
    else:
        half_span = args.half_span


    mu0 = 1.2566370614e-6
    ep0 = 8.8541878176e-12

    epmin = 1.0
    mumin = 1.0
    epmax = 15.0
    mumax = 1.0
    fmax = freq #Hz
    ### parameter end ###


    ### src & rec ###
    if args.dx_src > 0:
        dx_src = args.dx_src
    else:
        dx_src = float(dic_model['dx_src'])
    if args.dy_src > 0:
        dy_src = args.dy_src
    else:
        dy_src = float(dic_model['dy_src'])
    if args.dx_rec > 0:
        dx_rec = args.dx_rec
    else:
        dx_rec = float(dic_model['dx_rec'])
    if args.dy_rec > 0:
        dy_rec = args.dy_rec
    else:
        dy_rec = float(dic_model['dy_rec'])
    ### src & rec end ###

    ### init workdir ###
    dir_suffix = ''
    if forward_method == 'pstd':
        dir_suffix += '_pstd_%d'%pnum
    elif forward_method == 'fdtd':
        dir_suffix += '_fdtd_%d'%pnum
    if mode=='z':
        dir_suffix += '_0o'

    dirname = '%s_%dMHz_%gm_%gm%s'%(modelname,freq/1e6,dx_src,dx_rec if dx_rec!=dx_src else dy_src, dir_suffix)
    workdir = os.path.join('tasks',dirname)
    logger.info('workdir="%s"'%(workdir))
    ### init workdir end ###

    ### directories ###
    indir = os.path.join(workdir,'Input')
    outdir = os.path.join(workdir,'Output')
    std_dir = os.path.join(workdir,'STD')
    std_indir = os.path.join(std_dir,'Input')
    std_outdir = os.path.join(std_dir,'Output')
    rtm_dir = os.path.join(workdir,'RTM')
    rtm_indir = os.path.join(rtm_dir,'Input')
    rtm_outdir = os.path.join(rtm_dir,'Output')
    rtm0_dir = os.path.join(workdir,'RTM0')
    rtm0_indir = os.path.join(rtm0_dir,'Input')
    rtm0_outdir = os.path.join(rtm0_dir,'Output')
    dirlist1 = [workdir,std_dir]
    dirlist2 = [indir,std_indir]
    dirlist3 = [outdir,std_outdir]
    if 'm' in mode:
        dirlist1 += [rtm_dir]
        dirlist2 += [rtm_indir]
        dirlist3 += [rtm_outdir]
    if 'z' in  mode:
        dirlist1 += [rtm0_dir]
        dirlist2 += [rtm0_indir]
        dirlist3 += [rtm0_outdir]

    cleanlist = []
    for dirlist in [dirlist1,dirlist2,dirlist3]:
        for idir in dirlist:
            if not os.path.exists(idir):
                os.mkdir(idir)
            elif not idir in dirlist1:
                cleanlist += [idir]
    cleanfiles(cleanlist,noprompt)

    if forward_method == 'fdtd':
        forward_fn = "FDTD_MPI.exe"
    elif forward_method == 'pstd':
        forward_fn = "PSTD.exe"
        shutil.copy(forward_fn,workdir)
        shutil.copy(forward_fn,std_dir)
        if 'm' in mode:
            shutil.copy(forward_fn,rtm_dir)
        if 'z' in mode:
            shutil.copy(forward_fn,rtm0_dir)
    ### directories end ###


    ### gird  parameter ###
    npmlx = dic_model['npmlx']
    npmly = dic_model['npmly']
    npmlz = dic_model['npmlz']

    outstep_t_wavefield = dic_model['outstep_t_wavefield']
    outstep_x_wavefield = dic_model['outstep_x_wavefield']
    outstep_slice = dic_model['outstep_slice']

    dx_max = finddx(epmax, mumax, fmax)
    dx = float(dic_model['dx'])
    dy = float(dic_model['dy'])
    dz = float(dic_model['dz'])
    logger.info("dx=%f, dy=%f, dz=%f"%(dx, dy, dz))
    assert np.max([dx,dy,dz]) < dx_max, 'dx,dy,dz too big!!!'

    nx = int(dic_model['nx'])
    ny = int(dic_model['ny'])
    nx0 = nx/2 # middle nx
    ny0 = ny/2 # middle ny
    nz_air = int(dic_model['nz_air'])
    nz = int(dic_model['nz'])
    logger.info('nx=%d, ny=%d, nz=%d(include nz_air =%d)'%(nx, ny, nz, nz_air))

    dt_max = finddt(epmin, mumin, dx, dy, dz)
    dt = float(dic_model['dt'])
    nt = round(float(dic_model['T'])/dt)
    # nt better be k*outstep_t_wavefield+1(k is integer), so that forward wavefield and 
    # backward wavefield will coincide perfectly when doing cross-correlation.
    nt += outstep_t_wavefield + 1 - nt%outstep_t_wavefield
    logger.info("dt: %fns"%(dt/1e-9)) 
    assert dt < dt_max, 'dt too big!!! (%f>%f)'%(dt,dt_max)
    ### gird paraeter end ###

    ### source ###
    srcpulse = blackharrispulse(fmax, dt)
    # srcpulse = ricker(fmax,4*1/fmax,dt)
    nt_src = len(srcpulse)
    logger.info("nt=%d, nt_src=%d"%(nt, nt_src)) 
    dx_max = check_dx(srcpulse)
    assert np.max([dx,dy,dz]) < dx_max, 'dx,dy,dz too big!!!'
    ### source end ###

    ### generate src & rec ###
    dnx_src = round(dx_src/dx)
    dny_src = round(dy_src/dy)
    dnx_rec = round(dx_rec/dx)
    dny_rec = round(dy_rec/dx)
    mx = float(dic_model['src_margin_nx'])
    my = float(dic_model['src_margin_ny'])
    mxr = float(dic_model['rec_margin_nx'])
    myr = float(dic_model['rec_margin_ny'])
    [nsrc,nrec] = src_rec(dnx_src,dny_src,dnx_rec,dny_rec,marginx=mx, marginy=my, marginx_rec=mxr,marginy_rec=myr,half_span=half_span)
    ### generate src & rec end ###

    ### generate model ###
    if gen_model:
        if forward_method == 'fdtd':
            NUM_OF_PROCESS = pnum
        elif forward_method == 'pstd':
            NUM_OF_PROCESS = 1
        order = 2 # num of interchange layers of each process
        logger.info("NUM_OF_PROCESS: %d"%NUM_OF_PROCESS) 
        rangex, nxSize = X_partition(nx, NUM_OF_PROCESS)
        eps_sig_mu(dic_model['ep'],dic_model['ep_bg'],pstd=forward_method=='pstd')
    ### generate model end ###

    slx,sly,slz = islice(dic_model['slicex'][0].tolist(),dic_model['slicey'][0].tolist(),dic_model['slicez'][0].tolist())
    par()

    # backup this file
    cp('model_em.py', workdir)
    # cp(os.path.join('Model',model[:-4]+'_sr'+model[-4:]), workdir) # backup model_sr
    cp(model, workdir) # backup model
    cp(os.path.join('Model','make_model.m'), workdir) # backup make_model

    
    if forward_method == 'pstd':
        forward_method = '--pstd'
    elif forward_method == 'fdtd':
        forward_method = '--fdtd'
    subtxt = 'python batchgen.py -d %s -s %d -p %d -m %s -c %d --server %s --steps %s %s %s'%(dirname,nsrc,pnum,mode,job_cap,server,args.steps,forward_method,'-'+noprompt if noprompt else '')
    logger.info(subtxt)
    os.system(subtxt)