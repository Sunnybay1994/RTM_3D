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

<<<<<<< HEAD:model_em.py
def src_rec(dnx_src,dny_src=False,dnx_rec=False,dny_rec=False,nzp_src=False,nzp_rec=False,nx_src=False,ny_src=False,nx_rec=False,ny_rec=False,nshift=0,marginx=0,marginy=False,marginx_rec=False,marginy_rec=False):
    logger.info('adding source and receiver...')
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
=======
def cleanfiles(paths):
>>>>>>> gly:model_em_gly.py

    def cleanChoice(path):
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

def src_rec(Dum_RTM=1):
    global src, rec
<<<<<<< HEAD:model_em.py
    def ant_pos(dnx_ant,dny_ant,nz_pos,nx_ant,ny_ant,nshift,marginx,marginy):
        ant = []
        dnx_ant = int(dnx_ant)
        dny_ant = int(dny_ant)
        nz_pos = int(nz_pos)

        if not nx_ant:
            nx_ant = (nx - 2*marginx) // dnx_ant
            if nx_ant % 2 == 0:
                nx_ant += 1
        if not ny_ant:
            ny_ant = (ny - 2*marginy) // dny_ant
            if ny_ant % 2 == 0:
                ny_ant += 1
        
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
=======

    nshift = 0

    # nx_src = round(nx/4.0)
    # ny_src = round(ny/2.0)
    # nx_src_real = nx_src * dx
    # ny_src_real = ny_src * dx + yaxis2[0]
    # nz_src_real = find_z1(nx_src_real, ny_src_real) - z_axis[0]
    # nz_src = round(nz_src_real/dx)

    # print nx_src_real, ny_src_real, nz_src_real
    # print nx_src, ny_src, nz_src

    # if nz_src_real + z_axis[0]> 0:
    #     src.append([nx_src, ny_src, nz_src, 'Ey',1])
    # nsrc = len(src)
    # print "nsrc: ",nsrc
    # nsrc = 10
    # src = []
    # nx_src =  float(nx-2*npmlx)/(nsrc+1)
    # ny_src =  float(ny-2*npmly)/(nsrc+1)
    # dumx = npmlx
    # dumy = npmly
    # for i in range(nsrc):
    #     dumx += nx_src
    #     dumy = npmly
    #     for i in range(nsrc):
    #         dumy += ny_src
    #         nx_src_real = dumx * dx
    #         ny_src_real = dumy * dx + yaxis2[0]
    #         nz_src_real = find_z1(nx_src_real, ny_src_real) - z_axis[0]
    #         nz_src = round(nz_src_real/dx) + 1
    # if nz_src_real + z_axis[0] -0.5> 0: ##  -0.5 for remove the sources near horizon.
    # mlab.points3d(nx_src_real, ny_src_real, nz_src_real)
    #             src.append([dumx+nshift,dumy+nshift,nz_src,'Ey',1])
    src = []
    delta_x_src = 1.0 # 0.3 # mbw 20180628
    delta_y_src = delta_x_src
    delta_nx_src = delta_x_src / dx
    delta_ny_src = delta_y_src / dy
    nx_src = int(fix((nx - 2*npmlx) * dx / delta_x_src)) - 1
    ny_src = int(fix((ny - 2*npmly) * dy / delta_y_src)) - 1
#    ny_rec = nx_rec
    dumx = npmlx
    dumy = npmly
#    nx_src = ny_src
    for i in range(nx_src):
        dumx += delta_nx_src
        dumy = npmly
        for i in range(ny_src):
            dumy += delta_ny_src
            nx_src_real = dumx * dx
            ny_src_real = dumy * dy
            nz_src_real = (20 - 2) * dz # surface: 0.4m
            nz_src = 20 - 2
            # -0.5 for remove the sources near horizon.
            # if nz_src_real + z_axis[0] - 0.5 > 0:
                # mlab.points3d(nx_src_real, ny_src_real, nz_src_real + z_axis[0],scale_factor = 0.1)
            src.append([dumx + nshift, dumy + nshift, nz_src, 'Ey', 1])     

    nsrc = len(src)
    logger.info("nsrc: %d"%nsrc) 

    # nrec = 50
    # rec = []
    # nx_rec =  float(nx-2*npmlx)/(nrec+1)
    # ny_rec =  float(ny-2*npmly)/(nrec+1)
    # dumx = npmlx
    # dumy = npmly
    # for i in range(nrec):
    #     dumx += nx_rec
    #     dumy = npmly
    #     for i in range(nrec):
    #         dumy += ny_rec
    #         nx_rec_real = dumx * dx
    #         ny_rec_real = dumy * dx + yaxis2[0]
    #         nz_rec_real = find_z1(nx_rec_real, ny_rec_real) - z_axis[0]
    #         nz_rec = round(nz_rec_real/dx) + 1
    #         if nz_rec_real + z_axis[0]> 0:
    # mlab.points3d(nx_rec_real, ny_rec_real, nz_rec_real)
    #             rec.append([dumx+nshift,dumy+nshift,nz_rec,'Ey'])

    rec = []
#    nx_rec =  float(nx-2*npmlx)/(nrec+1)
#    ny_rec =  float(ny-2*npmly)/(nrec+1)
#    nx_rec = ny_rec
    delta_x_rec = 0.5 # 0.3 # mbw 20180628
    delta_y_rec = delta_x_rec
    delta_nx_rec = delta_x_rec / dx
    delta_ny_rec = delta_y_rec / dy
    nx_rec = int(fix((nx - 2*npmlx) * dx / delta_x_rec)) - 1
    ny_rec = int(fix((ny - 2*npmly) * dy / delta_y_rec)) - 1
#    ny_rec = nx_rec
    dumx = npmlx
    dumy = npmly
    for i in range(nx_rec):
        dumx += delta_nx_rec
        dumy = npmly
        for i in range(ny_rec):
            dumy += delta_ny_rec
            nx_rec_real = dumx * dx
            ny_rec_real = dumy * dy
            nz_rec_real = (20 - 2) * dz
            nz_rec = 20 - 2
            # if nz_rec_real + z_axis[0] > 0:
                # mlab.points3d(nx_rec_real, ny_rec_real, nz_rec_real + z_axis[0],scale_factor = 0.05)
            rec.append([dumx + nshift, dumy + nshift, nz_rec, 'Ey'])

>>>>>>> gly:model_em_gly.py
    nrec = len(rec)
    logger.info("nrec: %d"%nrec) 

#    nrec = len(rec)
#    print "nrec: ",nrec
    #%%write data
    for i in range(nsrc):
        # Modified by mbw at 20180607
        # with open('./Input/src.in_' + str(i + 1).zfill(4), 'w') as fsrc:
        with open('./Input/src.in_' + str(i).zfill(4), 'w') as fsrc:
        # END Modify
            fsrc.write("%d %d\n" % (1, nt_src))
            # for i in range(nsrc):
            fsrc.write("%d %d %d %s\n" %
                       (src[i][0], src[i][1], src[i][2], src[i][3]))
            # for i in range(nsrc):
            savetxt(fsrc, array([srcpulse]) * src[i][4])

    with open('./Input/rec.in', 'w') as frec:
        frec.write("%d\n" % (nrec))
        for i in range(nrec):
            frec.write("%d %d %d %s\n" %
                       (rec[i][0], rec[i][1], rec[i][2], rec[i][3]))

<<<<<<< HEAD:model_em.py
    cp(os.path.join(indir,'src.in*'),std_indir)
    cp(os.path.join(indir,'rec.in'),os.path.join(std_indir,'rec.in'))
    if is_zRTM==1 or is_zRTM==2:
        with open(os.path.join(rtm0_indir,'rec.in'),'w+') as fo:
            fo.write('1\n')
            fo.write('%d %d %d Ey\n'%(nx0,ny0,nz_air-2))
    if is_zRTM==0 or is_zRTM==2:
        with open(os.path.join(rtm_indir,'rec.in'),'w+') as fo:
            fo.write('1\n')
            fo.write('%d %d %d Ey\n'%(nx0,ny0,nz_air-2))
            
    return nsrc,nrec
=======
    if Dum_RTM == 1:
      # os.system("cp ./src.in ./RTM/src.in")
        os.system("cp ./Input/rec.in ./RTM/Input/rec.in")
        os.system("cp ./Input/src.in* ./STD/Input/")
        os.system("cp ./Input/rec.in ./STD/Input/rec.in")
    return 0
>>>>>>> gly:model_em_gly.py

def eps_sig_mu(meps=1,meps_bg=False,msig=1e-5,msig_bg=False,mmiu=1,mmiu_bg=False):
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

<<<<<<< HEAD:model_em.py
    if is_zRTM==1 or is_zRTM==2:
        cp(os.path.join(std_indir, 'eps.in*'), rtm0_indir)
        cp(os.path.join(std_indir, 'sig.in*'), rtm0_indir)
        cp(os.path.join(std_indir, 'mu.in*'), rtm0_indir)
    if is_zRTM==0 or is_zRTM==2:
        cp(os.path.join(std_indir, 'eps.in*'), rtm_indir)
        cp(os.path.join(std_indir, 'sig.in*'), rtm_indir)
        cp(os.path.join(std_indir, 'mu.in*'), rtm_indir)
=======
    # cp(os.path.join(std_indir, 'eps.in*'), rtm_indir)
    # cp(os.path.join(std_indir, 'sig.in*'), rtm_indir)
    # cp(os.path.join(std_indir, 'mu.in*'), rtm_indir)
>>>>>>> gly:model_em_gly.py

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
    n = int(round(2 ** np.ceil(np.log2(len(srcpulse)))))
    freqs = np.linspace(0, int(1 / dt / 2), int(n / 2) + 1)
    sp = np.fft.rfft(srcpulse, n) / n
    W = abs(sp)
    fmax2 = max(freqs[W > max(W) / 10.0])
<<<<<<< HEAD:model_em.py
    logger.info("Src's main frequency: %fMHz" % (freqs[np.argmax(W)]/1e6))
    logger.info("!!check dx again (src_fmax(within 90%% of max amplitude)=%fMHz):"%(fmax2/1e6)) 
=======
    logger.info("Src's max frequency: %fMHz" % (freqs[np.argmax(W)]/1e6))
    # logger.info("!!check dx again (src_fmax(within 90%% of max amplitude)=%fMHz):"%(fmax2/1e6)) 
>>>>>>> gly:model_em_gly.py
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

<<<<<<< HEAD:model_em.py
    cp(os.path.join(indir, 'par.in'), std_indir)
    if is_zRTM==1 or is_zRTM==2:
        content[1] = "%e %e %e %e\n" % (dx, dy, dz, dt/2)
        with open(os.path.join(rtm0_indir, 'par.in'), 'w') as fpar:
            fpar.write(''.join(content))
    if is_zRTM==0 or is_zRTM==2:
        cp(os.path.join(indir, 'par.in'), rtm_indir)
=======
    # cp(os.path.join(indir, 'par.in'), std_indir)
    # if is_zRTM:
    #     content[1] = "%e %e %e %e\n" % (dx, dy, dz, dt/2)
    #     with open(os.path.join(rtm_indir, 'par.in'), 'w') as fpar:
    #         fpar.write(''.join(content))
    # else:
    #     cp(os.path.join(indir, 'par.in'), rtm_indir)
>>>>>>> gly:model_em_gly.py


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
    
<<<<<<< HEAD:model_em.py
    cp(fn_slice, std_indir)
    if is_zRTM==1 or is_zRTM==2:
        cp(fn_slice, rtm0_indir)
    if is_zRTM==0 or is_zRTM==2:
        cp(fn_slice, rtm_indir)


def cleanfiles(paths):

    def cleanChoice(path):
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


=======
    # cp(fn_slice, std_indir)
    # cp(fn_slice, rtm_indir)
    # if Dum_RTM == 1:
    #     # updated by mbw at 20190415
    #     shutil.copy('./Input/slice.in', './STD/Input/')
    #     shutil.copy('./Input/slice.in', './RTM/Input/')
>>>>>>> gly:model_em_gly.py

##############################################################################
if __name__ == '__main__':
    try:
<<<<<<< HEAD:model_em.py
        opts, args = getopt.getopt(sys.argv[1:], "z:m:f:", ["zero-offset=","model=","freq=","dx_src=","dx_rec=","dy_rec=","no_gen_model"])
=======
        opts, args = getopt.getopt(sys.argv[1:], "d:", ["workdir=","no_gen_model"])
>>>>>>> gly:model_em_gly.py
    except getopt.GetoptError as err:
        # print help information and exit:
        logger.error(err)  # will print something like "option -a not recognized"
        # usage()
        sys.exit(2)

<<<<<<< HEAD:model_em.py
    is_zRTM = 0
    model = "model.mat"
    gen_model = True
    workdir = os.path.join('tasks','default')
    freq = 300  #MHz
    dx_src = 0.4
    dx_rec = 0.2
    for o, a in opts:
        if o in ('-z','--zero-offset'):
            is_zRTM = int(a)
            logger.info('Zero-offset Mode: %d'%is_zRTM)
            if is_zRTM == 1:
                dx_src = 0.2
                dy_src = 0.5
        elif o in ('-m','--model'):
            model = a
        elif o in ('-f','--freq'):
            freq = int(a)  #MHz
        elif o in ('--dx_src'):
            dx_src = float(a)
        elif o in ('--dy_src'):
            dy_src = float(a)
        elif o in ('--dx_rec'):
            dx_rec = float(a)
=======
    gen_model = True
    workdir = os.path.join('tasks','default')
    for o, a in opts:
        if o in ('-d','--workdir'):
            workdir = os.path.join('tasks',a)
>>>>>>> gly:model_em_gly.py
        elif o in ('--no_gen_model'):
            gen_model = False
        else:
            assert False, "unhandled option"
    logger.info('workdir="%s"'%(workdir))

    indir = os.path.join(workdir,'Input')
    outdir = os.path.join(workdir,'Output')
    dirlist = [workdir,indir,outdir]
    for idir in dirlist:
        if not os.path.exists(idir):
            os.mkdir(idir)
        else:
            cleanfiles([idir])

<<<<<<< HEAD:model_em.py
    ### load model ###
    try:
        dic_model = sio.loadmat(os.path.join('Model',model))
    except Exception as e:
        logger.error(e)
        raise e
    else:
        modelname = str(dic_model['modelname'][0])
    finally:
        pass
    ### load model end ###
=======
>>>>>>> gly:model_em_gly.py

    # added by mbw at 20190415
    shutil.copy("FDTD_MPI_geop",workdir)

    #####################################################################
    mu0 = 1.2566370614e-6
    ep0 = 8.8541878176e-12

    epmin = 1.0 * 4
    mumin = 1.0
    epmax = 8.0 * 4
    mumax = 1.0
<<<<<<< HEAD:model_em.py
    fmax = freq * 1e6 #Hz
    ### parameter end ###

    ### init workdir ###
    if is_zRTM==1:
        dirname = '%s_%dMHz_0o_%.1fm_%.1fm'%(modelname,freq,dx_src,dy_src)
    else:
        dirname = '%s_%dMHz_%.1fm_%.1fm'%(modelname,freq,dx_src,dx_rec)
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
    if is_zRTM == 0 or is_zRTM == 2:
        dirlist1 += [rtm_dir]
        dirlist2 += [rtm_indir]
        dirlist3 += [rtm_outdir]
    if is_zRTM == 1 or is_zRTM == 2:
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
    cleanfiles(cleanlist)

    shutil.copy("FDTD_MPI",workdir)
    shutil.copy("FDTD_MPI",std_dir)
    if is_zRTM == 0 or is_zRTM == 2:
        shutil.copy("FDTD_MPI",rtm_dir)
    if is_zRTM == 1 or is_zRTM == 2:
        shutil.copy("FDTD_MPI",rtm0_dir)
    ### directories end ###


    ### gird  parameter ###
    npmlx = dic_model['npmlx']
    npmly = dic_model['npmly']
    npmlz = dic_model['npmlz']
=======
    fmax = 800e6  # Hz
>>>>>>> gly:model_em_gly.py

    outstep_t_wavefield = dic_model['outstep_t_wavefield']
    outstep_x_wavefield = dic_model['outstep_x_wavefield']
    outstep_slice = dic_model['outstep_slice']

    dx_max = finddx(epmax, mumax, fmax)
<<<<<<< HEAD:model_em.py
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
    if is_zRTM == 1:
        dnx_src = round(dx_src/dx)
        dny_src = round(dy_src/dy)
        [nsrc,nrec] = src_rec(dnx_src,dny_src, marginx=round(0.4/dx), marginy=round(0.4/dy))
    else:
        dnx_src = round(dx_src/dx)
        dnx_rec = round(dx_rec/dx)
        [nsrc,nrec] = src_rec(dnx_src,dnx_rec=dnx_rec, marginx=round(0.4/dx), marginx_rec=round(0.2/dy))
    ### generate src & rec end ###

    ### generate model ###
    NUM_OF_PROCESS = 4
    if gen_model:
        order = 2 # num of interchange layers of each process
        logger.info("NUM_OF_PROCESS: %d"%NUM_OF_PROCESS) 
=======
    dx = 0.05
    dy = 0.05
    dz = 0.05
    logger.info("dx=%f, dy=%f, dz=%f"%(dx, dy, dz)) 
    assert np.max([dx,dy,dz]) < dx_max, 'dx,dy,dz too big!!!'

    # y = arange(yaxis2[0], yaxis1[1], dx)
    # y_axis = y
    # x_axis = arange(0, 5, dx)
    # z_axis = arange(-1, 2.5, dx)

    nx = round(61/dx) # 60m
    ny = round(11//dx)  # 11m
    nz_air = 12
    nz = round(6//dx)+ nz_air# 6m
    logger.info('nx=%d, ny=%d, nz=%d'%(nx, ny, nz)) 

    npmlx = 8
    npmly = 8
    npmlz = 8
    nt_src = 50 # useless

    dt_max = finddt(epmin, mumin, dx, dy, dz)
    dt = 0.0587e-9
    nt = 952 + 5  # a little more steps, just in case.
    logger.info("dt=%g, nt=%d"%(dt, nt)) 
    assert dt < dt_max, 'dt too big!!! (%f>%f)'%(dt,dt_max)

    outstep_t_wavefield = 4
    outstep_x_wavefield = 2
    outstep_slice = 4


    NUM_OF_PROCESS = 16
    order = 2 # num of interchange layers of each process
    logger.info("NUM_OF_PROCESS: %d"%(NUM_OF_PROCESS)) 
    if gen_model:
>>>>>>> gly:model_em_gly.py
        rangex, nxSize = X_partition(nx, NUM_OF_PROCESS)
        eps_sig_mu(epmax)

    shutil.copy(os.path.join('src_rec','rec.in'), os.path.join(workdir,'Input'))
    shutil.copy(os.path.join('src_rec','8src_400MHz.in_0000'), os.path.join(workdir,'Input','src.in_0000'))
    par()
    islice(round(nx/2),round(ny/2),round(1/dz)+nz_air)

<<<<<<< HEAD:model_em.py
    # backup this file
    cp('model_em.py', workdir)
    cp(os.path.join('Model',model), workdir) # backup model
    cp(os.path.join('Model','make_model.m'), workdir) # backup make_model

    subtxt = 'python subgeop.py -d %s -s %d -p %d -z %d'%(dirname,nsrc,NUM_OF_PROCESS,is_zRTM)
    os.system(subtxt)
=======
    cp('model_em.py', workdir)

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

~/software/openmpi-4.0.1/bin/mpiexec -np $NSLOTS -wdir $WORKPATH $WORKPATH/FDTD_MPI_geop 0

echo "Computing is stopped at $(date)."

exit 0
''')
>>>>>>> gly:model_em_gly.py