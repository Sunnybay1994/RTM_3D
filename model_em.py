#!/usr/bin/env python
from numpy import *
import numpy as np
#import mayavi.mlab as mlab
from scipy.linalg import *
import matplotlib.pyplot as plt
import os,sys,shutil,logging,getopt
import glob

#logger
logger = logging.getLogger('model_em')
logger.setLevel(logging.INFO) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
fh = logging.FileHandler(os.path.join('log','model_em.log'))
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

def src_rec(dnx_src,dny_src=False,nx_src=False,ny_src=False,nzp_src=False,dnx_rec=False,dny_rec=False,nx_rec=False,ny_rec=False,nzp_rec=False,nshift=0):
    logger.info('adding source and receiver...')
    if not dny_src:
        dny_src = dnx_src
    if not nzp_src:
        nzp_src = nz_air - 1
    if not dnx_rec:
        dnx_rec = dnx_src
        if not dny_rec:
            dny_rec = dny_src
    else:
        if not dny_rec:
            dny_rec = dnx_rec
    if not nzp_rec:
        nzp_rec = nzp_src

    global src, rec
    def ant_pos(dnx_ant,dny_ant,nz_pos,nx_ant=False,ny_ant=False,nshift=0):
        ant = []
        dnx_ant = int(dnx_ant)
        dny_ant = int(dny_ant)
        nz_pos = int(nz_pos)

        if not nx_ant:
            nx_ant = (nx - 2*npmlx) // dnx_ant
            if nx_ant % 2 == 0:
                nx_ant -= 1
        if not ny_ant:
            ny_ant = (ny - 2*npmly) // dny_ant
            if ny_ant % 2 == 0:
                ny_ant -= 1
        
        for i in range(-(nx_ant-1)//2,(nx_ant+1)//2):
            dumx = nx0 + i*dnx_ant
            for j in range(-(ny_ant-1)//2,(ny_ant+1)//2):
                dumy = ny0 + j*dny_ant
                ant.append([dumx + nshift, dumy + nshift, nz_pos, 'Ey', 1])
        return ant

    src = ant_pos(dnx_src,dny_src,nzp_src,nx_src,ny_src,nshift)
    nsrc = len(src)
    logger.info("nsrc: %d"%nsrc) 

    rec = ant_pos(dnx_rec,dny_rec,nzp_rec,nx_rec,ny_rec,nshift)
    nrec = len(rec)
    logger.info("nrec: %d"%nrec) 

#    nrec = len(rec)
#    print "nrec: ",nrec
    # write data
    for i in range(nsrc):
        # Modified by mbw at 20180607
        # with open('./Input/src.in_' + str(i + 1).zfill(4), 'w') as fsrc:
        with open(os.path.join(indir,'src.in_' + str(i).zfill(4)), 'w') as fsrc:
        # END Modify
            fsrc.write("%d %d\n" % (1, nt_src))
            # for i in range(nsrc):
            fsrc.write("%d %d %d %s\n" %
                       (src[i][0], src[i][1], src[i][2], src[i][3]))
            # for i in range(nsrc):
            savetxt(fsrc, array([srcpulse]) * src[i][4])

    with open(os.path.join(indir,'rec.in'), 'w') as frec:
        frec.write("%d\n" % (nrec))
        for i in range(nrec):
            frec.write("%d %d %d %s\n" %
                       (rec[i][0], rec[i][1], rec[i][2], rec[i][3]))

    cp(os.path.join(indir,'src.in*'),std_indir)
    cp(os.path.join(indir,'rec.in'),os.path.join(std_indir,'rec.in'))
    # if is_zRTM:
    with open(os.path.join(rtm_indir,'rec.in'),'w+') as fo:
        fo.write('1\n')
        fo.write('%d %d %d Ey\n'%(nx0,ny0,nz_air-2))
    # else:
        # os.system("cp ./Input/rec.in ./RTM/Input/rec.in")
    return 0


def eps_sig_mu(model=False,modelname=False,depths=False,half_thickness=5,atype='ep'):
    logger.info('Generating model...')

    if model=='pic':
        model_em = plt.imread(modelname)[:, :, 0]
        model_em = model_em.T

    for ii in range(NUM_OF_PROCESS):

        fep_STD = open(os.path.join(std_indir,'eps.in_' + str(ii).zfill(3)), 'w')
        fsig_STD = open(os.path.join(std_indir,'sig.in_' + str(ii).zfill(3)), 'w')

        feps = open(os.path.join(indir,'eps.in_' + str(ii).zfill(3)), 'w')
        fsig = open(os.path.join(indir,'sig.in_' + str(ii).zfill(3)), 'w')

        logger.info('file: %s, %s'%(feps,fsig))
        logger.info('rangex: %d~%d,%d'%(rangex[ii], rangex[ii + 1], rangex[ii + 1] - rangex[ii]))
        for dumx in arange(rangex[ii] - order, rangex[ii + 1] + order):

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

            eps = ones((ny, nz))
            eps_STD = ones((ny, nz))
            sig = ones((ny, nz))
            sig_STD = ones((ny, nz))


            for j in range(ny):
                eps[j, 0:nz_air] = 1.0
                eps[j, nz_air:] = 9.0

                sig[j, 0:nz_air] = 1e-5
                sig[j, nz_air:] = 1e-5

                eps_STD[j, 0:nz_air] = 1.0
                eps_STD[j, nz_air:] = 9.0

                sig_STD[j, 0:nz_air] = 1e-5
                sig_STD[j, nz_air:] = 1e-5

                if model=='pic':
                    if 'ep' in atype:
                        for depth in depths:
                            bottom = depth - half_thickness
                            top =depth + half_thickness
                            eps[model_em[dumx,:] != 1, bottom:top] = 15

                    if 'sig' in atype:
                        for depth in depths:
                            bottom = depth - half_thickness
                            top =depth + half_thickness
                            sig[model_em[dumx,:] != 1, bottom:top] = 10
                else:
                    i = dumx
                    depth = depths[0]
                    if model=='dots':
                        ## 6 dots model
                        r = int(0.1/dx)
                        nx1 = nx0 - int(0.64/dx)
                        nx2 = nx0
                        nx3 = nx0 + int(0.64/dx)
                        ny1 = ny0 - int(dy_src/dy)
                        ny2 = ny0
                        ny3 = ny0 + int(dy_src/dy)
                        npdh = int(0.16/dx) 
                        for k in range(nz_air+5,nz):
                            if distance(i,j,k,nx1-npdh,ny1,depth) < r or distance(i,j,k,nx2-npdh,ny2,depth) < r or distance(i,j,k,nx2-npdh,ny2,depth) < r:
                                eps[j,k] = 3
                                logger.debug('ep(%.2f,%.2f,%.2f)=3'%((i-nx0)*dx,(j-ny0)*dy,(k-nz_air)*dz))
                            elif distance(i,j,k,nx1+npdh,ny1,depth) < r or distance(i,j,k,nx2+npdh,ny2,depth) < r or distance(i,j,k,nx2+npdh,ny2,depth) < r:
                                eps[j,k] = 15
                                logger.debug('ep(%.2f,%.2f,%.2f)=15'%((i-nx0)*dx,(j-ny0)*dy,(k-nz_air)*dz))
                    elif model == 'surface':
                        ## model: x/6 + y/6 + z/2 = 1
                        for k in range(nz_air+5,nz):
                            if abs(i/300 + j/300 + (k-nz_air)/100 - 1) < 0.01:
                                logger.debug('eps(%d,%d,%d)=15, %f'%(i,j,k,i/300 + j/300 + (k-nz_air)/100))
                                eps[j,k] = 15
                    elif model == 'checkbox':
                        ## checkbox model
                        nlbox = 30 #0.6m
                        ep1 = 6
                        ep2 = 12
                        for k in range(nz_air,nz):
                            if (i // nlbox + j // nlbox + (k - nz_air) // nlbox) % 2 == 0:
                                eps[j,k] = ep1
                                # logger.debug('eps(%d,%d,%d)=%d'%(i,j,k,ep1))
                            else:
                                eps[j,k] = ep2
                                # logger.debug('eps(%d,%d,%d)=%d'%(i,j,k,ep2))



            # if dumx == 0:
            #     plt.figure()
            #     plt.imshow(fliplr(eps).T)
            #     plt.show()
            # savetxt(feps, fliplr(eps))
            # savetxt(fep_STD, fliplr(eps_STD))
            savetxt(feps, eps, fmt='%.4e')
            savetxt(fep_STD, eps_STD, fmt='%.4e')
            savetxt(fsig, sig, fmt='%.4e')
            savetxt(fsig_STD, sig_STD, fmt='%.4e')
        feps.close()
        fep_STD.close()
        fsig.close()
        fsig_STD.close()

        # fsig = open('./Input/sig.in_' + str(ii).zfill(3), 'w')
        # dum_sig = zeros((nxSize[ii], ny, nz))
        # sig_out = dum_sig.reshape(nxSize[ii] * ny * nz)
        # savetxt(fsig, sig_out)
        # fsig.close()

        fmu = open(os.path.join(indir,'mu.in_' + str(ii).zfill(3)), 'w')
        # Modified by mbw at 20180607
        # dum_mu = ones((nxSize[ii], ny, nz))
        dum_mu = ones((int(nxSize[ii]), ny, nz))
        # mu_out = dum_mu.reshape(nxSize[ii] * ny * nz)
        mu_out = dum_mu.reshape(int(nxSize[ii]) * ny * nz)
        # End modify
        savetxt(fmu, mu_out, fmt='%.4e')
        fmu.close()

    
    cp(os.path.join(std_indir, 'eps.in*'), rtm_indir)
    cp(os.path.join(std_indir, 'sig.in*'), rtm_indir)

    cp(os.path.join(indir, 'mu.in*'), std_indir)
    cp(os.path.join(indir, 'mu.in*'), rtm_indir)

    return 0


##############################################################################
def finddt(epmin, mumin, dx, dy, dz):
    epmin = epmin * ep0
    mumin = mumin * mu0
    dtmax = 6.0 / 7.0 * \
        sqrt(epmin * mumin / (1.0 / dx ** 2 + 1.0 / dy ** 2 + 1.0 / dz ** 2))
    logger.info("dt max = %fns"%(dtmax/1e-9)) 
    return dtmax


def finddx(epmax, mumax, fmax):
    epmax = epmax * ep0
    mumax = mumax * mu0
    wlmin = 1 / (fmax * sqrt(epmax * mumax))
    dxmax = wlmin
    logger.info("max dx = %f"%dxmax) 
    return dxmax


def blackharrispulse(fmax, dt):
    a = [0.35322222, -0.488, 0.145, -0.010222222]
    T = 1.14 / fmax
    t = arange(0, T, dt)
    window = zeros(size(t))
    for n in range(4):
        window = window + a[n] * cos(2 * n * pi * t / T)

    window[t >= T] = 0
    p = window
    p = window[:] - append(0, window[:-1])
    p = p / max(abs(p))
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
    plt.figure()
    plt.plot(t, y)
    plt.savefig(os.path.join(workdir,'ricker.png'))
    return y


def check_dx(srcpulse):
    n = 2 ** int(ceil(log2(len(srcpulse))))
    freqs = np.linspace(0, 1 / dt / 2, n / 2 + 1)
    sp = np.fft.rfft(srcpulse, n) / n
    W = abs(sp)
    fmax2 = max(freqs[W > max(W) / 10.0])
    logger.info("!!check dx again:") 
    finddx(epmax, mumax, fmax2)
    logger.info("Src's max frequency: %fMHz" % (freqs[argmax(W)]/1e6))
    plt.figure()
    plt.plot(freqs, W)
    plt.title('frequency spectrum of source')
    plt.savefig(os.path.join(workdir,'spectral_src.png'))


def distance(x, y, z, x0, y0, z0):
    return sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2)


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

    cp(os.path.join(indir, 'par.in'), std_indir)
    if is_zRTM:
        content[1] = "%e %e %e %e\n" % (dx, dy, dz, dt/2)
        with open(os.path.join(rtm_indir, 'par.in'), 'w') as fpar:
            fpar.write(''.join(content))
    else:
        cp(os.path.join(indir, 'par.in'), rtm_indir)


def X_partition(nx, NUM_OF_PROCESS):
    # global
    rangex = zeros(NUM_OF_PROCESS + 1)
    nxSize = zeros(NUM_OF_PROCESS)
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
    
    cp(fn_slice, std_indir)
    cp(fn_slice, rtm_indir)


def cleanfiles(paths):

    def cleanChoice(path):
        choice = input('Clear %s? A/Y/(N)'%path).lower()
        if choice in ['all', 'a']:
            confirm = input("Clear all files in '%s'? Y/(N)"%"''&'".join(paths+[path])).lower()
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



##############################################################################
if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "zm:p:d:", ["zero-offset","model=","modelpic=","workdir="])
    except getopt.GetoptError as err:
        # print help information and exit:
        logger.error(err)  # will print something like "option -a not recognized"
        # usage()
        sys.exit(2)

    is_zRTM = False
    model = False # os.path.join('Model','EM300.png')
    modelpic = False
    workdir = os.path.join('tasks','default')
    for o, a in opts:
        if o in ('-z','--zero-offset'):
            logger.info('Zero-offset Mode.')
            is_zRTM = True
        elif o in ('-m','--model'):
            model = a
            if model == '0':
                logger.info("Don't generate model.")
        elif o in ('-p','--modelpic'):
            modelpic = os.path.join('Model',a)
        elif o in ('-d','--workdir'):
            workdir = os.path.join('tasks',a)
        else:
            assert False, "unhandled option"

    if not os.path.exists(workdir):
        os.mkdir(workdir)
    logger.info('model="%s", workdir="%s"'%(model,workdir))


    ### parameter ###
    mu0 = 1.2566370614e-6
    ep0 = 8.8541878176e-12

    epmin = 1.0
    mumin = 1.0
    epmax = 15.0
    mumax = 1.0
    fmax = 300*1e6  #Hz
    ### parameter ###

    ### gird  parameter ###
    npmlx = 9
    npmly = npmlx
    npmlz = npmlx

    dx_max = finddx(epmax, mumax, fmax)
    dx = 0.1
    dy = dx
    dz = dx
    logger.info("dx=%f, dy=%f, dz=%f"%(dx, dy, dz)) 

    nx = 100
    ny = 100
    nx0 = nx/2
    ny0 = ny/2
    nz_air = 10
    nz0 = nz_air 
    nz = nz0 + 60
    logger.info('nx=%d, ny=%d, nz=%d'%(nx, ny, nz))

    dt_max = finddt(epmin, mumin, dx, dy, dz)
    dt = 1e-10
    nt = 1000 #100ns
    logger.info("dt: %fns"%(dt/1e-9)) 
    ### gird paraeter ###

    ### source ###
    # srcpulse = blackharrispulse(fmax, dt)
    srcpulse = ricker(fmax,4*1/fmax,dt)
    nt_src = len(srcpulse)
    logger.info("nt=%d, nt_src=%d"%(nt, nt_src)) 
    check_dx(srcpulse)
    ### source ###

    outstep_t_wavefield = 5
    outstep_x_wavefield = 2
    outstep_slice = 5

    ### directories ###
    indir = os.path.join(workdir,'Input')
    outdir = os.path.join(workdir,'Output')
    std_dir = os.path.join(workdir,'STD')
    std_indir = os.path.join(std_dir,'Input')
    std_outdir = os.path.join(std_dir,'Output')
    rtm_dir = os.path.join(workdir,'RTM')
    rtm_indir = os.path.join(rtm_dir,'Input')
    rtm_outdir = os.path.join(rtm_dir,'Output')
    dirlist = [indir,outdir,std_dir,std_indir,std_outdir,rtm_dir,rtm_indir,rtm_outdir]

    for idir in dirlist:
        if not os.path.exists(idir):
            os.mkdir(idir)

    cleanfiles([indir,std_indir,rtm_indir])
    shutil.copy("FDTD_MPI_geop",workdir)
    shutil.copy("FDTD_MPI_geop",std_dir)
    shutil.copy("FDTD_MPI_geop",rtm_dir)
    ### directories ###

    obstacle_depth = [1]
    ob_nz = [int(nz_air+z/dz) for z in obstacle_depth]
    logger.info('obstacle_depth: %s (nz=%s)'%(obstacle_depth,ob_nz))

    if is_zRTM:
        dx_src = 0.2
        dy_src = 0.5
        dnx_src = int(dx_src/dx)
        dny_src = int(dy_src/dy)
        nx_src = 25
        ny_src = 7
        src_rec(dnx_src,dny_src,nx_src=nx_src,ny_src=ny_src,nx_rec=nx_src,ny_rec=ny_src)
    else:
        dx_src = 0.6
        dx_rec = 0.3
        dnx_src = int(dx_src/dx)
        dnx_rec = int(dx_rec/dx)
        src_rec(dnx_src=dnx_src, dnx_rec=dnx_rec)

    if not model == '0':
        NUM_OF_PROCESS = 8
        order = 2 # num of interchange layers of each process
        logger.info("NUM_OF_PROCESS: %d"%NUM_OF_PROCESS) 
        rangex, nxSize = X_partition(nx, NUM_OF_PROCESS)
        eps_sig_mu(model,modelname=modelpic,depths=ob_nz)


    par()
    islice(nx0-1,ny0-1,ob_nz)

    # clean file
    cleanfiles([outdir,std_outdir,rtm_outdir])

    # backup this file
    cp('model_em.py', workdir)