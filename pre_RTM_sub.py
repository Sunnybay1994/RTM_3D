#!/usr/bin/env python
import numpy as np
import matplotlib
import scipy.io as sio
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os,sys,logging,getopt
from par_RTM import *
from writesource import *
from normal_moveout import *

#logger
logger = logging.getLogger('pre_RTM_sub')
logger.setLevel(logging.DEBUG) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
fh = logging.FileHandler('log/pre_RTM_sub.log')
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s(%(process)d-%(processName)s): (%(levelname)s) %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)

def merge_gather(idir, isrc):
    # global isum
    # global iloc
    if 'win' in sys.platform:
        ilist = os.popen('dir/b/on '+ os.path.join(idir,'gather*dat'+'_'+str(isrc).zfill(4))).readlines()
        logger.debug(ilist)
        ilist = [os.path.join(idir,f) for f in ilist]
        logger.debug(ilist)
    elif 'linux' in sys.platform:
        ilist = os.popen('ls '+ os.path.join(idir,'gather*dat'+'_'+str(isrc).zfill(4))).readlines()
    else:
        logger.error('unknown platform: %s'%sys.platform)
        return 0
    isum = []
    iloc = []
    if ilist:
        logger.info("merging gather for %d files: %s-src%d"%(len(ilist),idir,isrc)) 
        for name in ilist:
            gather = np.loadtxt(name.strip('\n'))
            if np.shape(gather)[0] == 0:
                continue
            if np.shape(np.shape(gather))[0] == 1:
                iloc.append((int(x) for x in gather[:3]))
                isum.append(gather[3:])
            else:
                for i in range(len(gather[:,0])):
                    iloc.append(gather[i,:3])
                    isum.append(gather[i,3:])
        plt.figure()
        for i in range(len(isum)):
            plt.plot(isum[i]/max(abs(isum[i]))+i)
        plt.savefig(os.path.join(idir, 'merge_gather_'+str(isrc).zfill(4)+'.png'))
        with open(os.path.join(idir, 'merge_gather_'+str(isrc).zfill(4)+'.dat'),'w') as fp:
            np.savetxt(fp,isum)
        with open(os.path.join(idir, 'merge_gather_loc_'+str(isrc).zfill(4)+'.dat'),'w') as fp:
            np.savetxt(fp,iloc,'%d')
        for fname in ilist:
            os.remove(fname.strip('\n'))
    else:
        logger.info("merging gather: No files. Loading merge_gather.") 
        with open(os.path.join(idir, 'merge_gather_'+str(isrc).zfill(4)+'.dat')) as fp:
            isum = np.loadtxt(fp)
        with open(os.path.join(idir, 'merge_gather_loc_'+str(isrc).zfill(4)+'.dat')) as fp:
            iloc = np.loadtxt(fp)
    return isum,iloc


def remove_STD(gather,gather_std,isrc):
    # global gather_RTM
    logger.info("removing source: src%d"%isrc) 
    gather_RTM = np.fliplr(np.array(gather) - np.array(gather_std))
    # for i in range(len(gather_RTM[:,0])):
    #     gather_RTM[i,:] = gather_RTM[i,:]/max(abs(gather_RTM[i,:]))
    nt = gather_RTM.shape[1]
    plt.figure()
    plt.imshow(gather_RTM,cmap='gray',origin='lower',extent=(0,nt,0,nt/2))
    plt.savefig(os.path.join(rtmdir,'Input','gather_without_src'+'_'+str(isrc).zfill(4)+'.png'))
    return gather_RTM


def prepare_RTM(isum,iloc,isum_std,iloc_std,gather_rtm,isrc):
    logger.info("preparing RTM: src%d"%isrc) 
    with open("./Input/rec.in","r") as fp:
        nrec0 = int(fp.readline().split(',')[0])
        component = fp.readline().split(',')[3]
        # print(fp.readline().split())
    nrec = len(iloc)
    nrec_std = len(iloc_std)
    if nrec != nrec_std or nrec0!=nrec:
        logger.warning("nrec not equal! %d(%d)-%d"%(nrec,nrec0,nrec_std)) 
    nt = len(isum[0])
    nt_std = len(isum_std[0])
    if int(nt) != nt_std:
        logger.warning("nt not equal! %d-%d"%(nt,nt_std) )

    if mode == 1:
        fn = os.path.join(rtmdir,'Input','src.in'+'_'+str(isrc).zfill(4))
        srcinfos = [item.tolist()+[component] for item in iloc]
        extend_and_write_sources(fn,srcinfos,gather_rtm)

    if mode == 0:
        
        tt = []
        for ig in range(len(gather_rtm)):
            tt += [np.array(range(len(gather_rtm[ig])))*dt]

        logger.info('Doing NMO...')
        nmo_results = gen_nmo_gathers(srclocs[isrc], iloc, max_offset, gather_rtm, tt, v[:,:,z>0],z[z>0],dx,dy,dz)
        return nt,component,nmo_results

def bin_gathers(locs,dats,offsets=False,decay_fac=0.5,source_span=5):
    assert len(locs)==len(dats), 'number of gathers(%d) and its locations(%d) not equal.'%(len(dats),len(locs))
    ndat = len(locs)
    logger.info('bin gathers into grids...')
    locx_max,locy_max,locz = max(locs)
    ilocx_max = int(np.ceil(locx_max/source_span)) + 1
    ilocy_max = int(np.ceil(locy_max/source_span)) + 1
    bin_ga = np.zeros([ilocx_max,ilocy_max,nt])
    bin_ga_sum = np.zeros([ilocx_max,ilocy_max])
    # put into bin
    logger.debug('put into bin')
    for ig in range(ndat):
        loc = locs[ig]
        ga = dats[ig]
        if offsets:
            offset = offsets[ig]
            weight = np.exp(-decay_fac * offset)
        else:
            weight = 1
        ixbin = int(round(loc[0]/source_span))
        iybin = int(round(loc[1]/source_span))
        bin_ga[ixbin,iybin,:] += ga * weight
        bin_ga_sum[ixbin,iybin] += weight
    # merge gather in the same location
    logger.debug('merge gather in the same location')
    for ix in range(ilocx_max):
        for iy in range(ilocy_max):
            if bin_ga_sum[ix,iy] == 0:
                # logger.debug('no data at: (%d,%d)*%d'%(ix,iy,source_span))
                pass
            else:
                yield (ix*source_span,iy*source_span,locz),bin_ga[ix][iy]/bin_ga_sum[ix][iy]


def pre_RTM(list_src):
    if mode == 0:
        global dt,dx,dy,dz,v,z,max_offset,srclocs
        locs = []
        dats = []
        offsets = []
        logger.info('loading para')
        dic_model = sio.loadmat(os.path.join('model.mat'))
        dict_sr = sio.loadmat(os.path.join('model_sr.mat'))
        epr = np.array(dic_model['ep'])
        dt = float(dic_model['dt'])
        dx = float(dic_model['dx'])
        dy = float(dic_model['dy'])
        dz = float(dic_model['dz'])
        nx = int(dic_model['nx'])
        ny = int(dic_model['ny'])
        nz_air = int(dic_model['nz_air'])
        nz = int(dic_model['nz'])
        v = 299792458/np.sqrt(epr)
        z = (np.array(range(nz)) - nz_air) * dz
        srcxs = [int(round(x/dx)) for x in dict_sr['srcx'][0].tolist()]
        srcys = [int(round(y/dy)) for y in dict_sr['srcy'][0].tolist()]
        srcz = int(round(float(dict_sr['srcz'])/dz + nz_air))
        srclocs = [locxy+(srcz,) for locxy in zip(srcxs,srcys)]

        dsrc_grid = np.linalg.norm(np.array(srclocs[1])-np.array(srclocs[0]))
        max_offset = np.sqrt(2)*dsrc_grid
        if no_nmo:
            max_offset = 1
        logger.debug('max_offset:%d'%max_offset)

    for i in list_src:#range(nsrc)
        isum,iloc = merge_gather(os.path.join(workdir,'Output'),i)
        isum_std,iloc_std = merge_gather(os.path.join(workdir,'STD','Output'),i)
        gather_rtm = remove_STD(isum,isum_std,i)
        
        if mode == 0:
            nt,component,nmo_results = prepare_RTM(isum,iloc,isum_std,iloc_std,gather_rtm,i)
            nmo_locs,nmo_gathers,nmo_offsets = zip(*nmo_results)
            locs += nmo_locs
            dats += nmo_gathers
            offsets += nmo_offsets

            if i == list_src[-1]:
                locs_binned,dats_binned = zip(*bin_gathers(locs,dats,offsets))
                locinfos = [(loc[0],loc[1],loc[2],component) for loc in locs_binned]
                logger.info('Writing zero-offset source data...')
                fn = os.path.join(rtmdir,'Input','src.in_0000')
                extend_and_write_sources(fn, locinfos, dats_binned)
                #     fsrc.write("%d %d\n" % (len(list_src), nt))
                #     fsrc.write(locs)
                #     np.savetxt(fsrc,dats)
                logger.info('All Done.')
        else:
            prepare_RTM(isum,iloc,isum_std,iloc_std,gather_rtm,i)

        logger.info('Done: src%d\n'%i)


if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "m:o:d:", ["mode=","outdir=","workdir=","no_nmo"])
    except getopt.GetoptError as err:
        # print help information and exit:
        logger.error(err)  # will print something like "option -a not recognized"
        # usage()
        sys.exit(2)

    mode = 1
    workdir = ''
    rtmdir_name = 'RTM'
    no_nmo = False
    for o, a in opts:
        if o in ('-m','--mode'):
            if a == '0':
                mode = 0
                nsrc = len([f for f in os.listdir('Input') if os.path.isfile(os.path.join('Input',f)) and 'src.in_' in f])
                rtmdir_name = 'RTM0'
                logger.info('Zero-offset mode. We have %d sources to process.'%nsrc)
            elif a == '1':
                mode = 1
                logger.info('normal mode. src%d'%isrc)
            else:
                logger.error('mode parameter must be 0,1')
        elif o in ('-d','--workdir'):
            workdir = a
        elif o in ('-o','--outdir'):
            rtmdir_name = a
            logger.info('Output dir = "%s"'%rtmdir_name)
        
        elif o in ("--no_nmo"):
            no_nmo = True
        else:
            assert False, "unhandled option"

    rtmdir = os.path.join(workdir,rtmdir_name)
    read_par
    if mode == 0:
        pre_RTM(range(nsrc))
        # pre_RTM([60])
    elif mode == 1:
        isrc = int(args[0])
        pre_RTM([isrc])

