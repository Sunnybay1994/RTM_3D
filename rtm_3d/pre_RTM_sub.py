#!/usr/bin/env python
import numpy as np
import matplotlib
import scipy.io as sio
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os,sys,logging,argparse,struct
from common import *

def merge_gather(wdir, isrc):
    # global isum
    # global iloc
    idir = os.path.join(wdir,'Output')
    fn_data = os.path.join(idir, 'merge_gather_'+str(isrc).zfill(4)+'.bin')
    fn_data_ext = os.path.join(idir, 'ext','merge_gather_'+str(isrc).zfill(4)+'.bin')
    
    try:
        logger.info("merging gather: Try loading external merged gathers first.")
        with open(os.path.join(idir, 'merge_gather_loc_'+str(isrc).zfill(4)+'.dat')) as fp:
            iloc = np.loadtxt(fp)
        nrec = np.array(iloc).shape[0]
        with open(fn_data_ext,'rb') as fp:
            isum_ext = np.fromfile(fp,dtype='float32')
            isum_ext = isum_ext.reshape(nrec,-1)
        return isum_ext,iloc
    except Exception as e:
        logger.info('merging external gather WRONG: %s'%e)

    if 'win' in sys.platform:
        ilist = os.popen('dir/b/on '+ os.path.join(idir,'gather'+'_'+str(isrc).zfill(4)+'*')).readlines()
        logger.debug(ilist)
        ilist = [os.path.join(idir,f) for f in ilist]
        logger.debug(ilist)
    elif 'linux' in sys.platform:
        ilist = os.popen('ls '+ os.path.join(idir,'gather'+'_'+str(isrc).zfill(4)+'*')).readlines()
    else:
        logger.error('unknown platform: %s'%sys.platform)
        return 0
    isum = []
    iloc = []
    if ilist:
        if forward_method == 'fdtd':
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
        elif forward_method == 'pstd':
            logger.info("Reading gather from bin (%d files): %s-src%d"%(len(ilist),idir,isrc))
            with open(os.path.join(wdir,'Input','rec.in'),'r') as fo:
                nrec = int(fo.readline())
                for line in fo.readlines():
                    iloc.append([int(x) for x in line.split(',')[:3]])
                iloc = np.array(iloc)
            with open(ilist[0].strip('\n'),'rb') as fo:
                data_raw = struct.unpack('f'*nt*nrec,fo.read())
                isum = np.reshape(data_raw,(nrec,nt))
        plt.figure()
        for i in range(len(isum)):
            plt.plot(isum[i]/max(abs(isum[i]))+i)
        plt.savefig(os.path.join(idir, 'merge_gather_'+str(isrc).zfill(4)+'.png'))
        np.array(isum).astype('float32').tofile(fn_data)
        with open(os.path.join(idir, 'merge_gather_loc_'+str(isrc).zfill(4)+'.dat'),'w') as fp:
            np.savetxt(fp,iloc,'%d')
        for fname in ilist:
            os.remove(fname.strip('\n'))
    else:
        logger.info("merging gather: No files. Loading merged gathers.") 
        with open(os.path.join(idir, 'merge_gather_loc_'+str(isrc).zfill(4)+'.dat')) as fp:
            iloc = np.loadtxt(fp)
        nrec = iloc.shape[0]
        isum = np.fromfile(fn_data,dtype='float32')
        isum = isum.reshape(nrec,-1)
    return isum,iloc


def remove_STD(gather,gather_std,isrc,remove_std=True):
    # global gather_RTM
    logger.info("removing source: src%d"%isrc) 
    
    gather_RTM = np.fliplr(np.array(gather))
    nt = gather_RTM.shape[1]
    plt.figure()
    plt.imshow(gather_RTM,cmap='gray',origin='lower',extent=(0,nt,0,nt/2))
    plt.savefig(os.path.join(rtmdir,'Input','gather'+'_'+str(isrc).zfill(4)+'.png'))
    if remove_std:
        gather_RTM = np.fliplr(np.array(gather) - np.array(gather_std))
        nt = gather_RTM.shape[1]
        plt.figure()
        plt.imshow(gather_RTM,cmap='gray',origin='lower',extent=(0,nt,0,nt/2))
        plt.savefig(os.path.join(rtmdir,'Input','gather_without_src'+'_'+str(isrc).zfill(4)+'.png'))
    # for i in range(len(gather_RTM[:,0])):
    #     gather_RTM[i,:] = gather_RTM[i,:]/max(abs(gather_RTM[i,:]))
    return gather_RTM


def prepare_RTM(isum,iloc,isum_std,iloc_std,gather_rtm,isrc):
    logger.info("preparing RTM: src%d"%isrc) 
    nrec = len(iloc)
    nrec_std = len(iloc_std)
    if nrec != nrec_std or nrec0!=nrec:
        logger.warning("nrec not equal! %d(%d)-%d"%(nrec,nrec0,nrec_std)) 
    nt = len(isum[0])
    nt_std = len(isum_std[0])
    if int(nt) != nt_std:
        logger.warning("nt not equal! %d-%d"%(nt,nt_std) )

    if 'm' in mode:
        fn = os.path.join(rtmdir,'Input','src.in'+'_'+str(isrc).zfill(4))
        srcinfos = [item.tolist()+[component] for item in iloc_std]
        extend_and_write_sources(fn,srcinfos,gather_rtm)

    if 'z' in mode:
        
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
    global nrec0,component
    with open("./Input/rec.in","r") as fp:
        nrec0 = int(fp.readline().split(',')[0])
        component = fp.readline().strip().split(',')[3]
        # print(fp.readline().split())
    if 'z' in mode:
        global dt,dx,dy,dz,v,z,max_offset,srclocs
        locs = []
        dats = []
        offsets = []
        logger.info('loading para')
        dic_model = sio.loadmat(os.path.join('model.mat'))
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
        try:
            srcinfo = np.round(dic_model['srcinfo']).astype('int')
            srclocs = [(srcinfo[0][i],srcinfo[1][i],srcinfo[2][i]) for i in range(np.size(srcinfo,1))]
            logger.info('Loading src locations from srcx,srcy,srxz.')
        except Exception as e:
            logger.info('%s:Loading src locations from srcx,srcy,srxz.'%e)
            srcxs = [int(round(x/dx)) for x in dic_model['srcx'][0].tolist()]
            srcys = [int(round(y/dy)) for y in dic_model['srcy'][0].tolist()]
            srcz = int(round(float(dic_model['srcz'])/dz + nz_air))
            srclocs = [locxy+(srcz,) for locxy in zip(srcxs,srcys)]

        dsrc_grid = np.linalg.norm(np.array(srclocs[1])-np.array(srclocs[0]))
        max_offset = np.sqrt(2)*dsrc_grid
        if no_nmo:
            max_offset = 1
        logger.debug('max_offset:%d'%max_offset)

    for i in list_src:#range(nsrc)
        logger.info('Prepare RTM for src%d.'%isrc)
        isum,iloc = merge_gather(workdir,i)
        isum_std,iloc_std = merge_gather(os.path.join(workdir,'STD'),i)
        gather_rtm = remove_STD(isum,isum_std,i, rm_std)
        
        if 'z' in mode and len(list_src)==nsrc:
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
                status_info = extend_and_write_sources(fn, locinfos, dats_binned, autospan=True)
                logger.info(status_info)
                #     fsrc.write("%d %d\n" % (len(list_src), nt))
                #     fsrc.write(locs)
                #     np.savetxt(fsrc,dats)
                logger.info('All Done.')
        if 'm' in mode:
            prepare_RTM(isum,iloc,isum_std,iloc_std,gather_rtm,i)

        logger.info('Done: src%d'%i)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Prepare data before RTM',conflict_handler='resolve')
    parser.add_argument('isrc',type=int,default=-1,help="The source No. to be processed in multi-offset mode; Number of sources in zero-offset mode.")
    parser.add_argument('-m','--mode',choices=['m','z'],default='m',help="Mode: 'm' for multi-offset, 'z' for zero-offset")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--forward_method',choices=['fdtd','pstd'],default='fdtd',help="Forward method used in RTM.")
    group.add_argument('--fdtd',action='store_const',const='fdtd',dest='forward_method',help='Use finite difference time domain as the forward method.')
    group.add_argument('--pstd',action='store_const',const='pstd',dest='forward_method',help='Use pseudo spectral time domain as the forward method.')
    parser.add_argument('--workdir',type=str,default='',help="The parent directory of 'OUTPUT' and 'STD'.")
    # parser.add_argument('--half_span',type=int,default=-1,help='Span the point source to a (1+2*half_span)x(1+2*half_span) source while forwarding. It helps to reduce the sideslobe in FDTD forwarding.')
    parser.add_argument('--nmo',action='store_const',const=False,dest='no_nmo',default=True,help="Use nmo to traansfer some CMP to zero_offset traces.")
    parser.add_argument('--no_rm_std',action='store_const',const=False,dest='rm_std',default=True,help="Do not remove std(eg. direct waves) from gathers.")
    
    args = parser.parse_args()
    isrc = args.isrc
    mode = args.mode
    forward_method = args.forward_method
    workdir = args.workdir
    no_nmo = args.no_nmo
    rm_std = args.rm_std

    # logger
    logger=addlogger(os.path.basename(sys.argv[0]))

    if 'm' in mode:
        rtmdir_name = 'RTM'
        rtmdir = os.path.join(workdir,rtmdir_name)
        logger.info('Multi-offset mode. src%d'%isrc)
        if isrc>=0:
            pre_RTM([isrc])

    if 'z' in mode:
        rtmdir_name = 'RTM0'
        rtmdir = os.path.join(workdir,rtmdir_name)
        statusdir = os.path.join(rtmdir,'status')
        nsrc = len([f for f in os.listdir('Input') if os.path.isfile(os.path.join('Input',f)) and 'src.in_' in f])
        logger.info('Zero-offset mode. src%d/%d'%(isrc,nsrc))
        if isrc == 0 and os.path.isfile(os.path.join(statusdir,'force_run')):
            logger.info('(Force run)Zero-offset mode. We have %d sources to process.'%nsrc)
            pre_RTM(range(nsrc))
            exit(100)
        else:
            pre_RTM([isrc])
            logger.info('Zero-offset src%d/%d done.'%(isrc,nsrc))
            ifn = os.path.join(statusdir,'src{isrc}_done')
            with open(ifn.format(isrc=isrc),'w') as fo:
                pass
            logger.info('Checking if all done.')
            for i in range(nsrc):
                if not os.path.isfile(ifn.format(isrc=i)):
                    logger.info('src%d not done: progress %d/%d'%(isrc,i+1,nsrc))
                    exit()
            print('All src done, prepareing rtm.')
            pre_RTM(range(nsrc))
            cleanfiles(statusdir,'y')
            exit(100)

