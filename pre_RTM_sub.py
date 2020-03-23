#!/usr/bin/env python
from numpy import *
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import *
import os,sys,logging,getopt
from par_RTM import *

#logger
logger = logging.getLogger('pre_RTM_sub')
logger.setLevel(logging.INFO) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
fh = logging.FileHandler('log/pre_RTM_sub.log')
# fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
# ch.setLevel(logging.DEBUG)
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
            gather = loadtxt(name.strip('\n'))
            if shape(gather)[0] == 0:
                continue
            if shape(shape(gather))[0] == 1:
                iloc.append(gather[:3])
                isum.append(gather[3:])
            else:
                for i in range(len(gather[:,0])):
                    iloc.append(gather[i,:3])
                    isum.append(gather[i,3:])
        figure()
        for i in range(len(isum)):
            plot(isum[i]/max(abs(isum[i]))+i)
        savefig(os.path.join(idir, 'merge_gather_'+str(isrc).zfill(4)+'.png'))
        with open(os.path.join(idir, 'merge_gather_'+str(isrc).zfill(4)+'.dat'),'w') as fp:
            savetxt(fp,isum)
        with open(os.path.join(idir, 'merge_gather_loc_'+str(isrc).zfill(4)+'.dat'),'w') as fp:
            savetxt(fp,iloc)
        for fname in ilist:
            os.remove(fname.strip('\n'))
    else:
        logger.info("merging gather: No files. Loading merge_gather.") 
        with open(os.path.join(idir, 'merge_gather_'+str(isrc).zfill(4)+'.dat'),'w') as fp:
            isum = loadtxt(fp)
        with open(os.path.join(idir, 'merge_gather_loc_'+str(isrc).zfill(4)+'.dat'),'w') as fp:
            iloc = loadtxt(fp)
    return isum,iloc


def remove_STD(gather,gather_std,isrc):
    # global gather_RTM
    logger.info("removing source: src%d"%isrc) 
    gather_RTM = fliplr(array(gather) - array(gather_std))
    # for i in range(len(gather_RTM[:,0])):
    #     gather_RTM[i,:] = gather_RTM[i,:]/max(abs(gather_RTM[i,:]))
    figure()
    imshow(gather_RTM,cmap='gray',origin='lower',extent=(0,nt,0,nt/2))
    savefig(os.path.join(rtmdir,'Input','gather_without_src'+'_'+str(isrc).zfill(4)+'.png'))
    return gather_RTM


def prepare_RTM(isum,iloc,isum_std,iloc_std,gather_rtm,isrc):
    logger.info("preparing RTM: src%d"%isrc) 
    with open("./Input/rec.in","r") as fp:
        nrec0 = int(fp.readline().split()[0])
        component = fp.readline().split()[3]
    nrec = len(iloc)
    nrec_std = len(iloc_std)
    if nrec != nrec_std or nrec0!=nrec:
        logger.warning("nrec not equal! %d(%d)-%d"%(nrec,nrec0,nrec_std)) 
    nt = len(isum[0])
    nt_std = len(isum_std[0])
    if int(nt) != nt_std:
        logger.warning("nt not equal! %d-%d"%(nt,nt_std) )

    if mode == 1:
        with open(os.path.join(rtmdir,'Input','src.in'+'_'+str(isrc).zfill(4)),'w') as fsrc:
            fsrc.write("%d %d\n" % (nrec_std, nt_std))
            for i in range(nrec_std):
                fsrc.write("%d %d %d %s\n"%(iloc[i][0],iloc[i][1],iloc[i][2],component))
            savetxt(fsrc,gather_rtm)

    if mode == 0:
        fn_src = os.path.join('Input','src.in_'+str(isrc).zfill(4))
        with open(fn_src) as fo:
            fo.readline()
            srcinfo = fo.readline()
            [srcx, srcy] = [int(n) for n in srcinfo.split(' ')[:-2]]
        for irec in range(nrec):
            assert (iloc[irec] == iloc_std[irec]).all(), logger.warning("Position NOT correct for 'OUTPUT'(%s) and 'STD/OUTPUT'(%s)"%(iloc(irec),iloc_std(irec)))
            if iloc[irec][0] == srcx and iloc[irec][1] == srcy:
                logger.debug('src%d:(%d,%d); rec%d:(%d,%d)'%(isrc,srcx,srcy,irec,iloc[irec][0],iloc[irec][1]))
                gather0 = gather_rtm[irec,:]
                return nrec,nt,iloc,component,gather0
                break


def pre_RTM(list_src):
    if mode == 0:
        locs = ''
        dats = []

    for i in list_src:#range(nsrc)
        isum,iloc=merge_gather(os.path.join(workdir,'Output'),i)
        isum_std,iloc_std=merge_gather(os.path.join(workdir,'STD','Output'),i)
        gather_rtm=remove_STD(isum,isum_std,i)
        
        if mode == 0:
            nrec,nt,iloc,component,gather0 = prepare_RTM(isum,iloc,isum_std,iloc_std,gather_rtm,i)
            locs += "%d %d %d %s\n"%(iloc[i][0],iloc[i][1],iloc[i][2],component)
            dats.append(gather0)
            if i == list_src[-1]:
                logger.info('Writing zero-offset source data...')
                with open(os.path.join(rtmdir,'Input','src.in_0000'),'w') as fsrc:
                    fsrc.write("%d %d\n" % (nrec, nt))
                    fsrc.write(locs)
                    savetxt(fsrc,dats)
                logger.info('All Done.')
        else:
            prepare_RTM(isum,iloc,isum_std,iloc_std,gather_rtm,i)

        logger.info('Done: src%d\n'%i)


if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "m:o:d:", ["mode=","outdir=","workdir="])
    except getopt.GetoptError as err:
        # print help information and exit:
        logger.error(err)  # will print something like "option -a not recognized"
        # usage()
        sys.exit(2)

    mode = 1
    workdir = ''
    rtmdir_name = 'RTM'
    for o, a in opts:
        if o in ('-d','--workdir'):
            workdir = a
        elif o in ('-o','--outdir'):
            rtmdir_name = a
            logger.info('Output dir = "%s"'%rtmdir_name)
        elif o in ('-m','--mode'):
            if a == '0':
                mode = 0
                nsrc = len([f for f in os.listdir('Input') if os.path.isfile(os.path.join('Input',f)) and 'src.in_' in f])
                logger.info('Zero-offset mode. We have %d sources to process.'%nsrc)
            elif a == '1':
                mode = 1
                logger.info('normal mode. src%d'%isrc)
            else:
                logger.error('mode parameter must be 0,1')
        else:
            assert False, "unhandled option"


    rtmdir = os.path.join(workdir,rtmdir_name)
    read_par
    if mode == 0:
        pre_RTM(range(nsrc))
    elif mode == 1:
        isrc = int(args[0])
        pre_RTM([isrc])

