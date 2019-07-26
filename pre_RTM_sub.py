#!/usr/bin/env python
from numpy import *
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import *
import os,sys,logging,getopt
from par_RTM import *

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
    logger.info("merging gather: %s-src%d"%(idir,isrc)) 
    global isum
    global iloc
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
    savefig(os.path.join(idir, 'merge_gather'+'_'+str(isrc).zfill(4)+'.png'))
    with open(os.path.join(idir, 'merge_gather'+'_'+str(isrc).zfill(4)+'.dat'),'w') as fp:
        savetxt(fp,isum)


def remove_STD(idir,idir_STD,isrc):
    logger.info("removing source: src%d"%isrc) 
    global gather_RTM
    gather = loadtxt(os.path.join(idir,'merge_gather'+'_'+str(isrc).zfill(4)+'.dat'))
    gather_STD = loadtxt(os.path.join(idir_STD,'merge_gather'+'_'+str(isrc).zfill(4)+'.dat'))
    gather_RTM = gather-gather_STD
    gather_RTM = fliplr(gather_RTM)
    # for i in range(len(gather_RTM[:,0])):
    #     gather_RTM[i,:] = gather_RTM[i,:]/max(abs(gather_RTM[i,:]))
    figure()
    imshow(gather_RTM,cmap='gray',origin='lower',extent=(0,nt,0,nt/2))
    savefig(os.path.join(workdir,'gather_without_src'+'_'+str(isrc).zfill(4)+'.png'))


def prepare_RTM(isrc):
    logger.info("preparing RTM: src%d"%isrc) 
    # with file('rec.in','r') as frec:
    with open(os.path.join("Input","rec.in"),"r") as fp:
        nrec = int(fp.readline().split()[0])
        component = fp.readline().split()[3]
    nrec_rtm = len(iloc)
    # nrec = int(frec.readline().split()[0])
    if nrec != nrec_rtm:
        logger.warning("nrec not equal! %d-%d"%(nrec,nrec_rtm)) 
    nt_rtm = len(isum[0])
    if int(nt) != nt_rtm:
        logger.warning("nt not equal! %d-%d"%(nt,nt_rtm) )

    if mode == 1:
        with open(os.path.join(rtmdir,'Input','src.in'+'_'+str(isrc).zfill(4)),'w') as fsrc:
            fsrc.write("%d %d\n" % (nrec_rtm, nt_rtm))
            # fsrc.close()
            # fsrc=open('./RTM/src.in','a')
            for i in range(nrec_rtm):
                fsrc.write("%d %d %d %s\n"%(iloc[i][0],iloc[i][1],iloc[i][2],component))
            savetxt(fsrc,gather_RTM)

    if mode == 0:
        fn_src = os.path.join('Input','src.in_'+str(isrc).zfill(4))
        with open(fn_src) as fo:
            fo.readline()
            srcinfo = fo.readline()
            [srcx,srcy] = [int(n) for n in srcinfo.split(' ')[:-2]]
        for irec in range(nrec_rtm):
            if iloc[irec][0] == srcx and iloc[irec][1] == srcy:
                logger.debug('src%d:(%d,%d); rec%d:(%d,%d)'%(isrc,srcx,srcy,irec,iloc[irec][0],iloc[irec][1]))
                gather0 = gather_RTM[irec,:]
                return nrec_rtm,nt_rtm,iloc,component,gather0
    #             with open(os.path.join('RTM','Input','srcloc.in'),'a+') as floc:
    #                 if isrc == 0:
    #                     floc.write("%d %d\n" % (nrec_rtm, nt_rtm))
    #                 floc.write("%d %d %d %s\n"%(iloc[irec][0],iloc[irec][1],iloc[irec][2],component))
    #             with open(os.path.join('RTM','Input','srcdat.in'),'a+') as fsrc:
    #                 savetxt(fsrc,transpose([gather_RTM[irec,:]]), newline=' ')
    #                 fsrc.write('\n')
    # ## muilti-process writing file may cause unpredicted problem ##
                break


def pre_RTM(list_src):
    if mode == 0:
        locs = ''
        dats = []

    for i in list_src:#range(nsrc)
        merge_gather(os.path.join(workdir,'Output'),i)
        merge_gather(os.path.join(workdir,'STD','Output'),i)
        remove_STD(os.path.join(workdir,'Output'),os.path.join(workdir,'STD','Output'),i)
        
        if mode == 0:
            nrec_rtm,nt_rtm,iloc,component,gather0 = prepare_RTM(i)
            locs += "%d %d %d %s\n"%(iloc[i][0],iloc[i][1],iloc[i][2],component)
            dats.append(gather0)
            if i == list_src[-1]:
                logger.info('Writing zero-offset source data...')
                with open(os.path.join(rtmdir,'Input','src.in_0000'),'w') as fsrc:
                    fsrc.write("%d %d\n" % (nrec_rtm, nt_rtm))
                    fsrc.write(locs)
                    savetxt(fsrc,dats)
                logger.info('All Done.')
        else:
            prepare_RTM(i)

        logger.info('Done prepareing RTM: src%d'%i)



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

