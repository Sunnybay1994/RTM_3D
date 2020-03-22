import os,sys,logging
from numpy import *

#logger
logger = logging.getLogger('extract_zero_offset')
logger.setLevel(logging.DEBUG) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
fh = logging.FileHandler('log/extract_zero_offset.log')
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s(%(process)d-%(processName)s): (%(levelname)s) %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)

def merge_gather(idir, isrc):
    logger.info('Mering gather: %s - src%d'%(idir,isrc))
    if 'win' in sys.platform:
        ilist = os.popen('dir/b/on '+ os.path.join(idir,'gather*dat'+'_'+str(isrc).zfill(4))).readlines()
        ilist = [os.path.join(idir,f) for f in ilist]
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
    # figure()
    # for i in range(len(isum)):
    #     plot(isum[i]/max(abs(isum[i]))+i)
    # savefig(idir+'/merge_gather'+'_'+str(isrc).zfill(4)+'.png')
    with open(idir + '/merge_gather'+'_'+str(isrc).zfill(4)+'.dat','w+') as fp:
        savetxt(fp,isum)
    logger.debug('done.')

def extract_zero_offset(Outdir='Zero_offset',datadir='Output'):
    logger.info('extracting zero offset data.')
    indir = 'Input'
    fn_out = os.path.join(Outdir,'Input','src.in_0000')

    # load rec pos
    logger.info('Loading rec loc info...')
    with open(os.path.join(indir,'rec.in')) as fo:
        nrec = int(fo.readline())
        # src and rec pos should both be y increaing first.
        recpos = [[int(x) for x in line.split(' ')[:-2]] for line in fo.readlines()]
    
    nsrc = len([f for f in os.listdir(indir) if os.path.isfile(os.path.join(indir,f)) and 'src.in_' in f])
    irec = -1
    srcinfos = ''
    datas = ''
    nt = 0
    for isrc in range(nsrc):
        logger.info('Extracting src%d of %d...'%(isrc,nsrc))
        fsrc = os.path.join(indir,'src.in_'+str(isrc).zfill(4))
        # logger.debug(fsrc)
        # match src pos
        irec = 0
        with open(fsrc) as fo:
            fo.readline()
            srcinfo = fo.readline()
            srcinfos += srcinfo
            [x,y] = [int(n) for n in srcinfo.split(' ')[:-2]]
            [xr,yr] = [int(n) for n in recpos[irec]]
            while xr<x or yr<y: # rec loc shouble be exactly the same as src loc
                irec += 1
                [xr,yr] = [int(n) for n in recpos[irec]]
            logger.debug('src%d:(%d,%d); rec%d:(%d,%d)'%(isrc,x,y,irec,xr,yr))

        # read zero offset data and write file
        fgather = os.path.join(datadir,'merge_gather_'+str(isrc).zfill(4)+'.dat')
        # logger.debug(fgather)
        if not os.path.isfile(fgather):
            merge_gather(datadir,isrc)

        with open(fgather) as fgo:
            lines = fgo.readlines()
            data = lines[irec].strip().split()
            data.reverse() # reverse data for rtm input
            if isrc == 0:
                nt = len(data)
            datas += ' '.join(data)+'\n'

    # ourput file objct for zero_offset_rtm input
    logger.info('Output file:%s'%fn_out)
    with open(fn_out,'w') as fo_out:
        fo_out.write(str(nsrc)+' '+str(nt)+'\n')
        fo_out.write(srcinfos)
        fo_out.write(datas)   

if __name__ == '__main__':
    extract_zero_offset('RTM')