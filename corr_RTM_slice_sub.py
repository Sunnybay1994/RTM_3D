#!/usr/bin/env python
from numpy import *
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import *
import sys,os,re,logging,getopt
from par_RTM import *

#logger
logger = logging.getLogger('corr_slice')
logger.setLevel(logging.DEBUG) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
fh = logging.FileHandler('log/corr_RTM_slice_sub.log')
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s(%(process)d-%(processName)s): (%(levelname)s) %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)

def corr_slice(isrc, workdir, dir1 = os.path.join('STD','Output'), dir2 = os.path.join('RTM','Output'), dir3 = 'Result'):
    dir1 = os.path.join(workdir,dir1)
    dir2 = os.path.join(workdir,dir2)
    dir3 = os.path.join(workdir,dir3)
    logger.info('corr_slice src%d'%isrc)
    if 'win' in sys.platform:
        xlist1 = os.popen('dir/b/on '+os.path.join(dir1, '*xSlice*dat'+'_'+str(isrc).zfill(4))).readlines()
        ylist1 = os.popen('dir/b/on '+os.path.join(dir1, '*ySlice*dat'+'_'+str(isrc).zfill(4))).readlines()
        zlist1 = os.popen('dir/b/on '+os.path.join(dir1, '*zSlice*dat'+'_'+str(isrc).zfill(4))).readlines()
        xlist2 = os.popen('dir/b/on '+os.path.join(dir2, '*xSlice*dat'+'_'+str(isrc).zfill(4))).readlines()[::-1]
        ylist2 = os.popen('dir/b/on '+os.path.join(dir2, '*ySlice*dat'+'_'+str(isrc).zfill(4))).readlines()[::-1]
        zlist2 = os.popen('dir/b/on '+os.path.join(dir2, '*zSlice*dat'+'_'+str(isrc).zfill(4))).readlines()[::-1]
    elif 'linux' in sys.platform:
        xlist1 = os.popen('ls '+os.path.join(dir1,'*xSlice*dat'+'_'+str(isrc).zfill(4))).readlines()
        ylist1 = os.popen('ls '+os.path.join(dir1,'*ySlice*dat'+'_'+str(isrc).zfill(4))).readlines()
        zlist1 = os.popen('ls '+os.path.join(dir1,'*zSlice*dat'+'_'+str(isrc).zfill(4))).readlines()
        xlist2 = os.popen('ls '+os.path.join(dir2,'*xSlice*dat'+'_'+str(isrc).zfill(4))).readlines()[::-1]
        ylist2 = os.popen('ls '+os.path.join(dir2,'*ySlice*dat'+'_'+str(isrc).zfill(4))).readlines()[::-1]
        zlist2 = os.popen('ls '+os.path.join(dir2,'*zSlice*dat'+'_'+str(isrc).zfill(4))).readlines()[::-1]
    else:
        logger.error('unknown platform: %s'%sys.platform)
        return 0
    if len(xlist1) != len(xlist2) or len(ylist1) != len(ylist2) or len(zlist1) != len(zlist2):
        logger.error("src%d: The number of slice are different!(std:%d %d %d, rtm:%d %d %d)"%(isrc,len(xlist1),len(ylist1),len(zlist1),len(xlist2),len(ylist2),len(zlist2))) 
        return 0

    xdata1 = []
    ydata1 = []
    zdata1 = []
    xdata2 = []
    ydata2 = []
    zdata2 = []
    dum = 0
    for i in range(len(xlist1)):
        xname1 = xlist1[i].strip('\n')
        yname1 = ylist1[i].strip('\n')
        zname1 = zlist1[i].strip('\n')
        xname2 = xlist2[i].strip('\n')
        yname2 = ylist2[i].strip('\n')
        zname2 = zlist2[i].strip('\n')
        if 'win' in sys.platform:
            xname1 = os.path.join(dir1,xname1)
            yname1 = os.path.join(dir1,yname1)
            zname1 = os.path.join(dir1,zname1)
            xname2 = os.path.join(dir2,xname2)
            yname2 = os.path.join(dir2,yname2)
            zname2 = os.path.join(dir2,zname2)
        logger.debug('corr_slice:%s,%s'%(os.path.basename(xname1),os.path.basename(xname2)))
        xdata1 = loadtxt(xname1)
        ydata1 = loadtxt(yname1)
        zdata1 = loadtxt(zname1)
        xdata2 = loadtxt(xname2)
        ydata2 = loadtxt(yname2)
        zdata2 = loadtxt(zname2)
        if i == 0:
            xcorr_data = xdata1 * xdata2
            ycorr_data = ydata1 * ydata2
            zcorr_data = zdata1 * zdata2

            xnormal_forward = xdata1 * xdata1
            ynormal_forward = ydata1 * ydata1
            znormal_forward = zdata1 * zdata1

            xnormal_backward = xdata2 * xdata2
            ynormal_backward = ydata2 * ydata2
            znormal_backward = zdata2 * zdata2    
        else:
            xcorr_data += xdata1 * xdata2
            ycorr_data += ydata1 * ydata2
            zcorr_data += zdata1 * zdata2

            xnormal_forward += xdata1 * xdata1
            ynormal_forward += ydata1 * ydata1
            znormal_forward += zdata1 * zdata1

            xnormal_backward += xdata2 * xdata2
            ynormal_backward += ydata2 * ydata2
            znormal_backward += zdata2 * zdata2
    try:
        xcorr_normal_forward = xcorr_data / xnormal_forward
    except Warning as wa:
        logger.warning(wa+'(xf:%d/%d)'%(xcorr_data,xnormal_forward))
    try:
        ycorr_normal_forward = ycorr_data / ynormal_forward
    except Warning as wa:
        logger.warning(wa+'(yf:%d/%d)'%(ycorr_data,ynormal_forward))
    try:
        zcorr_normal_forward = zcorr_data / znormal_forward
    except Warning as wa:
        logger.warning(wa+'(zf:%d/%d)'%(zcorr_data,znormal_forward))
    try:
        xcorr_normal_backward = xcorr_data / xnormal_backward
    except Warning as wa:
        logger.warning(wa+'(xb:%d/%d)'%(xcorr_data,xnormal_backward))
    try:
        ycorr_normal_backward = ycorr_data / ynormal_backward
    except Warning as wa:
        logger.warning(wa+'(yb:%d/%d)'%(ycorr_data,ynormal_backward))
    try:
        zcorr_normal_backward = zcorr_data / znormal_backward
    except Warning as wa:
        logger.warning(wa+'(zb:%d/%d)'%(zcorr_data,znormal_backward))


    savetxt(os.path.join(dir3, 'result_xcorr.dat'+'_'+str(isrc).zfill(4)),xcorr_data)
    savetxt(os.path.join(dir3, 'result_ycorr.dat'+'_'+str(isrc).zfill(4)),ycorr_data)
    savetxt(os.path.join(dir3, 'result_zcorr.dat'+'_'+str(isrc).zfill(4)),zcorr_data)
    xslice = reshape(xcorr_data,(slice_nx,ny,nz))
    yslice = reshape(ycorr_data,(slice_ny,nx,nz))
    zslice = reshape(zcorr_data,(slice_nz,nx,ny))

    savetxt(os.path.join(dir3, 'result_xcorr_normal_forward.dat'+'_'+str(isrc).zfill(4)),xcorr_normal_forward)
    savetxt(os.path.join(dir3, 'result_ycorr_normal_forward.dat'+'_'+str(isrc).zfill(4)),ycorr_normal_forward)
    savetxt(os.path.join(dir3, 'result_zcorr_normal_forward.dat'+'_'+str(isrc).zfill(4)),zcorr_normal_forward)
    xslice_normal_forward = reshape(xcorr_normal_forward,(slice_nx,ny,nz))
    yslice_normal_forward = reshape(ycorr_normal_forward,(slice_ny,nx,nz))
    zslice_normal_forward = reshape(zcorr_normal_forward,(slice_nz,nx,ny))

    savetxt(os.path.join(dir3, 'result_xcorr_normal_backward.dat'+'_'+str(isrc).zfill(4)),xcorr_normal_backward)
    savetxt(os.path.join(dir3, 'result_ycorr_normal_backward.dat'+'_'+str(isrc).zfill(4)),ycorr_normal_backward)
    savetxt(os.path.join(dir3, 'result_zcorr_normal_backward.dat'+'_'+str(isrc).zfill(4)),zcorr_normal_backward)
    xslice_normal_backward = reshape(xcorr_normal_backward,(slice_nx,ny,nz))
    yslice_normal_backward = reshape(ycorr_normal_backward,(slice_ny,nx,nz))
    zslice_normal_backward = reshape(zcorr_normal_backward,(slice_nz,nx,ny))

    clf()
    imshow(xslice[0].T, interpolation='none')
    colorbar()
    savefig(os.path.join(dir3, 'xResult_RTM'+'_'+str(isrc).zfill(4)+'.png'))
    clf()
    imshow(yslice[0].T, interpolation='none')
    colorbar()
    savefig(os.path.join(dir3, 'yResult_RTM'+'_'+str(isrc).zfill(4)+'.png'))
    clf()
    imshow(zslice[0].T, interpolation='none')
    colorbar()
    savefig(os.path.join(dir3, 'zResult_RTM'+'_'+str(isrc).zfill(4)+'.png'))

    clf()
    imshow(xslice_normal_forward[0].T, interpolation='none')
    colorbar()
    savefig(os.path.join(dir3, 'xResult_RTM_normal_forward'+'_'+str(isrc).zfill(4)+'.png'))
    clf()
    imshow(yslice_normal_forward[0].T, interpolation='none')
    colorbar()
    savefig(os.path.join(dir3, 'yResult_RTM_normal_forward'+'_'+str(isrc).zfill(4)+'.png'))
    clf()
    imshow(zslice_normal_forward[0].T, interpolation='none')
    colorbar()
    savefig(os.path.join(dir3, 'zResult_RTM_normal_forward'+'_'+str(isrc).zfill(4)+'.png'))

    clf()
    imshow(xslice_normal_backward[0].T, interpolation='none')
    colorbar()
    savefig(os.path.join(dir3, 'xResult_RTM_normal_backward'+'_'+str(isrc).zfill(4)+'.png'))
    clf()
    imshow(yslice_normal_backward[0].T, interpolation='none')
    colorbar()
    savefig(os.path.join(dir3, 'yResult_RTM_normal_backward'+'_'+str(isrc).zfill(4)+'.png'))
    clf()
    imshow(zslice_normal_backward[0].T, interpolation='none')
    colorbar()
    savefig(os.path.join(dir3, 'zResult_RTM_normal_backward'+'_'+str(isrc).zfill(4)+'.png'))
    return 1


    # if(raw_input('Clear Slices?(Y)') == 'Y'):
    #     os.system('rm -f '+dir1+'/*Slice*dat')
    #     os.system('rm -f '+dir2+'/*Slice*dat')

if __name__ == "__main__":
    # nsrc = int(raw_input("How many sources? "))
    # list_src = arange(200,210)
    # #list_src = 56
    # read_par()
    # ipool = Pool(12)
    # ipool.map(corr_slice, list_src)
    # #corr_slice(list_src)
    # ipool.close()
    # ipool.join()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:", ["workdir="])
    except getopt.GetoptError as err:
        # print help information and exit:
        logger.error(err)  # will print something like "option -a not recognized"
        # usage()
        sys.exit(2)

    workdir = os.getcwd()
    for o, a in opts:
        if o in ('-d','--workdir'):
            workdir = a
        else:
            assert False, "unhandled option"

    isrc = int(args[0])

    idir = 'Result'
    if not os.path.exists(idir):
        os.mkdir(idir)

    corr_slice(isrc,workdir)