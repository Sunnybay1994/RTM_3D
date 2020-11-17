#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys,os,re,logging,getopt,struct
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

def get_fn_from_dir(fn, order=1):
    logger.info('loading files with name "' + os.path.basename(fn) + '" from "' + os.path.dirname(fn) + '"')
    if 'win' in sys.platform:
        list1 = os.popen('dir/b/on ' + fn).readlines()
        dirname = os.path.dirname(fn)
        list1 = list(map(lambda x: os.path.join(dirname, x), list1))
    elif 'linux' in sys.platform:
        list1 = os.popen('ls ' + fn).readlines()
    else:
        logger.error('unknown platform: %s'%sys.platform)
        return 0
    if order == -1:
        list1.reverse()

    return list(map(lambda x: x.strip('\n'), list1))

def read_bin_data(fn, dims):
    with open(fn,'rb') as fo:
        data_raw = struct.unpack('f'*dims[0]*dims[1]*dims[2],fo.read())
    dims.reverse()
    dims = np.array(dims)
    data = np.reshape(data_raw,dims[dims!=1])
    return data


def corr_slice(isrc, workdir, dir1 = 'Output', dir2 = os.path.join('RTM','Output'), dir3 = 'Result'):
    logger.info('corr_slice src%d'%isrc)

    dir1 = os.path.join(workdir,dir1)
    dir2 = os.path.join(workdir,dir2)
    dir3 = os.path.join(workdir,dir3)
    
    xlist1 = get_fn_from_dir(os.path.join(dir1, 'slx_Ey_'+str(0).zfill(4)+'*.bin'))
    ylist1 = get_fn_from_dir(os.path.join(dir1, 'sly_Ey_'+str(0).zfill(4)+'*.bin'))
    zlist1 = get_fn_from_dir(os.path.join(dir1, 'slz_Ey_'+str(0).zfill(4)+'*.bin'))

    xlist2 = get_fn_from_dir(os.path.join(dir2, 'slx_Ey_'+str(0).zfill(4)+'*.bin'), -1)
    ylist2 = get_fn_from_dir(os.path.join(dir2, 'sly_Ey_'+str(0).zfill(4)+'*.bin'), -1)
    zlist2 = get_fn_from_dir(os.path.join(dir2, 'slz_Ey_'+str(0).zfill(4)+'*.bin'), -1)

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
        xname1 = xlist1[i]
        yname1 = ylist1[i]
        zname1 = zlist1[i]
        xname2 = xlist2[i]
        yname2 = ylist2[i]
        zname2 = zlist2[i]
        logger.debug('corr_slice:%s,%s'%(os.path.basename(xname1),os.path.basename(xname2)))
        xdata1 = read_bin_data(xname1,[slice_nx,ny,nz])
        ydata1 = read_bin_data(yname1,[nx,slice_ny,nz])
        zdata1 = read_bin_data(zname1,[nx,ny,slice_nz])
        xdata2 = read_bin_data(xname2,[slice_nx,ny,nz])
        ydata2 = read_bin_data(yname2,[nx,slice_ny,nz])
        zdata2 = read_bin_data(zname2,[nx,ny,slice_nz])
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


    np.savetxt(os.path.join(dir3, 'result_xcorr.dat'+'_'+str(isrc).zfill(4)),np.reshape(xcorr_data,[slice_nx*ny,nz]))
    np.savetxt(os.path.join(dir3, 'result_ycorr.dat'+'_'+str(isrc).zfill(4)),np.reshape(ycorr_data,[nx*slice_ny,nz]))
    np.savetxt(os.path.join(dir3, 'result_zcorr.dat'+'_'+str(isrc).zfill(4)),np.reshape(zcorr_data,[nx*ny,slice_nz]))
    # xslice = np.reshape(xcorr_data,(slice_nx,ny,nz))
    # yslice = np.reshape(ycorr_data,(slice_ny,nx,nz))
    # zslice = np.reshape(zcorr_data,(slice_nz,nx,ny))
    xslice = xcorr_data
    yslice = ycorr_data
    zslice = zcorr_data

    np.savetxt(os.path.join(dir3, 'result_xcorr_normal_forward.dat'+'_'+str(isrc).zfill(4)),np.reshape(xcorr_normal_forward,[slice_nx*ny,nz]))
    np.savetxt(os.path.join(dir3, 'result_ycorr_normal_forward.dat'+'_'+str(isrc).zfill(4)),np.reshape(ycorr_normal_forward,[nx*slice_ny,nz]))
    np.savetxt(os.path.join(dir3, 'result_zcorr_normal_forward.dat'+'_'+str(isrc).zfill(4)),np.reshape(zcorr_normal_forward,[nx*ny,slice_nz]))
    # xslice_normal_forward = np.reshape(xcorr_normal_forward,(slice_nx,ny,nz))
    # yslice_normal_forward = np.reshape(ycorr_normal_forward,(slice_ny,nx,nz))
    # zslice_normal_forward = np.reshape(zcorr_normal_forward,(slice_nz,nx,ny))
    xslice_normal_forward = xcorr_normal_forward
    yslice_normal_forward = ycorr_normal_forward
    zslice_normal_forward = zcorr_normal_forward

    np.savetxt(os.path.join(dir3, 'result_xcorr_normal_backward.dat'+'_'+str(isrc).zfill(4)),np.reshape(xcorr_normal_backward,[slice_nx*ny,nz]))
    np.savetxt(os.path.join(dir3, 'result_ycorr_normal_backward.dat'+'_'+str(isrc).zfill(4)),np.reshape(ycorr_normal_backward,[nx*slice_ny,nz]))
    np.savetxt(os.path.join(dir3, 'result_zcorr_normal_backward.dat'+'_'+str(isrc).zfill(4)),np.reshape(zcorr_normal_backward,[nx*ny,slice_nz]))
    # xslice_normal_backward = np.reshape(xcorr_normal_backward,(slice_nx,ny,nz))
    # yslice_normal_backward = np.reshape(ycorr_normal_backward,(slice_ny,nx,nz))
    # zslice_normal_backward = np.reshape(zcorr_normal_backward,(slice_nz,nx,ny))
    xslice_normal_backward = xcorr_normal_backward
    yslice_normal_backward = ycorr_normal_backward
    zslice_normal_backward = zcorr_normal_backward

    plt.clf()
    plt.imshow(xslice[0].T, interpolation='none')
    plt.colorbar()
    plt.savefig(os.path.join(dir3, 'xResult_RTM'+'_'+str(isrc).zfill(4)+'.png'))
    plt.clf()
    plt.imshow(yslice[0].T, interpolation='none')
    plt.colorbar()
    plt.savefig(os.path.join(dir3, 'yResult_RTM'+'_'+str(isrc).zfill(4)+'.png'))
    plt.clf()
    plt.imshow(zslice[0].T, interpolation='none')
    plt.colorbar()
    plt.savefig(os.path.join(dir3, 'zResult_RTM'+'_'+str(isrc).zfill(4)+'.png'))

    plt.clf()
    plt.imshow(xslice_normal_forward[0].T, interpolation='none')
    plt.colorbar()
    plt.savefig(os.path.join(dir3, 'xResult_RTM_normal_forward'+'_'+str(isrc).zfill(4)+'.png'))
    plt.clf()
    plt.imshow(yslice_normal_forward[0].T, interpolation='none')
    plt.colorbar()
    plt.savefig(os.path.join(dir3, 'yResult_RTM_normal_forward'+'_'+str(isrc).zfill(4)+'.png'))
    plt.clf()
    plt.imshow(zslice_normal_forward[0].T, interpolation='none')
    plt.colorbar()
    plt.savefig(os.path.join(dir3, 'zResult_RTM_normal_forward'+'_'+str(isrc).zfill(4)+'.png'))

    plt.clf()
    plt.imshow(xslice_normal_backward[0].T, interpolation='none')
    plt.colorbar()
    plt.savefig(os.path.join(dir3, 'xResult_RTM_normal_backward'+'_'+str(isrc).zfill(4)+'.png'))
    plt.clf()
    plt.imshow(yslice_normal_backward[0].T, interpolation='none')
    plt.colorbar()
    plt.savefig(os.path.join(dir3, 'yResult_RTM_normal_backward'+'_'+str(isrc).zfill(4)+'.png'))
    plt.clf()
    plt.imshow(zslice_normal_backward[0].T, interpolation='none')
    plt.colorbar()
    plt.savefig(os.path.join(dir3, 'zResult_RTM_normal_backward'+'_'+str(isrc).zfill(4)+'.png'))
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

    idir = os.path.join(workdir,'Result')
    if not os.path.exists(idir):
        os.mkdir(idir)

    corr_slice(isrc,workdir)