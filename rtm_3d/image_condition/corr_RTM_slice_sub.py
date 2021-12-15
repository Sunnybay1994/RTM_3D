#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys,os,re,logging,getopt
import image_condition.get_filename_list as rtmfn
sys.path.append('..')
from common import *

def corr_slice(isrc, workdir, path1,path2,xlist1,xlist2,ylist1,ylist2,zlist1,zlist2, outdir,logger):
    logger.info('corr_slice src%d'%isrc)

    if len(xlist1) != len(xlist2) or len(ylist1) != len(ylist2) or len(zlist1) != len(zlist2):
        logger.error("src%d: The number of slice are different!(std:%d %d %d, rtm:%d %d %d)"%(isrc,len(xlist1),len(ylist1),len(zlist1),len(xlist2),len(ylist2),len(zlist2))) 
        return 0

    xdata1 = []
    ydata1 = []
    zdata1 = []
    xdata2 = []
    ydata2 = []
    zdata2 = []
    for i in range(len(xlist1)):
        xname1 = os.path.join(path1,xlist1[i])
        yname1 = os.path.join(path1,ylist1[i])
        zname1 = os.path.join(path1,zlist1[i])
        xname2 = os.path.join(path2,xlist2[i])
        yname2 = os.path.join(path2,ylist2[i])
        zname2 = os.path.join(path2,zlist2[i])
        logger.debug('corr_slice:%s,%s'%(os.path.basename(xname1),os.path.basename(xname2)))
        logger.debug('corr_slice:%s,%s'%(os.path.basename(yname1),os.path.basename(yname2)))
        logger.debug('corr_slice:%s,%s'%(os.path.basename(zname1),os.path.basename(zname2)))
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

    def savefiles(xyz):
        if xyz == 'x':
            nslice = slice_nx
        elif xyz == 'y':
            nslice = slice_ny
        elif xyz == 'z':
            nslice = slice_nz
        for i in range(nslice):
            logger.info('Saving slice%s%d of src%d'%(xyz,i,isrc))
            if xyz == 'x':
                slice_data = xcorr_data[i,:,:]
                slice_normal_forward = xcorr_normal_forward[i,:,:]
                slice_normal_backward = xcorr_normal_backward[i,:,:]
            elif xyz == 'y':
                slice_data = ycorr_data[:,i,:]
                slice_normal_forward = ycorr_normal_forward[:,i,:]
                slice_normal_backward = ycorr_normal_backward[:,i,:]
            elif xyz == 'z':
                slice_data = zcorr_data[:,:,i]
                slice_normal_forward = zcorr_normal_forward[:,:,i]
                slice_normal_backward = zcorr_normal_backward[:,:,i]
            # np.savetxt(os.path.join(outdir, rtmfn.result_slice_fn.format(xyz=xyz,isrc=isrc,islice=i)), slice_data)
            # np.savetxt(os.path.join(outdir, rtmfn.result_slice_nf_fn.format(xyz=xyz,isrc=isrc,islice=i)), slice_normal_forward)
            # np.savetxt(os.path.join(outdir, rtmfn.result_slice_nb_fn.format(xyz=xyz,isrc=isrc,islice=i)), slice_normal_backward)
            slice_data.astype('float32').tofile(os.path.join(outdir, rtmfn.result_slice_fn.format(xyz=xyz,isrc=isrc,islice=i,ext='bin')))
            slice_normal_forward.astype('float32').tofile(os.path.join(outdir, rtmfn.result_slice_nf_fn.format(xyz=xyz,isrc=isrc,islice=i,ext='bin')))
            slice_normal_backward.astype('float32').tofile(os.path.join(outdir, rtmfn.result_slice_nb_fn.format(xyz=xyz,isrc=isrc,islice=i,ext='bin')))
            logger.info('Figureing slice%s%d of src%d'%(xyz,i,isrc))
            plt.clf()
            plt.imshow(slice_data.T, interpolation='none')
            plt.colorbar()
            plt.savefig(os.path.join(outdir, rtmfn.result_slice_fn.format(xyz=xyz,isrc=isrc,islice=i,ext='png')))
            plt.clf()
            plt.imshow(slice_normal_forward.T, interpolation='none')
            plt.colorbar()
            plt.savefig(os.path.join(outdir, rtmfn.result_slice_nf_fn.format(xyz=xyz,isrc=isrc,islice=i,ext='png')))
            plt.clf()
            plt.imshow(slice_normal_backward.T, interpolation='none')
            plt.colorbar()
            plt.savefig(os.path.join(outdir,  rtmfn.result_slice_nb_fn.format(xyz=xyz,isrc=isrc,islice=i,ext='png')))
    for xyz in ['x','y','z']:
        savefiles(xyz)
    logger.info('corr_slice src%d done.'%isrc)
    return 1

def corr_RTM_slice_sub(isrc, workdir, path1,path2,xlist1,xlist2,ylist1,ylist2,zlist1,zlist2, outdir='Result'):
    # logger
    logger=addlogger('corr_RTM_slice_sub',path=os.path.join(workdir,'log'))
    # mkdir
    outdir = os.path.join(workdir,outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    # main
    corr_slice(isrc, workdir, path1,path2,xlist1,xlist2,ylist1,ylist2,zlist1,zlist2, outdir,logger)

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:", ["workdir="])
    except getopt.GetoptError as err:
        # print help information and exit:
        raise err  # will print something like "option -a not recognized"
        # usage()
        sys.exit(2)
    # parse arg
    workdir = os.getcwd()
    for o, a in opts:
        if o in ('-d','--workdir'):
            workdir = a
        else:
            assert False, "unhandled option"
    isrc = int(args[0])
    # get files
    path1,path2,listsx1,listsx2,listy1,listy2,listz1,listz2,listw1,listw2 = rtmfn.get_isrc_filenames_rtm(isrc,workdir)
    # main
    corr_RTM_slice_sub(isrc,workdir,path1,path2,listsx1,listsx2,listy1,listy2,listz1,listz2)