#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys,os,re,logging,getopt
import image_condition.get_filename_list as rtmfn
sys.path.append('..')
from common import *

def corr_wavefield(isrc, workdir, path1,path2,list1,list2, outdir,logger):
    logger.info('corr_wavefield src%d'%isrc)
    
    if len(list1) != len(list2):
        logger.error("src%d: The number of wavefiled are different!(std:%d,rtm:%d)"%(isrc,len(list1),len(list2))) 
        return 0
    
    wvf_nx = nx//step_x_wavefield
    wvf_ny = ny//step_x_wavefield
    wvf_nz = nz//step_x_wavefield
    for i in range(len(list1)):
        name1 = os.path.join(path1,list1[i])
        name2 = os.path.join(path2,list2[i])
        logger.debug('corr_wavefield:%s,%s'%(os.path.basename(name1),os.path.basename(name2)))
        data1 = read_bin_data(name1,[wvf_nx,wvf_ny,wvf_nz])
        data2 = read_bin_data(name2,[wvf_nx,wvf_ny,wvf_nz])
        if i == 0:
            corr_data = data1 * data2
            data_forward = data1 * data1
            data_backward = data2 * data2
        else:
            corr_data += data1 * data2
            data_forward += data1 * data1
            data_backward += data2 * data2

    logger.info('Saving wavefield of src%d'%(isrc))
    # np.savetxt(os.path.join(outdir,rtmfn.result_wavefield_fn.format(isrc=isrc)),np.reshape(corr_data,[wvf_nx*wvf_ny,wvf_nz]))
    # np.savetxt(os.path.join(outdir,rtmfn.result_wavefield_f_fn.format(isrc=isrc)),np.reshape(data_forward,[wvf_nx*wvf_ny,wvf_nz]))
    # np.savetxt(os.path.join(outdir,rtmfn.result_wavefield_b_fn.format(isrc=isrc)),np.reshape(data_backward,[wvf_nx*wvf_ny,wvf_nz]))
    corr_data.astype('float32').tofile(os.path.join(outdir,rtmfn.result_wavefield_fn.format(isrc=isrc,ext='bin')))
    data_forward.astype('float32').tofile(os.path.join(outdir,rtmfn.result_wavefield_f_fn.format(isrc=isrc,ext='bin')))
    data_backward.astype('float32').tofile(os.path.join(outdir,rtmfn.result_wavefield_b_fn.format(isrc=isrc,ext='bin')))
    logger.info('Figureing wavefield of src%d'%(isrc))
    for i in range(slice_nx):
        plt.clf()
        plt.imshow(corr_data[int(slice_x[i]/step_x_wavefield),:,:].T, interpolation='none')
        plt.colorbar()
        plt.savefig(os.path.join(outdir,rtmfn.result_wavefield_fn.format(isrc=isrc,ext='slx%d.png'%i)))
    for i in range(slice_ny):
        plt.clf()
        plt.imshow(corr_data[:,int(slice_y[i]/step_x_wavefield),:].T, interpolation='none')#,vmin = -0.002, vmax = 0.002)
        plt.colorbar()
        plt.savefig(os.path.join(outdir,rtmfn.result_wavefield_fn.format(isrc=isrc,ext='sly%d.png'%i)))
    for i in range(slice_nz):
        plt.clf()
        plt.imshow(corr_data[:,:,int(slice_z[i]/step_x_wavefield)].T, interpolation='none')
        plt.colorbar()
        plt.savefig(os.path.join(outdir,rtmfn.result_wavefield_fn.format(isrc=isrc,ext='slz%d.png'%i)))

    logger.info('corr_wavefield src%d done.'%isrc)
    return 0

def corr_RTM_wavefield_sub(isrc, workdir, path1,path2,list1,list2,outdir='Result'):
    # logger
    logger=addlogger('corr_RTM_wavefield_sub',path=os.path.join(workdir,'log'))
    # mkdir
    outdir = os.path.join(workdir,outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    # main
    corr_wavefield(isrc, workdir, path1,path2,list1,list2,outdir,logger)

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
    corr_RTM_wavefield_sub(isrc,workdir,path1,path2,listsw1,listsw2)