#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys,os,re,logging,getopt
from par_RTM import *
from corr_RTM_slice_sub import get_fn_from_dir,read_bin_data

#logger
logger = logging.getLogger('corr_wavefield')
logger.setLevel(logging.DEBUG) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
fh = logging.FileHandler('log/corr_RTM_wavefield_sub.log')
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s(%(process)d-%(processName)s): (%(levelname)s) %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)

def corr_wavefield(isrc, workdir, dir1 = os.path.join('STD','Output'), dir2 = os.path.join('RTM','Output'), dir3 = 'Result'):
    dir1 = os.path.join(workdir,dir1)
    dir2 = os.path.join(workdir,dir2)
    dir3 = os.path.join(workdir,dir3)
    logger.info('corr_wavefield src%d'%isrc)
    
    list1 = get_fn_from_dir(os.path.join(dir1, 'wvf_Ey_'+str(isrc).zfill(4)+'*.bin'))
    list2 = get_fn_from_dir(os.path.join(dir2, 'wvf_Ey_'+str(isrc).zfill(4)+'*.bin'),-1)
    
    if len(list1) != len(list2):
        logger.error("src%d: The number of wavefiled are different!(std:%d,rtm:%d)"%(isrc,len(list1),len(list2))) 
        return 0
    
    wvf_nx = nx//step_x_wavefield
    wvf_ny = ny//step_x_wavefield
    wvf_nz = nz//step_x_wavefield
    for i in range(len(list1)):
        name1 = list1[i]
        name2 = list2[i]
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


    np.savetxt(os.path.join(dir3,'result_wavefield_corr_' + str(isrc).zfill(4) + '.dat'),np.reshape(corr_data,[wvf_nx*wvf_ny,wvf_nz]))
    np.savetxt(os.path.join(dir3,'result_wavefield_forward_' + str(isrc).zfill(4) + '.dat'),np.reshape(data_forward,[wvf_nx*wvf_ny,wvf_nz]))
    np.savetxt(os.path.join(dir3,'result_wavefield_backward_' + str(isrc).zfill(4) + '.dat'),np.reshape(data_backward,[wvf_nx*wvf_ny,wvf_nz]))

    # print len(corr_data)
    # print nx/step_x_wavefield, ny/step_x_wavefield, nz/step_x_wavefield
    # corr_data = reshape(corr_data,(int(ceil(float(nx)/float(step_x_wavefield))),
    #     int(ceil(float(ny)/float(step_x_wavefield))),
    #     int(ceil(float(nz)/float(step_x_wavefield)))))
    # savetxt('result_wavefield_corr.dat',result)

    # print('Figuring wavefield.')
    plt.clf()
    plt.imshow(corr_data[int(slice_x[0]/step_x_wavefield),:,:].T, interpolation='none')
    plt.colorbar()
    plt.savefig(os.path.join(dir3,'xResult_Wavefield'+'_'+str(isrc).zfill(4)+'.png'))
    plt.clf()
    plt.imshow(corr_data[:,int(slice_y[0]/step_x_wavefield),:].T, interpolation='none')#,vmin = -0.002, vmax = 0.002)
    plt.colorbar()
    plt.savefig(os.path.join(dir3,'yResult_Wavefield'+'_'+str(isrc).zfill(4)+'.png'))
    plt.clf()
    plt.imshow(corr_data[:,:,int(slice_z[0]/step_x_wavefield)].T, interpolation='none')
    plt.colorbar()
    plt.savefig(os.path.join(dir3,'zResult_Wavefield'+'_'+str(isrc).zfill(4)+'.png'))

    return 0

if __name__ == "__main__":
    # nsrc = int(raw_input("How many sources? "))
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
    
    corr_wavefield(isrc,workdir)

    # list_src = range(200,210)

    # #list_src = 172
    # print list_src

    # read_par()
    # ipool = Pool(12)
    # ipool.map(corr_wavefield, list_src)
    # #corr_wavefield(list_src)
    # ipool.close()
    # ipool.join()
