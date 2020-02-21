#!/usr/bin/env python
from numpy import *
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import *
import sys,os,re,logging,getopt
from par_RTM import *

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
    if 'win' in sys.platform:
        wavefield1 = os.popen('dir/b/on '+os.path.join(dir1,'Wavefield*dat'+'_'+str(isrc).zfill(4))).readlines()
        wavefield2 = os.popen('dir/b/on '+os.path.join(dir2,'Wavefield*dat'+'_'+str(isrc).zfill(4))).readlines()[::-1]
    elif 'linux' in sys.platform:
        wavefield1 = os.popen('ls ' + os.path.join(dir1,'Wavefield*dat'+'_'+str(isrc).zfill(4))).readlines()
        wavefield2 = os.popen('ls ' + os.path.join(dir2,'Wavefield*dat'+'_'+str(isrc).zfill(4))).readlines()[::-1]
    else:
        logger.error('unknown platform: %s'%sys.platform)
        return 0
    # xlist1 = os.popen('ls '+dir1+'/*xSlice*dat').readlines()
    # ylist1 = os.popen('ls '+dir1+'/*ySlice*dat').readlines()
    # zlist1 = os.popen('ls '+dir1+'/*zSlice*dat').readlines()
    # xlist2 = os.popen('ls '+dir2+'/*xSlice*dat').readlines()[::-1]
    # ylist2 = os.popen('ls '+dir2+'/*ySlice*dat').readlines()[::-1]
    # zlist2 = os.popen('ls '+dir2+'/*zSlice*dat').readlines()[::-1]
    # if len(xlist1) != len(xlist2) or len(ylist1) != len(ylist2) or len(zlist1) != len(zlist2):
    #     print "The number of slices are different!"
    #     return 0
    if len(wavefield1) != len(wavefield2):
        logger.error("src%d: The number of wavefiled are different!(std:%d,rtm:%d)"%(isrc,len(wavefield1),len(wavefield2))) 
        return 0

    for i in range(len(wavefield1)):
        name1 = wavefield1[i].strip('\n')
        name2 = wavefield2[i].strip('\n')
        if 'win' in sys.platform:
            name1 = os.path.join(dir1,name1)
            name2 = os.path.join(dir2,name2)
        logger.debug('corr_wavefield:%s,%s'%(os.path.basename(name1),os.path.basename(name2)))
        data1 = loadtxt(name1)
        data2 = loadtxt(name2)
        if i == 0:
            corr_data = data1 * data2
            data_forward = data1 * data1
            data_backward = data2 * data2
        else:
            corr_data += data1 * data2
            data_forward += data1 * data1
            data_backward += data2 * data2


    savetxt(os.path.join(dir3,'result_wavefield_corr.dat'+'_'+str(isrc).zfill(4)),corr_data)
    savetxt(os.path.join(dir3,'result_wavefield_forward.dat'+'_'+str(isrc).zfill(4)),data_forward)
    savetxt(os.path.join(dir3,'result_wavefield_backward.dat'+'_'+str(isrc).zfill(4)),data_backward)
    return 1


    # print len(corr_data)
    # print nx/step_x_wavefield, ny/step_x_wavefield, nz/step_x_wavefield
    # corr_data = reshape(corr_data,(int(ceil(float(nx)/float(step_x_wavefield))),
    #     int(ceil(float(ny)/float(step_x_wavefield))),
    #     int(ceil(float(nz)/float(step_x_wavefield)))))
    # savetxt('result_wavefield_corr.dat',result)

    # clf()
    # imshow(corr_data[20,:,:].T, interpolation='none')
    # colorbar()
    # savefig('xResult_Wavefield'+'_'+str(isrc).zfill(4)+'.png')
    # clf()
    # imshow(corr_data[:,71,:].T, interpolation='none',vmin = -0.002, vmax = 0.002)
    # colorbar()
    # savefig('yResult_Wavefield'+'_'+str(isrc).zfill(4)+'.png')
    # clf()
    # imshow(corr_data[:,:,29].T, interpolation='none')
    # colorbar()
    # savefig('zResult_Wavefield'+'_'+str(isrc).zfill(4)+'.png')

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
