#!/usr/bin/env python
import sys,os,logging,getopt
import image_condition.get_filename_list as rtmfn
sys.path.append('..')
from common import *

def cleanfiles(filelist):
    for fn in filelist:
        try:
            os.remove(fn)
            logger.debug('Removed: %s'%fn)
        except Exception as e:
            logger.error(e)

def check_file(isrc,workdir,logger,result_dir='Result'):
    logger.info("Checking files of src%d"%isrc)
    dum_s = True
    dum_w = True
    result_dir = os.path.join(workdir,result_dir)
    for xyz in ['x','y','z']:
        if xyz == 'x':
            nslice = slice_nx
        elif xyz == 'y':
            nslice = slice_ny
        elif xyz == 'z':
            nslice = slice_nz
        if dum_s:
            for i in range(nslice):
                dum_s = os.path.isfile(os.path.join(result_dir, rtmfn.result_slice_fn.format(xyz=xyz,isrc=isrc,islice=i,ext='bin'))) and dum_s
                dum_s = os.path.isfile(os.path.join(result_dir, rtmfn.result_slice_nf_fn.format(xyz=xyz,isrc=isrc,islice=i,ext='bin'))) and dum_s
                dum_s = os.path.isfile(os.path.join(result_dir,  rtmfn.result_slice_nb_fn.format(xyz=xyz,isrc=isrc,islice=i,ext='bin'))) and dum_s
                if not dum_s:
                    break
        if dum_w:
            dum_w = os.path.isfile(os.path.join(result_dir, rtmfn.result_wavefield_fn.format(isrc=isrc,ext='bin'))) and dum_w
            dum_w = os.path.isfile(os.path.join(result_dir, rtmfn.result_wavefield_f_fn.format(isrc=isrc,ext='bin'))) and dum_w
            dum_w = os.path.isfile(os.path.join(result_dir, rtmfn.result_wavefield_f_fn.format(isrc=isrc,ext='bin'))) and dum_w
    return dum_s,dum_w

def clean(isrc,workdir,mode):
    # logger
    logger=addlogger('clean',path=os.path.join(workdir,'log'))

    # get files
    logger.info("Getting files of src%d"%isrc)
    path0,listx0,listy0,listz0,listw0 = rtmfn.get_isrc_filenames_data(isrc,workdir)
    list0 = [os.path.join(path0,fn) for fn in listx0+listy0+listz0]
    listw0 = [os.path.join(path0,fn) for fn in listw0]
    path1,listx1,listy1,listz1,listw1 = rtmfn.get_isrc_filenames(isrc,os.path.join(workdir,'STD'))
    list1 = [os.path.join(path1,fn) for fn in listx1+listy1+listz1]
    listw1 = [os.path.join(path1,fn) for fn in listw1]
    slicelist = list0+list1
    wvlist = listw0+listw1
    if 'm' in mode:
        path2,listx2,listy2,listz2,listw2 = rtmfn.get_isrc_filenames(isrc,os.path.join(workdir,'RTM'))
        slicelist += [os.path.join(path2,fn) for fn in listx2+listy2+listz2]
        wvlist += [os.path.join(path2,fn) for fn in listw2]
    # if 'z' in mode:
    #     path3,listx3,listy3,listz3,listw3 = rtmfn.get_isrc_filenames(isrc,os.path.join(workdir,'RTM0'))
    #     slicelist += [os.path.join(path3,fn) for fn in listx3+listy3+listz3]
    #     wvlist += [os.path.join(path3,fn) for fn in listw3]

    logger.info("Begin cleaning files of src%d"%isrc)
    if not mode == 'z':
        check_s,check_w = check_file(isrc,workdir,logger=logger)
        if check_s == True:
            logger.info('cleaning slices: src%d'%isrc)
            cleanfiles(slicelist)
        else:
            logger.warning("src%d: No slices deleted, as results are not complete!"%isrc)
        if check_w == True:
            logger.info('cleaning wavefields: src%d'%isrc)
            cleanfiles(wvlist)
        else:
            logger.warning("src%d: No wavefield deleted, as results are not complete!"%isrc)
    else:
        logger.info('Zero-offset mode, Cleaning all without check: src%d'%isrc)
        cleanfiles(slicelist+wvlist)
    logger.info('clean src%d done.'%isrc)

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "md:", ["mode=","workdir="])
    except getopt.GetoptError as err:
        # print help information and exit:
        raise err  # will print something like "option -a not recognized"
        # usage()
        sys.exit(2)
    # parse arg
    workdir = ''
    mode = 'm'
    for o, a in opts:
        if o in ('-m','--mode'):
            mode = a
        elif o in ('-d','--workdir'):
            workdir = a
        else:
            assert False, "unhandled option"
    isrc = int(args[0])

    # main
    clean(isrc,workdir,mode)


    