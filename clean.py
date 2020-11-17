#!/usr/bin/env python
import sys,os,logging,getopt
from corr_RTM_slice_sub import get_fn_from_dir

#logger
logger = logging.getLogger('clean')
logger.setLevel(logging.DEBUG) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
fh = logging.FileHandler('log/clean.log')
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s(%(process)d-%(processName)s): (%(levelname)s) %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)

def clean_slice(isrc, workdir, dir1 = os.path.join('STD','Output'), dir2 = os.path.join('RTM','Output'), dir3 = 'Output'):
    dir1 = os.path.join(workdir,dir1)
    dir2 = os.path.join(workdir,dir2)
    dir3 = os.path.join(workdir,dir3)
    logger.info("Begin cleaning slice of src%d"%isrc) 
    
    list1 = get_fn_from_dir(os.path.join(dir1, 'sl*_Ey_'+str(0).zfill(4)+'*.bin'))
    list2 = get_fn_from_dir(os.path.join(dir2, 'sl*_Ey_'+str(0).zfill(4)+'*.bin'))
    list3 = get_fn_from_dir(os.path.join(dir3, 'sl*_Ey_'+str(0).zfill(4)+'*.bin'))

    for flist in [list1,list2,list3]:
        for fname in flist:
            fname = fname.strip('\n')
            # logger.debug('rm '+ fname)
            try:
                os.remove(fname)
            except Exception as e:
                logger.error(e)

def clean_wavefield(isrc, workdir, dir1 = os.path.join('STD','Output'), dir2 = os.path.join('RTM','Output'), dir3 = 'Output'):
    dir1 = os.path.join(workdir,dir1)
    dir2 = os.path.join(workdir,dir2)
    dir3 = os.path.join(workdir,dir3)
    logger.info("Begin cleaning wavefield of src%d"%isrc) 
    
    list1 = get_fn_from_dir(os.path.join(dir1, 'wvf_Ey_'+str(0).zfill(4)+'*.bin'))
    list2 = get_fn_from_dir(os.path.join(dir2, 'wvf_Ey_'+str(0).zfill(4)+'*.bin'))
    list3 = get_fn_from_dir(os.path.join(dir3, 'wvf_Ey_'+str(0).zfill(4)+'*.bin'))

    for flist in [list1,list2,list3]:
        for fname in flist:
            fname = fname.strip('\n')
            # logger.debug('rm '+ fname)
            try:
                os.remove(fname)
            except Exception as e:
                logger.error(e)


def check_file(isrc,result_dir='Result'):
    logger.info("Checking files of src%d"%isrc)
    dum_s = True
    dum_w = True
    dum_s = os.path.isfile(os.path.join(result_dir, 'result_xcorr.dat'+'_'+str(isrc).zfill(4))) and dum_s
    dum_s = os.path.isfile(os.path.join(result_dir, 'result_ycorr.dat'+'_'+str(isrc).zfill(4))) and dum_s
    dum_s = os.path.isfile(os.path.join(result_dir, 'result_zcorr.dat'+'_'+str(isrc).zfill(4))) and dum_s
    dum_s = os.path.isfile(os.path.join(result_dir, 'result_xcorr_normal_forward.dat'+'_'+str(isrc).zfill(4))) and dum_s
    dum_s = os.path.isfile(os.path.join(result_dir, 'result_ycorr_normal_forward.dat'+'_'+str(isrc).zfill(4))) and dum_s
    dum_s = os.path.isfile(os.path.join(result_dir, 'result_zcorr_normal_forward.dat'+'_'+str(isrc).zfill(4))) and dum_s
    dum_s = os.path.isfile(os.path.join(result_dir, 'result_xcorr_normal_backward.dat'+'_'+str(isrc).zfill(4))) and dum_s
    dum_s = os.path.isfile(os.path.join(result_dir, 'result_ycorr_normal_backward.dat'+'_'+str(isrc).zfill(4))) and dum_s
    dum_s = os.path.isfile(os.path.join(result_dir, 'result_zcorr_normal_backward.dat'+'_'+str(isrc).zfill(4))) and dum_s
    dum_w = os.path.isfile(os.path.join(result_dir, 'result_wavefield_corr.dat'+'_'+str(isrc).zfill(4))) and dum_w
    dum_w = os.path.isfile(os.path.join(result_dir, 'result_wavefield_forward.dat'+'_'+str(isrc).zfill(4))) and dum_w
    dum_w = os.path.isfile(os.path.join(result_dir, 'result_wavefield_backward.dat'+'_'+str(isrc).zfill(4))) and dum_w
    return dum_s,dum_w




if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "fd:", ["force","workdir="])
    except getopt.GetoptError as err:
        # print help information and exit:
        logger.error(err)  # will print something like "option -a not recognized"
        # usage()
        sys.exit(2)

    isrc = int(args[0])
    logger.info("Begin cleaning files of src%d"%isrc)

    check = True
    workdir = ''
    for o, a in opts:
        if o in ('-f','--force'):
            logger.info('Zero-offset mode: clean without check.')
            check = False
        elif o in ('-d','--workdir'):
            workdir = a
        else:
            assert False, "unhandled option"


    if check:
        check_s,check_w = check_file(isrc)
        if check_s == True:
            clean_slice(isrc,workdir)
        else:
            logger.warning("src%d: No slices deleted, as results are not complete!"%isrc)
        if check_w == True:
            clean_wavefield(isrc,workdir)
        else:
            logger.warning("src%d: No wavefield deleted, as results are not complete!"%isrc)
    else:
        clean_slice(isrc,workdir)
        clean_wavefield(isrc,workdir)

    logger.info('src%d done.'%isrc)