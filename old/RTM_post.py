#!/usr/bin/env python
import sys,os,logging
import multiprocessing as mp

from corr_RTM_slice_sub import corr_slice
from corr_RTM_wavefield_sub import corr_wavefield
from clean import *

#logger
logger = logging.getLogger('RTM_post')
logger.setLevel(logging.DEBUG) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
fh = logging.FileHandler('log/RTM_post.log')
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s(%(process)d-%(processName)s): (%(levelname)s) %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)

def corr_rtm(isrc):
    logger.info('rtm post: src'+str(isrc))

    check_s,check_w = check_file(isrc)
    if not check_s:
        logger.info('corr_slice src%d'%isrc)
        status_s = corr_slice(isrc)
        if status_s:
            logger.info("Clean slice of src%d"%isrc) 
            clean_slice(isrc)
        else:
            logger.warning('src%d: corr_slice failed.'%isrc)
    else:
        logger.info("Clean slice of src%d"%isrc) 
        clean_slice(isrc)

    if not check_w:
        logger.info('corr_wavefield src%d'%isrc)
        status_w = corr_wavefield(isrc)
        if status_w:
            logger.info("Clean wavefield of src%d"%isrc) 
            clean_wavefield(isrc)
        else:
            logger.warning('src%d: corr_wavefield failed.'%isrc)
    else:
        logger.info("Clean wavefield of src%d"%isrc) 
        clean_wavefield(isrc)
    logger.info('rtm post done: src'+str(isrc))

if __name__ == "__main__":
    idir = 'Result'
    if not os.path.exists(idir):
        os.mkdir(idir)
    src_start = int(sys.argv[1])
    src_end = int(sys.argv[2])

    # p = mp.Pool()
    # p.map(corr_rtm, range(src_start,src_end+1))
    for isrc in range(src_start,src_end+1):
        corr_rtm(isrc)