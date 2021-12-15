#!/usr/bin/env python
import os,re,argparse
from image_condition.get_filename_list import *
from image_condition.corr_RTM_wavefield_sub import *
from image_condition.corr_RTM_slice_sub import *
from image_condition.clean import *
from multiprocessing import Process

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Apply cross-correlation imaging method for RTM and conduct the clean stuff.',conflict_handler='resolve')
    parser.add_argument('srcs', metavar='isrc', type=int, nargs='+',help='The src numbers of which shot gathers to be migrated.')
    parser.add_argument('-d','--dir',type=str,dest='workdir',default='.',help='The path of the work directory.')
    parser.add_argument('-m','--mode',choices=['m','z','mz'],default='m',help="Mode: 'm' for multi-offset, 'z' for zero-offset, 'mz' for both multi- and zero-offset.")
    parser.add_argument('--noclean',default=False,action='store_const',const=True,help='Donnot clean middle results.')
    parser.add_argument('--nocorr',default=False,action='store_const',const=True,help='Donnot conduct cross-correlation.')
    args = parser.parse_args()

    workdir = args.workdir
    for isrc in args.srcs:
        if ('m' in args.mode) and (not args.nocorr):
            path1,path2,listx1,listx2,listy1,listy2,listz1,listz2,listw1,listw2 = get_isrc_filenames_rtm(isrc,workdir)
            # corr slice
            corr_RTM_slice_sub(isrc,workdir,path1,path2,listx1,listx2,listy1,listy2,listz1,listz2)
            # p = Process(target=corr_RTM_slice_sub,args=(isrc,workdir,path1,path2,listx1,listx2,listy1,listy2,listz1,listz2))
            # p.start()
            # corr wavefield
            corr_RTM_wavefield_sub(isrc,workdir,path1,path2,listw1,listw2)
            # p.join()
        # clean
        if not args.noclean:
            clean(isrc,workdir,args.mode)