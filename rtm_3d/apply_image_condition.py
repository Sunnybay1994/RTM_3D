#!/usr/bin/env python
import os,re,argparse
from image_condition.get_filename_list import *
from image_condition.corr_RTM_wavefield_sub import *
from image_condition.corr_RTM_slice_sub import *
from image_condition.clean import *

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
        # get files
        path0,listx0,listy0,listz0,listw0 = get_isrc_filenames_data(isrc,workdir)
        path1,path2,listx1,listx2,listy1,listy2,listz1,listz2,listw1,listw2 = get_isrc_filenames_rtm(isrc,workdir)
        if args.mode == 'z':
            nocheck = True
        else:
            nocheck = False
            if not args.nocorr:
                # corr slice
                corr_RTM_slice_sub(isrc,workdir,path1,path2,listx1,listx2,listy1,listy2,listz1,listz2)
                # corr wavefield
                corr_RTM_wavefield_sub(isrc,workdir,path1,path2,listw1,listw2)
        # clean
        if not args.noclean:
            list0 = [os.path.join(path0,fn) for fn in listx0+listy0+listz0]
            listw0 = [os.path.join(path0,fn) for fn in listw0]
            list1 = [os.path.join(path1,fn) for fn in listx1+listy1+listz1]
            listw1 = [os.path.join(path1,fn) for fn in listw1]
            list2 = [os.path.join(path2,fn) for fn in listx2+listy2+listz2]
            listw2 = [os.path.join(path2,fn) for fn in listw2]
            clean(isrc,workdir,nocheck,list0+list1+list2,listw0+listw1+listw2)