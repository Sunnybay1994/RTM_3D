#!/usr/bin/env python
import os,re
from common import *

result_slice_fn = 'result_{xyz}corr_{isrc:0>4d}_{islice:0>2d}.{ext}'
result_slice_nf_fn = 'result_{xyz}corr_norm_forward_{isrc:0>4d}_{islice:0>2d}.{ext}'
result_slice_nb_fn = 'result_{xyz}corr_norm_backward_{isrc:0>4d}_{islice:0>2d}.{ext}'
result_wavefield_fn = 'result_wavefield_corr_{isrc:0>4d}.{ext}'
result_wavefield_f_fn = 'result_wavefield_forward_corr_{isrc:0>4d}.{ext}'
result_wavefield_b_fn = 'result_wavefield_backward_corr_{isrc:0>4d}.{ext}'

def gen_patterns(isrc):
    pattern_slx = re.compile(r'slx_Ey_%04d_[0-9]{5}\.bin'%isrc)
    pattern_sly = re.compile(r'sly_Ey_%04d_[0-9]{5}\.bin'%isrc)
    pattern_slz = re.compile(r'slz_Ey_%04d_[0-9]{5}\.bin'%isrc)
    pattern_wvf = re.compile(r'wvf_Ey_%04d_[0-9]{5}\.bin'%isrc)
    return pattern_slx,pattern_sly,pattern_slz,pattern_wvf

def get_isrc_filenames(isrc,wdir):
    path=os.path.join(wdir,'Output')
    listx,listy,listz,listw = filename_re_match(path,'file',*gen_patterns(isrc))
    return path,listx,listy,listz,listw

def get_isrc_filenames_data(isrc,wdir):
    path,listx,listy,listz,listw = get_isrc_filenames(isrc,wdir)
    return path,listx,listy,listz,listw

def get_isrc_filenames_rtm(isrc,wdir,dir1='STD',dir2='RTM'):
    dir1=os.path.join(wdir,dir1)
    dir2=os.path.join(wdir,dir2)
    path1,listx1,listy1,listz1,listw1 = get_isrc_filenames(isrc,dir1)
    path2,listx2,listy2,listz2,listw2 = get_isrc_filenames(isrc,dir2)
    assert len(listx1) == len(listx2), 'slicex number not equal (std%d-rtm%d)'
    assert len(listy1) == len(listy2), 'slicey number not equal'
    assert len(listz1) == len(listz2), 'slicez number not equal'
    assert len(listw1) == len(listw2), 'wavefield number not equal'
    listx1.sort()
    listx2.sort(reverse=True)
    listy1.sort()
    listy2.sort(reverse=True)
    listz1.sort()
    listz2.sort(reverse=True)
    listw1.sort()
    listw2.sort(reverse=True)
    return path1,path2,listx1,listx2,listy1,listy2,listz1,listz2,listw1,listw2
