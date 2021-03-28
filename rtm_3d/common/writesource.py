import numpy as np

def extend_source(srcinfo,srcpulse,xhalfspan=2,yhalfspan=2):
    #srcinfo: [x,y,z,component]
    x = srcinfo[0]
    y = srcinfo[1]
    z = srcinfo[2]
    others = tuple(srcinfo[3:])
    for i in range(-xhalfspan,xhalfspan+1):
        for j in range(-yhalfspan ,yhalfspan+1):
            yield (x+i,y+j,z)+others,srcpulse

def write_source(fn,srcinfos,srcpulses,nsrc_ori,nspan):
    nsrc = len(srcinfos)
    assert nsrc == nsrc_ori * nspan, 'Total num of src not correct %d!=%dx%d'%(nsrc, nsrc_ori, nspan)
    assert nsrc_ori == len(srcpulses), 'write_source:number of srcinfos(%d) and srcpulses(%d) not equal.'%(len(srcinfos),len(srcpulses))
    temp_bool = True
    nt_src = len(srcpulses[0])
    assert (temp_bool and nt_src==len(srcpulse) for srcpulse in srcpulses), 'length of each srcpulses not equal.'
    with open(fn,'w') as fo:
        fo.write("%d,%d,%d\n" % (nsrc_ori,nspan,nt_src))
        list(map(lambda srcpos:fo.write("%d,%d,%d,%s\n" %
                       (srcpos[0], srcpos[1], srcpos[2], srcpos[3])),srcinfos))
        np.savetxt(fo,srcpulses,'%g')

def extend_and_write_one_source(fn,srcinfo,srcpulse,xhalfspan=2,yhalfspan=2):
    srcinfos,srcpulses = zip(*extend_source(srcinfo,srcpulse,xhalfspan,yhalfspan))
    nspan = len(srcinfos)
    write_source(fn,srcinfos,[srcpulse],1,nspan)

def extend_and_write_sources(fn,srcinfos,srcpulses,xhalfspan=2,yhalfspan=2,autospan=False):
    status_info = 'extend_and_write_sources: number of original sources=%d\n'%len(srcinfos)
    if autospan:
        status_info += 'Autospan: '
        xs,ys,*_ = zip(*srcinfos)
        xs = list({}.fromkeys(xs).keys()) # remove repeat elements
        xs.sort() # sort
        xhalfspan = (xs[1] - xs[0])//2 # get difference
        ys = list({}.fromkeys(ys).keys())
        ys.sort()
        yhalfspan = (ys[1] - ys[0])//2
    assert len(srcinfos) == len(srcpulses), 'extend_and_write_sources:number of srcinfos(%d) and srcpulses(%d) not equal.'%(len(srcinfos),len(srcpulses))
    new_srcinfos = []
    new_srcpulses = []
    for i in range(len(srcinfos)):
        ext_srcinfos,ext_srcpulses = zip(*extend_source(srcinfos[i],srcpulses[i],xhalfspan,yhalfspan))
        new_srcinfos += ext_srcinfos
        new_srcpulses += ext_srcpulses
    nsrc = len(srcinfos)
    nspan = len(ext_srcinfos)
    status_info += 'xhalfspan=%d, yhalfspan=%d\n'%(xhalfspan,yhalfspan)
    status_info += 'total_span/source: %d, total sources after span: %d\n'%(nspan,len(new_srcinfos))
    write_source(fn,new_srcinfos,srcpulses,nsrc,nspan)
    return status_info