import numpy as np

def extend_source(srcinfo,srcpulse,xhalfspan=3,yhalfspan=3):
    #srcinfo: [x,y,z,component]
    x = srcinfo[0]
    y = srcinfo[1]
    z = srcinfo[2]
    others = tuple(srcinfo[3:])
    for i in range(-xhalfspan+1,xhalfspan):
        for j in range(-yhalfspan+1,yhalfspan):
            yield (x+i,y+j,z)+others,srcpulse

def write_source(fn,srcinfos,srcpulses):
    nsrc = len(srcinfos)
    assert len(srcinfos) == len(srcpulses), 'write_source:number of srcinfos(%d) and srcpulses(%d) not equal.'%(len(srcinfos),len(srcpulses))
    temp_bool = True
    nt_src = len(srcpulses[0])
    assert (temp_bool and nt_src==len(srcpulse) for srcpulse in srcpulses), 'length of each srcpulses not equal.'
    with open(fn,'w') as fo:
        fo.write("%d %d\n" % (nsrc,nt_src))
        list(map(lambda srcpos:fo.write("%d %d %d %s\n" %
                       (srcpos[0], srcpos[1], srcpos[2], srcpos[3])),srcinfos))
        np.savetxt(fo,srcpulses,'%g')

def extend_and_write_one_source(fn,srcinfo,srcpulse,xhalfspan=3,yhalfspan=3):
    srcinfos,srcpulses = zip(*extend_source(srcinfo,srcpulse))
    write_source(fn,srcinfos,srcpulses)

def extend_and_write_sources(fn,srcinfos,srcpulses,xhalfspan=3,yhalfspan=3):
    assert len(srcinfos) == len(srcpulses), 'extend_and_write_sources:number of srcinfos(%d) and srcpulses(%d) not equal.'%(len(srcinfos),len(srcpulses))
    new_srcinfos = []
    new_srcpulses = []
    for i in range(len(srcinfos)):
        ext_srcinfos,ext_srcpulses = zip(*extend_source(srcinfos[i],srcpulses[i]))
        new_srcinfos += ext_srcinfos
        new_srcpulses += ext_srcpulses
    write_source(fn,new_srcinfos,new_srcpulses)