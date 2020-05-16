import numpy as np
import scipy.io as sio
import os

def trace_normal_moveout_layered(gather,tt,offset,v,z_v):
    ###################################
    # gather: data of one trace
    # tt: traveltime corresponding to gather
    # midpos: (x,y,z) position of the gather, midpoint of source and receiver
    # offset: distance between source and receiver
    # v: 1D velocity model at pos
    # z_v: z-coodrinate corresponding to v, should be above zero
    ###################################

    z_v = np.array(z_v)
    assert z_v.all()>0, 'z_v should be always > 0'
    v = np.array(v)
    # extend v and z_v if necessary
    v_max = np.max(v)
    z_max = tt[-1]*v_max
    if z_max > z_v[-1]:
        dz = z_v[-1] - z_v[-2]
        z_ext = np.arange(z_v[-1]+dz,z_max,dz)
        v_ext = v[-1]*np.ones(np.shape(z_ext))
        z_v = np.append(z_v,z_ext)
        v = np.append(v,v_ext)
    # calculate v_rms
    z_vv = np.insert(z_v,0,0)
    dz_v = z_vv[1:] - z_vv[:-1] # space step
    dt_v = dz_v/v # delta_t
    v_rms = np.zeros(np.shape(v))
    for i in range(len(v)):
        v_rms[i] = np.sqrt(np.sum(v[:i+1]**2*dt_v[:i+1])/np.sum(dt_v[:i+1]))
    t0 = dt_v.copy()
    for i in range(1,len(dt_v)):
        t0[i] += t0[i-1]
    # calculate t correspond to t0
    t = np.sqrt(t0**2 + (offset/v_rms)**2)
    # change tt to zero-offset time tt0
    gather = np.array(gather)
    tt = np.array(tt)
    tt1 = tt[tt>t[0]]
    gather1 = gather[tt>t[0]]
    tt0 = np.interp(tt1,t,t0)
    # interpolate gather from tt0 to the ori time series tt
    gather0 = np.interp(tt,tt0,gather1,0,0)
    return gather0

if __name__ == '__main__':
    ### load model ###
    try:
        dic_model = sio.loadmat(os.path.join('example','model','model.mat'))
        dict_sr = sio.loadmat(os.path.join('example','model','model_sr.mat'))
    except Exception as e:
        raise e
    else:
        pass
    finally:
        pass
    ### load model end ###

    # load src and velocity info
    epr = np.array(dic_model['ep'])
    dt = float(dic_model['dt'])
    dx = float(dic_model['dx'])
    dy = float(dic_model['dy'])
    dz = float(dic_model['dz'])
    nx = int(dic_model['nx'])
    ny = int(dic_model['ny'])
    nz_air = int(dic_model['nz_air'])
    nz = int(dic_model['nz'])

    v = 299792458/np.sqrt(epr)
    z = (np.array(range(nz)) - nz_air) * dz
    
    isrc = 60
    srcx = int(round(dict_sr['srcx'][0][isrc]/dx))
    srcy = int(round(dict_sr['srcy'][0][isrc]/dy))

    # load gather info
    with open(os.path.join('example','gather', 'merge_gather_'+str(isrc).zfill(4)+'.dat')) as fp:
        isum = np.loadtxt(fp)
    with open(os.path.join('example','gather', 'merge_gather_loc_'+str(isrc).zfill(4)+'.dat')) as fp:
        iloc = np.loadtxt(fp)

    # do NMO
    locs = []
    gathers = []
    for i in range(len(isum)):
        loc = [int(round((iloc[i][0]+srcx)/2)),int(round((iloc[i][1]+srcy)/2)),int(iloc[i][2])]
        locs += [loc]
        offset = np.sqrt(((iloc[i][0]-srcx)*dx)**2 + ((iloc[i][1]-srcy)*dy)**2)
        gather = trace_normal_moveout_layered(isum[i],np.array(range(len(isum[i])))*dt,offset,v[srcx][srcy][z>0],z[z>0])
        gathers += [gather]
    with open(os.path.join('example','gather', 'merge_gather_NMO_'+str(isrc).zfill(4)+'.dat'),'w') as fp:
        np.savetxt(fp,gathers)
    with open(os.path.join('example','gather', 'merge_gather_loc_NMO_'+str(isrc).zfill(4)+'.dat'),'w') as fp:
        np.savetxt(fp,locs,'%d')
