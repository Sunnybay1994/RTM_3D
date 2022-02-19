import argparse,os,shutil
import numpy as np
from numpy.lib.function_base import interp

def resampleGather(indir,outdir,dt,dtout,ntr):
    if not os.path.isdir(outdir):
        print('Making dir: %s'%outdir)
        os.mkdir(outdir)
    for file in os.listdir(indir):
        if 'merge_gather' in file: # and '0000' in file:
            if 'loc' in file:
                print('Copying %s to %s'%(os.path.join(indir,file),os.path.join(outdir,file)))
                shutil.copy(os.path.join(indir,file),os.path.join(outdir,file))
            elif file.endswith('.bin'):
                print('Resampling %s'%os.path.join(indir,file))
                isum = np.fromfile(os.path.join(indir,file),dtype='float32')
                isum = isum.reshape(ntr,-1)
                nt = isum.shape[1]
                t = np.array(range(nt)) * dt
                t2 = np.arange(0,t[-1]+dtout,dtout)
                isum2 = np.zeros((ntr,len(t2)))
                for i in range(ntr):
                    isum2[i,:] = interp(t2,t,isum[i,:])
                print('Saving to %s'%os.path.join(outdir,file))
                isum2.astype('float32').tofile(os.path.join(outdir,file))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Resample the time-axis of the merged gathers and save to another task.',conflict_handler='resolve')
    parser.add_argument('indir',type=str,help="Input task-directory.")
    parser.add_argument('-o','--outdir',type=str, required=True,help="Output task-directory.")
    parser.add_argument('--dt',type=float, required=True,help="Time interval of input data.")
    parser.add_argument('--dtout',type=float, required=True,help="Time interval of output data.")
    parser.add_argument('--ntr',type=int, required=True,help="Number of traces of input data.")

    args = parser.parse_args()

    indir = args.indir
    outdir = args.outdir
    dt = args.dt
    dtout = args.dtout
    ntr = args.ntr

    resampleGather(os.path.join(indir,'Output'),os.path.join(outdir,'Output'),dt,dtout,ntr)
    resampleGather(os.path.join(indir,'STD','Output'),os.path.join(outdir,'STD','Output','ext'),dt,dtout,ntr)