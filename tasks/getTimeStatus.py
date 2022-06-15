import datetime,os

def DurationTotal(logdir,nsrc=-1):
    runtimes = []
    if nsrc > 0:
        pass
    else:
        for file in os.listdir(logdir):
            if file.startswith('script_0') and file.endswith('.out'):
                runtimes.append(Duration(os.path.join(logdir,file)))
    return sum(runtimes),runtimes

def Duration(logfn):
    starttext = '=====Computing started at '
    endtext = '=====Computing stopped at '
    for line in open(logfn):
        if starttext in line:
            rawtime = line.replace(starttext,'').strip('=\n').split()
            rawtime.pop(4)
            rawtime.pop(0)
            starttime = datetime.datetime.strptime('-'.join(rawtime),'%b-%d-%X-%Y')
        elif endtext in line:
            rawtime = line.replace(endtext,'').strip('=\n').split()
            rawtime.pop(4)
            rawtime.pop(0)
            endtime = datetime.datetime.strptime('-'.join(rawtime),'%b-%d-%X-%Y')
        else:
            continue
    return (endtime-starttime).total_seconds()

if __name__ == "__main__":
    print(DurationTotal(os.path.join('PDM','npdm2_g0.05m_400MHz_0.5x0.5_2_0.1x0.1_pstd_8','log','p5')))
