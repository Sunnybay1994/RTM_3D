import time, sys, queue, os, time, logging, datetime
from multiprocessing.managers import BaseManager
# from multiprocessing import Pool

#logger
today = datetime.date.today()
logger = logging.getLogger('post')
logger.setLevel(logging.DEBUG) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
fh = logging.FileHandler(os.path.join('log','post_worker-' + today.strftime('%Y%m%d') + '.log'))
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s(%(process)d-%(processName)s): (%(levelname)s) %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)

# 创建类似的QueueManager:
class QueueManager(BaseManager):
    pass
# 由于这个QueueManager只从网络上获取Queue，所以注册时只提供名字:
QueueManager.register('get_task_queue')
QueueManager.register('get_result_queue')

# 连接到服务器，也就是运行task_master.py的机器:
server_addr = 'pkugem.synology.me'
print('Connect to server %s...' % server_addr)
# 端口和验证码注意保持与task_master.py设置的完全一致:
m = QueueManager(address=(server_addr, 5678), authkey=b'rtm')
# 从网络连接:
m.connect()
# 获取Queue的对象:
task = m.get_task_queue()
result = m.get_result_queue()

if __name__ == '__main__':
    cwd = os.getcwd()

    pypath = cwd
    py1 = os.path.join(pypath,'corr_RTM_wavefield_sub.py')
    py2 = os.path.join(pypath,'corr_RTM_slice_sub.py')
    py3 = os.path.join(pypath,'clean.py')
    
    logger.info('Current folder: %s'%cwd)

############## jobs here ##############
    while 1:
        try:
            (taskname,isrc,is_zRTM) = task.get(timeout=1800)
            # end master
            if taskname.lower().strip() == 'end_master':
                result.put((taskname,isrc,is_zRTM,-1))
                continue
        except ConnectionError as e:
            logger.error('ConnectionError, worker will exit.(Error info: %s)'%e)
            break
        except Exception as e:
            print('(%s)Waiting for job... %s'%(time.strftime("%H:%M:%S", time.localtime(time.time())),e))
        else: 
            task_str = '%s-src%d'%(taskname,isrc)
            logger.info('Processing %s (Tasks waiting: %d)'%(task_str,task.qsize()))
            result.put((taskname,isrc,is_zRTM,-1)) # -1 represents jobs are processing.

            workpath = os.path.join(cwd,'tasks',taskname)
            start_time = time.time()
            if is_zRTM:
                p_status = os.system('cd %s;python %s -f %d'%(workpath,py3,isrc))
            else:
                p_status = os.system('cd %s;python %s %d;python %s %d;python %s %d'%(workpath,py1,isrc,py2,isrc,py3,isrc))
            end_time = time.time()
            logger.info('Job %s done, status code: %d, time cost: %.2fs.\n'%(task_str, p_status, end_time - start_time))

            result.put((taskname,isrc,is_zRTM,p_status)) # p_status==0 means succeed.

            
#######################################