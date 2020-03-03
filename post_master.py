# task_worker.py

import time, sys, queue, os, time, logging, datetime
from multiprocessing.managers import BaseManager
# from multiprocessing import Pool

#logger
today = datetime.date.today()
logger = logging.getLogger('post')
logger.setLevel(logging.DEBUG) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
fh = logging.FileHandler(os.path.join('log','post-' + today.strftime('%Y%m%d') + '.log'))
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s(%(process)d-%(processName)s): (%(levelname)s) %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)


# 发送任务的队列:
task_queue = queue.Queue()
# 接收结果的队列:
result_queue = queue.Queue()

# 创建类似的QueueManager:
class QueueManager(BaseManager):
    pass

# 把两个Queue都注册到网络上, callable参数关联了Queue对象:
QueueManager.register('get_task_queue', callable=lambda: task_queue)
QueueManager.register('get_result_queue', callable=lambda: result_queue)
# 绑定端口5000, 设置验证码'abc':
manager = QueueManager(address=('', 5555), authkey=b'rtm')



if __name__ == '__main__':
    cwd = os.getcwd()

    pypath = cwd
    py1 = os.path.join(pypath,'corr_RTM_wavefield_sub.py')
    py2 = os.path.join(pypath,'corr_RTM_slice_sub.py')
    py3 = os.path.join(pypath,'clean.py')
    
    logger.info('Current folder: %s'%cwd)
    # 启动Queue:
    manager.start()
    # 获得通过网络访问的Queue对象:
    task = manager.get_task_queue()
    result = manager.get_result_queue()

############## jobs here ##############
    while 1:
        try:
            [taskname,isrc,is_zRTM] = task.get(timeout=1800)
        except Exception as e:
            print('(%s)Listening... %s'%(time.strftime("%H:%M:%S", time.localtime(time.time())),e))
            continue

        if taskname.lower().strip() == 'end':
            logger.warning("'end' instruction detected, closing manager.")
            if task.empty():
                break
            else:
                logger.warning("Queue NOT empty, manager didn't close.")
                continue

        task_str = '%s-src%d'%(taskname,isrc)
        logger.info('Processing %s'%task_str)
        start_time = time.time()

        workpath = os.path.join(cwd,'tasks',taskname)
        if is_zRTM:
            p_status = os.system('cd %s;python %s -f %d'%(workpath,py3,isrc))
        else:
            p_status = os.system('cd %s;python %s %d;python %s %d;python %s %d'%(workpath,py1,isrc,py2,isrc,py3,isrc))
        logger.info('os.system status: %s'%str(p_status))

        end_time = time.time()
        logger.info('Job %s done (time cost: %.2fs).'%(task_str, end_time - start_time))
#######################################

    # 关闭:
    manager.shutdown()
    logger.info('master exit.')