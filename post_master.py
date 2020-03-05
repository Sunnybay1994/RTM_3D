import time, sys, queue, os, time, logging, datetime
from multiprocessing.managers import BaseManager
# from multiprocessing import Pool

#logger
today = datetime.date.today()
logger = logging.getLogger('post')
logger.setLevel(logging.DEBUG) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
fh = logging.FileHandler(os.path.join('log','post_master-' + today.strftime('%Y%m%d') + '.log'))
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
manager = QueueManager(address=('', 5678), authkey=b'rtm')



if __name__ == '__main__':
    cwd = os.getcwd()
    
    logger.info('Current folder: %s'%cwd)
    # 启动Queue:
    manager.start()
    # 获得通过网络访问的Queue对象:
    task = manager.get_task_queue()
    result = manager.get_result_queue()

############# master here #############
    jobs_processing = {} # processing list(dict): {[taskname,isrc,is_zRTM]:[p_status,start_time]}
    while 1:
        try:
            (taskname,isrc,is_zRTM,p_status) = result.get(timeout=1800)
            # print(taskname)
            # end master 
            if taskname.lower().strip() == 'end_master':
                logger.info("'end_master' instruction detected, closing manager.")
                if jobs_processing and not task.empty():
                    task.put((taskname,isrc,is_zRTM))
                    logger.warning("Queue NOT empty or still processing, waiting to close.")
                    continue
                else:
                    break
        except Exception as e:
            print('(%s)Waiting for results... %s'%(time.strftime("%H:%M:%S", time.localtime(time.time())),e))
        else:
            # useful variable
            task_str = '%s-src%d'%(taskname,isrc)
            time_now = time.time()
            # add or remove jobs from processing list
            if p_status == -1:
                # assign task
                logger.info('Task assigned: %s'%task_str)
                jobs_processing[(taskname,isrc,is_zRTM)] = (p_status,time_now)
            else:
                # task done: del item from processing list
                if p_status == 0:
                    # success
                    logger.info('Task complete: %s'%task_str)
                else:
                    # fail: put task back to task queue
                    logger.warning('Task not success: %s (return code: %d)'%(task_str,p_status))
                    task.put((taskname,isrc,is_zRTM))
                    logger.info('Put task back: %s'%task_str)
                try:
                    del(jobs_processing[(taskname,isrc,is_zRTM)])
                except Exception as e:
                    logger.warning('Error deleting task in processing list: %s'%e)
        finally:
            # put task back if timeout
            time_now = time.time()
            timeout = 3600*3 #3h
            for k in jobs_processing.keys():
                lasting_time = time_now - jobs_processing[k][1]
                if lasting_time > timeout:
                    task_str = '%s-src%d'%(k[0],k[1])
                    del(jobs_processing[k])
                    task.put(k)
                    logger.warning('task time out: %s (%.2fs)'%(task_str,lasting_time))


#######################################

    # 关闭:
    manager.shutdown()
    logger.info('master exit.')