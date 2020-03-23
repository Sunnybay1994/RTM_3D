import time, sys, queue, getopt, os, logging, datetime
from multiprocessing.managers import BaseManager


#logger
today = datetime.date.today()
logger = logging.getLogger('post')
logger.setLevel(logging.DEBUG) #CRITICAL>ERROR>WARNING>INFO>DEBUG》NOTSET
fh = logging.FileHandler(os.path.join('log','post_put-' + today.strftime('%Y%m%d') + '.log'))
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
logger.info('Connect to server %s...' % server_addr)
# 端口和验证码注意保持与task_master.py设置的完全一致:
m = QueueManager(address=(server_addr, 5678), authkey=b'rtm')
# 从网络连接:
m.connect()
# 获取Queue的对象:
task = m.get_task_queue()
result = m.get_result_queue()
if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "t:z", ["taskname=",'--zero-offset'])
    except getopt.GetoptError as err:
        # print help information and exit:
        logger.error(err)  # will print something like "option -a not recognized"
        # usage()
        sys.exit(2)

    is_zRTM = False
    for o, a in opts:
        if o in ('-t','--task'):
            taskname = a
        elif o in ('-z','--zero-offset'):
            logger.info('Zero-offset Mode.')
            is_zRTM = True
        else:
            assert False, "unhandled option"

    if len(args) == 1:
        isrcs = [int(args[0])]
    elif len(args) == 2:
        isrc_begin = int(args[0])
        isrc_end = int(args[1])
        isrcs = range(isrc_begin,isrc_end+1)
    else:
        isrcs = [int(i) for i in args]

    for isrc in isrcs:
        task.put((taskname,isrc,is_zRTM))
        task_str = '%s-src%d'%(taskname,isrc)
        logger.info('Put task: %s. (Tasks in queue: %d)'%(task_str,task.qsize()))