import logging
from datetime import datetime

def curr_time():
    return f'[{str(datetime.now().replace(microsecond=0))}]'

def calc_time(func):
    def wrapper(*args, **kwargs):
        start_time=curr_time()
        print(start_time)
        result=func(*args,**kwargs)
        end_time=curr_time()
        print(end_time)
        return result
    return wrapper

class Log:
    levelDict={
        'DEBUG':logging.DEBUG,
        'INFO':logging.INFO,
        'WARNING':logging.WARNING,
        'ERROR':logging.ERROR,
        'CRITICAL':logging.CRITICAL
    }
    def __init__(self,level='',logfile=''):
        self.level=level
        self.logfile=logfile
        logging.basicConfig(filename=self.logfile,level=self.levelDict[level])
    def __call__(self, func):
        def wrapper(*args, **kwargs):
            start_time=curr_time()
            logging.info(f"{start_time} {func.__name__}() start")
            result=func(*args,**kwargs)
            end_time=curr_time()
            logging.info(f"{end_time} {func.__name__}() end")
            return result
        return wrapper
