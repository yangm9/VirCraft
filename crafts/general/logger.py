import logging
from datetime import datetime




def curr_time():
    return datetime.now().replace(microsecond=0)

def calc_time(func):
    def wrapper(*args, **kwargs):
        start_time = curr_time()
        print(f"{start_time} - VirCraft")
        result=func(*args,**kwargs)
        end_time = curr_time()
        print(f"{end_time} - VirCraft end!")
        elapsed_time = end_time - start_time
        print(f"Total time taken: {elapsed_time}")
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
    def __init__(self,level=None,logfile=None):
        self.level=level
        self.logfile=logfile
        logging.basicConfig(filename=self.logfile,level=self.levelDict[level])
    def __call__(self, func):
        def wrapper(*args, **kwargs):
            start_time=curr_time()
            logging.info(f"{start_time} - VirCraft {func.__name__} start")
            result=func(*args,**kwargs)
            end_time=curr_time()
            logging.info(f"{end_time} - VirCraft {func.__name__} end")
            elapsed_time = end_time - start_time
            logging.info(f"VirCraft {func.__name__} elapsed time: {elapsed_time}")
            return result
        return wrapper
