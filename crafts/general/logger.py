import logging
from datetime import datetime
import sys

def curr_time():
    return datetime.now().replace(microsecond=0)

class Log:
    levelDict = {
        'DEBUG': logging.DEBUG,
        'INFO': logging.INFO,
        'WARNING': logging.WARNING,
        'ERROR': logging.ERROR,
        'CRITICAL': logging.CRITICAL
    }

    def __init__(self, level='INFO', logfile=None):
        self.level = level.upper()
        self.logfile = logfile

        # 检查日志级别是否有效
        if self.level not in self.levelDict:
            raise ValueError(f"Invalid log level: {self.level}. Choose from {list(self.levelDict.keys())}")

        # 创建日志记录器
        logger = logging.getLogger()
        logger.setLevel(self.levelDict[self.level])
        logger.handlers = []  # 清空已有的处理器

        # 配置标准输出处理器 (DEBUG 和 INFO)
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)  # 处理所有 DEBUG 及以上级别
        stdout_filter = logging.Filter()
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)  # 只处理 DEBUG 和 INFO
        stdout_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

        # 配置标准错误处理器 (WARNING 及以上)
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)  # 处理 WARNING 及以上级别
        stderr_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

        # 将处理器添加到记录器
        logger.addHandler(stdout_handler)
        logger.addHandler(stderr_handler)

        # if have a log file
        if self.logfile:
            file_handler = logging.FileHandler(self.logfile)
            file_handler.setLevel(self.levelDict[self.level])
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)

        self.logger = logger

    def __call__(self, func):
        def wrapper(*args, **kwargs):
            start_time = curr_time()
            self.logger.info(f"VirCraft {func.__name__} started at {start_time}")
            
            result = func(*args, **kwargs)
            
            end_time = curr_time()
            elapsed_time = end_time - start_time
            self.logger.info(f"VirCraft {func.__name__} ended at {end_time}. Elapsed time: {elapsed_time}")
            return result
        return wrapper
