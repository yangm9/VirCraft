import logging
from datetime import datetime
import sys

def curr_time():
    return datetime.now().replace(microsecond=0)

class Log:
    levelDict = {
        'DEBUG' : logging.DEBUG,
        'INFO' : logging.INFO,
        'WARNING' : logging.WARNING,
        'ERROR' : logging.ERROR,
        'CRITICAL' : logging.CRITICAL
    }
    def __init__(self, level = 'INFO', logfile = None):
        self.level = level.upper()
        self.logfile = logfile
        # check if the log level is valid
        if self.level not in self.levelDict:
            raise ValueError(f"Invalid log level: {self.level}. Choose from {list(self.levelDict.keys())}")
        # create the log recorder
        logger = logging.getLogger()
        logger.setLevel(self.levelDict[self.level])
        logger.handlers = []  # Clear the existing processors
        # Configure standard output processor (DEBUG & INFO)
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)  # Treat all levels higher than DEBUG
        stdout_filter = logging.Filter()
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)  # Only treat DEBUG and INFO
        stdout_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        # Configure standard error processor (WARNING++)
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)  # Treat all levels higher than WARNING
        stderr_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        # Add the processor to the recorder
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
            self.logger.info(f"VirCraft {func.__name__} start ...")
            result = func(*args, **kwargs)
            end_time = curr_time()
            elapsed_time = end_time - start_time
            self.logger.info(f"VirCraft {func.__name__} done!\nElapsed time: {elapsed_time}")
            return result
        return wrapper
