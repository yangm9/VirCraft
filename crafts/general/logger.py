import logging
from datetime import datetime
from functools import wraps
import sys

version = '0.0.17'
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
        
        # create a dedicated recorder to avoid side effects on the root logger
        logger = logging.getLogger('VirCraft')
        logger.setLevel(self.levelDict[self.level])
        logger.handlers = []  # Clear the existing handlers
        logger.propagate = False
        
        # Configure standard output processor (DEBUG & INFO)
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)  # Treat all levels higher than DEBUG
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
        @wraps(func)
        def wrapper(*args, **kwargs):
            start_time = curr_time()
            self.logger.info(f'VirCraft {version} -- A flexible pipeline for metaviromic data analysis')
            self.logger.info(f'VirCraft {func.__name__} start ...')
            
            # Capture the args parameter
            if args:
                arg_obj = args[0]
                # Check if there is an unrun attribute and it is True
                if hasattr(arg_obj, 'unrun') and getattr(arg_obj, 'unrun'):
                    self.logger.info(f'VirCraft {func.__name__} will only output the shell scripts without execution (-u: unrun flag)')
            
            try:
                result = func(*args, **kwargs)
                self.logger.info(f'VirCraft {func.__name__} done!')
            except Exception as e:
                self.logger.error(f'{func.__name__} failed: {e}', exc_info=True)
                raise

            end_time = curr_time()
            elapsed_time = end_time - start_time
            self.logger.info(f'Elapsed time: {elapsed_time}')
            return result
        return wrapper
