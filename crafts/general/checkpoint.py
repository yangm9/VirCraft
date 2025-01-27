class CheckPoint:
    '''
    '''
    def __init__(self, filename = 'Done'):
        self.filename = filename
    def __call__(self, func):
        def wrapper(*args, **kwargs):
            result = func(*args, **kwargs)
            with open(self.filename, "w") as file:
                pass
            return result
        return wrapper
