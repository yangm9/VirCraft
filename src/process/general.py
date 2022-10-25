import os

def mkdir(name):
    if not os.path.exists(name):
        os.makedirs(name)
    return 0
