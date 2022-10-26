import os

def mkdir(name: str):
    if not os.path.exists(name):
        os.makedirs(name)
    return 0

def printSH(sh_path: str, command: list):
    SHPATH = open(sh_path,'w')
    scripts = ' '.join(command).replace('\n ', '\n')
    scripts = scripts.replace(' \n ', '\n')
    SHPATH.write(scripts)
    SHPATH.close()
    return 0
    
