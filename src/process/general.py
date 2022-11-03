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

def insLable(file_name: str, label: str):
    "Insert a label before the extension for a file name."
    prefix, extension = os.path.splitext(file_name)
    return f'{prefix}.{label}{extension}
