import subprocess

def execute(command:list, silent=False):
    '''
    :ref: virmatcher
    :param command: Command suitable for running in subprocess, must use a ['ls', '-l'] format
    :param silent: Run silently
    :return: Response from command
    '''
    
    cmd_txt=' '.join(command).replace('\n ','\n')
    cmd_txt=cmd_txt.replace(' \n ','\n')
    print(f'Running command:\n{cmd_txt}')
    '''
    if silent:
        results = subprocess.run(command, shell=False, encoding='utf-8', check=True, stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    else:
        results = subprocess.run(command, shell=False, encoding='utf-8', check=True)

    if results.returncode != 0:
        print(f'Error running command: {" ".join(command)}. The error message was:\n{res.stderr}')
        exit(1)
    '''
    return '0'#results
