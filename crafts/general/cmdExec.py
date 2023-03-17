import subprocess

def execute(cmd:list, silent=False):
    '''
    :ref: virmatcher, https://github.com/bolduc/kb_virmatcher
    :param command: Command suitable for running in subprocess, must use a ['ls', '-l'] format
    :param silent: Run silently
    :return: Response from command
    '''
    print(f'Running command:\n{cmd_txt}') 
    cmd_txt=' '.join(cmd).replace('\n ','\n')
    cmd_txt=cmd_txt.replace(' \n ','\n')

    if silent:
        results=subprocess.run(cmd,shell=False,encoding='utf-8',check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    else:
        results=subprocess.run(cmd,shell=False,encoding='utf-8',check=True)
    
    if results.returncode != 0:
        print(f'Error running command: {cmd_txt}. The error message was:\n{results.stderr}')
        exit(1)

    return results

def show_cmd(cmd):
    cmd_txt=' '.join(cmd).replace('\n ','\n')
    cmd_txt=cmd_txt.replace(' \n ','\n')
    print(f'Pending command:\n{cmd_txt}')
    return 0 
