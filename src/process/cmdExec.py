import subprocess

def execute(command: list, silent=False):
    """
    :param command: Command suitable for running in subprocess, must use a ['ls', '-l'] format
    :param silent: Run silently
    :return: Response from command
    """

    print(f'Running command: {command}')
    if silent:
        results = subprocess.run(command, shell=False, encoding='utf-8', check=True, stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        results = subprocess.run(command, shell=False, encoding='utf-8', check=True)

    if res.returncode != 0:
        print(f'Error running command: {" ".join(command)}. The error message was:\n{res.stderr}')
        exit(1)

    return results
