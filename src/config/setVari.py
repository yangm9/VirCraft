#!/usr/bin/env python3

def selectENV(env: str):
    conda_path=isInstalled('conda')
    if conda_path:        
        return f'export PATH="{conda_path}:$PATH"\nconda activate {env}\n'
    else:
        return 'Error!!!'

def isInstalled(name):
    """
    Check whether `name` is on PATH and marked as executable.
    """
    from shutil import which
    return which(name)
