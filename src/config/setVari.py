#!/usr/bin/env python3
def selectENV(conda_sh,env_s):
    return f'. "{conda_sh}"\nconda activate {env_s}\n'
