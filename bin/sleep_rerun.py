#!/usr/bin/env
import time
import task_states as ts

def judge(dvf_sh):
    states=ts.check_sh_stat(dvf_sh)
    return ts.check_process_states(states)

def rerun(dvf_sh):
    try:
        judgement=judge(dvf_sh)
    if judgement:
        time.sleep(600)
    else:
        time.sleep(600)

