#!/usr/bin/env
import time
import task_states as ts

def rerun(dvf_sh):
    stat=ts.judge(dvf_sh)
    while True:
        time.sleep(300)
        if not stat:
            time.sleep(5)
            if not judge(script_path):

    return 0

