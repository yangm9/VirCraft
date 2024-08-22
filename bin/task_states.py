import os
import sys
import subprocess

#判断脚本是否在运行

def get_process_state(pid):
    with open(f'/proc/{pid}/status') as f:
        for line in f:
            if line.startswith('State:'):
                return line.split(':')[1].strip().split(' ')[0]

def get_children(pid):
    result = subprocess.run(f'pgrep -P {pid}', shell=True, stdout=subprocess.PIPE, text=True)
    return result.stdout.splitlines()

def check_state(pid, level=0):
    state = get_process_state(pid)
    #prefix = '    ' * level
    #print(f'{prefix}Process {pid}: {state}')
    result = f'{level}:{state}'
    children = get_children(pid)
    for child_pid in children:
        result += '\n' + check_state(child_pid, level+1)
    return result

def check_sh_stat(sh_name):
    pids = subprocess.run(f'pgrep -f {sh_name}', shell=True, stdout=subprocess.PIPE, text=True).stdout.splitlines()
    #for pid in pids:
    states=check_state(pids[0]).split('\n')
    return states

def check_process_stat(lists):
    result = []
    for process_states in lists:
        if process_states.endswith('S'):
            result.append(False)
        else:
            result.append(True)
    return any(result)

def judge(sh_name):
    states=check_sh_stat(sys.argv[1])
    return check_process_stat(states)

if __name__=='__main__':
    print(judge(sys.argv[1]))
