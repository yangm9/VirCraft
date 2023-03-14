#!/usr/bin/env python3

import sys
import re

inputf=sys.argv[1]
with open(inputf,'r') as f:
    for i in f.readlines():
        if i.startswith('>'):
            ID=i.strip('>').split(' ')[0]
            sp=i.split('[',1)[1].rsplit(']',1)[0]
            if re.search(r'\].*\[',sp):
                sp=sp.split('[')[1]
            print(ID+'\t'+sp)
