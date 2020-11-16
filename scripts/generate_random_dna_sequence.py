#!/usr/bin/env python

import sys
import random

bases = 'ACGT'
for _ in range(int(sys.argv[1])):
    print(random.choice(bases), end='')
