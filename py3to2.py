#!/usr/bin/env python3


import sys, optparse, json, os, glob, re, random, string
from subprocess import Popen, PIPE
from math import *


pattern1 = re.compile("print\(")
pattern2 = re.compile("str\.maketrans")
pattern3 = re.compile("python3")

outFile = open(sys.argv[2], 'w') 

with open(sys.argv[1], 'r') as inFile:

    for line in inFile:

        if re.search(pattern1, line):
            #print(line)        
            tmp_1 = re.sub(pattern1, "print ", line)
            tmp_2 = re.sub("\)$", " ", tmp_1)
            outFile.write(tmp_2)
            continue

        if re.search(pattern2, line):
            tmp_1 = re.sub(pattern2, "string.maketrans", line)
            outFile.write(tmp_1)
            continue

        if re.search(pattern3, line):
            tmp_1 = re.sub(pattern3, "python", line)
            outFile.write(tmp_1)
            continue

        outFile.write(line)

outFile.close()

