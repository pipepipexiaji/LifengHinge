#!/usr/bin/env python

import base
import readers
try:
    import numpy as np
except:
    import mynumpy as np
import sys
import os
import getopt


GROOVE_ENV_VAR = "OXDNA_GROOVE"
POS_BACK = -0.4
try:
    args, files = getopt.getopt(sys.argv[1:], "-c")
except:
    base.Logger.log("Wrong usage. Aborting", base.Logger.CRITICAL)
    sys.exit(-2)


if len(files) < 2:
    base.Logger.log("Usage is %s configuration topology" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

if len(files) > 2:
    output = files[2]
else:
    output = files[0] + ".pdb"

append = 'w'
l = readers.LorenzoReader(files[0], files[1])
s = l.get_system()
count = 1

while s:
    s._prepare(None)
    # This script only build a single frame
    result = "HEADER    frame t= " +str(s._time)  # + " \nMODEL     " + "1" + " \nREMARK ## 0,0\n"
    commands = []
    commands.append("set bg_color white")
    commands.append("~bond #0")
    sel_ref = ""
    color = ""
    for ss in s._strands:
        strid = ss.index + 1
        nid = 0
        result += "MODEL     " + str(strid) + " \n"
        atomoutput = ""
        n_start = ss._first -1
        sel_ref += str(n_start) + "\n"

        for nucleo in ss._nucleotides:
            nid += 1
            n_ind = nucleo.index - n_start
            # the s# holds the position vector of each nucleotide element
            s1 = nucleo.cm_pos_box + nucleo.get_pos_back_rel()
            """
            if os.environ.get(GROOVE_ENV_VAR) == '1':
                s2 = nucleo.cm_pos_box
                index_jump = 3
            else:
                index_jump = 2
            """
            index_jump = 2
            s3 = nucleo.cm_pos_box + (POS_BACK + 0.68) * nucleo._a1

            # print the backbone site
            atomoutput += "ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (
                index_jump * n_ind - 1, "A", "ALA", 'A', n_ind, ' ', s1[0], s1[1], s1[2], 1, 7.895)
            # print the centre of mass site (if grooving is on)
            #if os.environ.get(GROOVE_ENV_VAR) == '1':
                #atomoutput += "ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (
                    #index_jump * n_ind, "B", "ALA", 'B', n_ind, ' ', s2[0], s2[1], s2[2], 1,
                    #7.895)

            # print the base site
            atomoutput += "ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (
                index_jump * n_ind, 'N', "ALA", 'C', n_ind, ' ', s3[0], s3[1], s3[2], 1,
                6.316)

            # get command file output
            if os.environ.get(GROOVE_ENV_VAR) == '1':
                commands.append("bond #0.%i:%i.A:%i.B" % (strid, n_ind, n_ind))
                commands.append("bond #0.%i:%i.B:%i.C" % (strid, n_ind, n_ind))
            else:
                commands.append("bond #0.%i:%i" % (strid, n_ind))
                color += "\nbondcolor cyan #0.%i:%i" % (strid, n_ind)

            if n_ind != 1:
                commands.append("bond #0.%i:%i.A,%i.A" % (strid, n_ind - 1, n_ind))

        result += atomoutput + "TER \nENDMDL \n"
    commands.append("setattr m stickScale 0.6 #0")
    f = open(output,append)
    f.write(result)
    f.close()
    f = open("chimera.com", 'w')
    for c in commands:
        # print c
        print >> f, c
    f.close()
    s = l.get_system()
    append = 'a'
    base.Logger.log("Finished frame %i" % count, base.Logger.INFO)
    count += 1
f = open("sel_ref", "w")
f.write(sel_ref)
f.close()
f = open("colorbond.com", "w")
f.write(color)
f.close()
#TODO for single representations, split the system into sub systems,
#TODO several models in a pdb frame??
#TODO give options to build backbone only representation