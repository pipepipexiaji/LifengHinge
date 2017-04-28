#!/usr/bin/env python3

import sys, optparse, json, os, glob, re, random, string
from subprocess import Popen, PIPE
from math import *

###############################################
# subroutines 
###############################################

# fake nucleotide if sequence information is not given
fNt = 'A'

def next_nonskip_ih_ib(stap_or_scaf, ih, ib): 
    global vstrands

    while (vstrands[ih]['skip'][ib] == -1): 
        ih = vstrands[ih][stap_or_scaf][ib][2]
        ib = vstrands[ih][stap_or_scaf][ib][3]
        #last if $ih == -1;
    
    return (ih,ib)


def is_upstream(ih): 
    global vstrands

    for ib in range(len(vstrands[ih]['stap'])): 
        if vstrands[ih]['skip'][ib] == -1:
            continue
        nih = vstrands[ih]['stap'][ib][2]
        nib = vstrands[ih]['stap'][ib][3]
        if (nib == -1):
            continue

        if (ih == nih): 
            return 1 if nib > ib else -1
    
    sys.stderr.write("ERROR: Can't determine direction!\n")


def rotate_z(xyz, th):

    x, y, z = xyz
    th = th / 180.0 * pi
    return [x * cos(th) - y * sin(th), x * sin(th) + y * cos(th), z]

def flip_x(xyz):

    return [xyz[0], -xyz[1], -xyz[2]]

def wc(seq):
    if re.search("[atgcATGC]",seq):
        intab  = "ATGC"
        outtab = "TACG"
        trantab = str.maketrans(intab, outtab)
        wc = seq.upper().translate(trantab)

    else:
        raise Exception("wc(): wrong nucleotide '%s'." % seq)
    return wc

def na1to3(nt): 
    if nt == 'A': 
        return 'DA'
    
    elif nt == 'T': 
        return 'DT'
    
    elif nt == 'G': 
        return 'DG'
    
    elif nt == 'C': 
        return 'DC'
    


def get_direction_scaf(ih):
    Sum = vstrands[ih]['row'] + vstrands[ih]['col']
    
    return 1 if Sum % 2 == 0 else -1

def get_direction_stap(ih):
    Sum = vstrands[ih]['row'] + vstrands[ih]['col']
    
    return -1 if Sum % 2 == 0 else 1

def assign_seq_staple(ih,ib):

    global csv
    global vstrands
    global nseg

    # find the corresponding sequence in csv.
    for i in range(len(csv)): 
        if (csv[i]['5helix'] == ih and csv[i]['5base'] == ib):
        
            # now we found it.
            pseq = csv[i]['seq']
            nseg = i
            break
    try:
        pseq
    except NameError: 
        sys.stderr.write("Can't find a sequence for helix=%d,base=%d\n" % (ih, ib))

    ## trace the staple.
    while ib != -1:
        dir = get_direction_stap(ih)  # upstream=1, downstream=-1

        # next base
        pih, pib, nih, nib = vstrands[ih]['stap'][ib]
        if vstrands[ih]['skip'][ib] == -1:   # skip
            pass
        else:
            list_il = list(range(vstrands[ih]['loop'][ib] + 1)) # upstream
            if dir < 0: list_il = list(reversed(list_il)) # downstream

            vstrands[ih]['stap_nt'][ib] = {}
            vstrands[ih]['scaf_nt'][ib] = {}
            for il in list_il:
                nt = csv.pop(0)
                
                # Add sequence to vstrands.
                vstrands[ih]['stap_nt'][ib][il] = nt
                vstrands[ih]['scaf_nt'][ib][il] = wc(nt)

        (ih, ib) = (nih, nib)

    return

def assign_seq_scaffold(ih,ib):

    global vstrands
    global scaf_resi

    resi = 0

    ## trace the strand.
    while ib != -1:
        dir = get_direction_scaf(ih)  # upstream=1, downstream=-1

        # next base
        pih, pib, nih, nib = vstrands[ih]['scaf'][ib]
        if vstrands[ih]['skip'][ib] == -1:   # skip
            pass
        else:
            list_il = list(range(vstrands[ih]['loop'][ib] + 1)) # upstream
            if dir < 0: list_il = list(reversed(list_il)) # downstream

            vstrands[ih]['scaf_nt'][ib] = {}
            vstrands[ih]['stap_nt'][ib] = {}
            for il in list_il:
                resi += 1
                if (scaf_resi+resi-1) > len(scafseq[scafseqname]):
                    raise Exception("Error: %d is more than $scafseqname seq\n" % scaf_resi + resi - 1)
                
                else:
                    nt = scafseq[scafseqname][scaf_resi + resi - 1]
                    # Add sequence to vstrands.
                    vstrands[ih]['scaf_nt'][ib][il] = nt
                    vstrands[ih]['stap_nt'][ib][il] = wc(nt)

        (ih, ib) = (nih, nib)

    scaf_resi += resi

    return

def trace_scaffold(ih,ib):

    global scaf_resi
    global ch
    global vstrands
    global global_ri
    global fNt

    resi = 0
    ## trace the strand.
    while ib != -1:
        dir = get_direction_scaf(ih)      # upstream=1, downstream=-1
        pih, pib, nih, nib = vstrands[ih]['scaf'][ib]

        if vstrands[ih]['skip'][ib] == -1:   # skip
            pass
        else:
            # coordinates
            hrow = vstrands[ih]['row']     # helix row
            hcol = vstrands[ih]['col']     # helix column

            list_il = list(range(vstrands[ih]['loop'][ib] + 1))
            if dir < 0: list_il = list(reversed(list_il)) # downstream
            vstrands[ih]['scaf_ri'][ib] = {}
            vstrands[ih]['scaf_atoms'][ib] = {}
            for il in list_il:
                resi += 1

                if ib in vstrands[ih]['scaf_nt'] and il in vstrands[ih]['scaf_nt'][ib]:
                    nt = vstrands[ih]['scaf_nt'][ib][il]
                else:
                    nt = fNt
                # global residue count
                global_ri += 1
                vstrands[ih]['scaf_ri'][ib][il] = global_ri

                theta = vstrands[ih]['theta'][ib][il] 
                z = vstrands[ih]['z'][ib][il] 

                resn = na1to3(nt)
                if resi == 1: 
                    resn += '5';  # 5-ter
                
                elif nib == -1:
                    resn += '3';  # 3-ter
                
                atoms = print_pdb(hrow, hcol, ih, ib, il, nt, dir, "SCAF", resi, resn, theta, z)
                vstrands[ih]['scaf_atoms'][ib][il] = atoms
        # next base
        (ih, ib) = (nih, nib)

    print("TER")
    if ch == 90:
        ch = 65
    else:
        ch += 1
    scaf_resi += resi
    return

def trace_staple(ih,ib):

    global ch
    global vstrands
    global global_ri
    global nseg
    global fNt

    resi = 1
    ## trace the strand.
    while ib != -1:
        dir = get_direction_stap(ih)      # upstream=1, downstream=-1
        pih, pib, nih, nib = vstrands[ih]['stap'][ib]

        if vstrands[ih]['skip'][ib] == -1:   # skip

            vstrands[ih]['stap_nt'][ib] = {0:'X'}
            
        else:
            # coordinates
            hrow = vstrands[ih]['row']     # helix row
            hcol = vstrands[ih]['col']     # helix column

            list_il = list(range(vstrands[ih]['loop'][ib] + 1))
            if dir < 0: list_il = list(reversed(list_il)) # downstream
            vstrands[ih]['stap_ri'][ib] = {}
            vstrands[ih]['stap_atoms'][ib] = {}
            for il in list_il:

                if ib in vstrands[ih]['stap_nt'] and il in vstrands[ih]['stap_nt'][ib]:
                    nt = vstrands[ih]['stap_nt'][ib][il]
                else:
                    nt = fNt
                # global residue count
                global_ri += 1
                vstrands[ih]['stap_ri'][ib][il] = global_ri

                theta = vstrands[ih]['theta'][ib][il] 
                z = vstrands[ih]['z'][ib][il] 

                resn = na1to3(nt)
                if resi == 1: 
                    resn += '5';  # 5-ter
                
                elif nib == -1:
                    resn += '3';  # 3-ter

                seg = "P%03d" % nseg
                
                atoms = print_pdb(hrow, hcol, ih, ib, il, nt, dir, seg, resi, resn, theta, z)
                vstrands[ih]['stap_atoms'][ib][il] = atoms

                resi += 1

        # next base
        (ih, ib) = (nih, nib)

    if ch == 90:
        ch = 65
    else:
        ch += 1
    return

def print_pdb(ix, iy, ih, ib, il, seq, dir, seg, resid, resn, theta, z):

    global forcefield
    global lattice 
    global pdb
    global global_aid

    #### CHARMM 3TER rearrange atoms... dirty.
    #### Always doublecheck atom indices made by this script & the PSF file!!
    if (forcefield == "CHARMM" and re.search("3", resn)): 
        seq += '3'
    atoms = {}
    for i in pdb[seq].items():
        aname, xyz = i
        #print(xyz)
        
        # 5' terminal
        if resid == 1: 
            if re.search("P", aname): 
                if forcefield == "AMBER": 
                    aname = " H5T"
                
                else:
                    ## in case of CHARMM, put H5T before O5'. Dirty...
                    continue
            
            elif re.search("O[12]P", aname): 
                continue
            
        # flip PDB if downstream.
        if dir == -1: xyz = flip_x(xyz) 

        if lattice == "HONEYCOMB": 
            # rotate wrt z
            xyz = rotate_z(xyz, theta)
            ########################################################
            # 180 +- ??. Fine tuning needed.
            xyz = rotate_z(xyz, 170)

            # translate helix
            if ((iy + ix) % 2) == 0: 
                xyz[0] -= ix * ( 1 + 1 / sqrt(3)) * dd
            
            else: 
                xyz[0] -= (ix + (ix + 1) / sqrt(3)) * dd
            
            xyz[1] += iy * sqrt(3) / 2.0 * dd
            xyz[2] += z
        
        elif lattice == "SQUARE":
            # rotate wrt z
            xyz = rotate_z(xyz, theta)
            ########################################################
            # 180 +- ??. Fine tuning needed.
            xyz = rotate_z(xyz, 15)

            # translate helix
            xyz[0] -= ix * dd
            xyz[1] += iy * dd
            xyz[2] += z
        
        # CHARMM H5T before O5'
        if (forcefield == "CHARMM" and re.search("D.5", resn) and re.search("O5'", aname)):
            global_aid += 1
            atoms[aname] = global_aid
            if global_aid <= 99999:
                print("%-6s%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s" % ("ATOM", global_aid, " H5T", resn, chr(ch), resid, xyz[0]+random.random()/2, xyz[1]+random.random()/2, xyz[2]+random.random()/2,1,0, seg))
            else:
                print("%-6s%5s %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s" % ("ATOM", "*****", " H5T", resn, chr(ch), resid, xyz[0]+random.random()/2, xyz[1]+random.random()/2, xyz[2]+random.random()/2,1,0, seg))

        global_aid += 1
        atoms[aname] = global_aid
        if global_aid <= 99999:
            print("%-6s%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s" % ("ATOM", global_aid, aname, resn, chr(ch), resid, xyz[0], xyz[1], xyz[2], 1,0, seg))
        else:
            print("%-6s%5s %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s" % ("ATOM", "*****", aname, resn, chr(ch), resid, xyz[0], xyz[1], xyz[2], 1,0, seg))

        
        # N1|N3 atom position
        if ((re.search("(DA|DG)", resn) and re.search("N1", aname))  or
            (re.search("(DT|DC)", resn) and re.search("N3", aname))):
        
            vstrands[ih]["x%s" % dir][ib] = {}
            vstrands[ih]["y%s" % dir][ib] = {}
            vstrands[ih]["z%s" % dir][ib] = {}
            vstrands[ih]["x%s" % dir][ib][il] = xyz[0]
            vstrands[ih]["y%s" % dir][ib][il] = xyz[1]
            vstrands[ih]["z%s" % dir][ib][il] = xyz[2]
        
        # H3T
        if (re.search("D.3", resn) and re.search("O3'", aname)): 
            global_aid += 1
            atoms[aname] = global_aid
            if global_aid <= 99999:
                print("%-6s%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s" % ("ATOM", global_aid, " H3T", resn, chr(ch), resid, xyz[0]+random.random()/2, xyz[1]+random.random()/2, xyz[2]+random.random()/2, 1, 0, seg))
            else:
                print("%-6s%5s %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s" % ("ATOM", "*****", " H3T", resn, chr(ch), resid, xyz[0]+random.random()/2, xyz[1]+random.random()/2, xyz[2]+random.random()/2, 1, 0, seg))

    return atoms

def init_pdb_charmm():
    pdb = {
    "A": {
            " P  ": [   0.288, -9.220,  -1.848 ],
	    " C4'": [   3.212, -6.864,  -1.355 ],
	    " H4'": [   4.069, -6.713,  -1.848 ],
	    " O4'": [   2.387, -5.664,  -1.352 ],
	    " C1'": [   2.281, -5.198,  -0.016 ],
	    " H1'": [   3.030, -4.578,   0.218 ],
	    " C2'": [   2.304, -6.454,   0.850 ],
	    " H2'": [   1.451, -6.973,   0.796 ],
	    "H2''": [   2.514, -6.256,   1.807 ],
	    " O1P": [   0.421,-10.485,  -2.605 ],
	    " O2P": [  -0.692, -9.226,  -0.740 ],
	    " O5'": [   1.721, -8.770,  -1.295 ],
	    " C5'": [   2.585, -7.995,  -2.146 ],
	    " H5'": [   2.038, -7.632,  -2.900 ],
	    "H5''": [   3.293, -8.601,  -2.509 ],
	    " N9 ": [   1.021, -4.412,   0.101 ],
	    " C5 ": [  -0.409, -2.735,   0.190 ],
	    " N7 ": [  -1.161, -3.897,   0.354 ],
	    " C8 ": [  -0.268, -4.856,   0.294 ],
	    " H8 ": [  -0.503, -5.824,   0.381 ],
	    " N1 ": [   0.269, -0.516,  -0.025 ],
	    " C2 ": [   1.584, -0.914,  -0.172 ],
	    " H2 ": [   2.282, -0.210,  -0.302 ],
	    " N3 ": [   1.969, -2.189,  -0.149 ],
	    " C4 ": [   0.923, -3.037,   0.035 ],
	    " C6 ": [  -0.820, -1.378,   0.165 ],
	    " N6 ": [  -1.951, -0.913,   0.284 ],
	    " H61": [  -2.731, -1.524,   0.420 ],
	    " H62": [  -2.091,  0.076,   0.245 ],
	    " C3'": [   3.458, -7.146,   0.127 ],
	    " H3'": [   3.568, -8.115,   0.350 ],
	    " O3'": [   4.677, -6.520,   0.518 ]},
    "T": {  
	    " P  ": [   0.249, -9.221,  -1.866 ],
	    " C4'": [   3.176, -6.866,  -1.373 ],
	    " H4'": [   4.034, -6.717,  -1.864 ],
	    " O4'": [   2.348, -5.666,  -1.370 ],
	    " C1'": [   2.243, -5.199,  -0.034 ],
	    " H1'": [   2.994, -4.581,   0.200 ],
	    " C2'": [   2.265, -6.453,   0.832 ],
	    " H2'": [   1.411, -6.971,   0.778 ],
	    "H2''": [   2.474, -6.254,   1.789 ],
	    " O1P": [   0.383,-10.486,  -2.623 ],
	    " O2P": [  -0.730, -9.227,  -0.758 ],
	    " O5'": [   1.683, -8.771,  -1.313 ],
	    " C5'": [   2.547, -7.996,  -2.164 ],
	    " H5'": [   2.000, -7.632,  -2.918 ],
	    "H5''": [   3.254, -8.602,  -2.528 ],
	    " N1 ": [   0.982, -4.412,   0.083 ],
	    " C6 ": [  -0.215, -5.030,   0.274 ],
	    " H6 ": [  -0.250, -6.027,   0.336 ],
	    " C2 ": [   1.072, -3.026,  -0.008 ],
	    " O2 ": [   2.183, -2.509,  -0.181 ],
	    " N3 ": [  -0.065, -2.290,   0.097 ],
	    " H3 ": [  -0.018, -1.293,   0.034 ],
	    " C4 ": [  -1.249, -2.888,   0.284 ],
	    " O4 ": [  -2.329, -2.129,   0.380 ],
	    " C5 ": [  -1.361, -4.314,   0.381 ],
	    " C5M": [  -2.245, -4.760,   0.521 ],
	    " H51": [  -2.964, -4.067,   0.574 ],
	    " H52": [  -2.432, -5.376,  -0.244 ],
	    " H53": [  -2.223, -5.283,   1.373 ],
	    " C3'": [   3.421, -7.147,   0.109 ],
	    " H3'": [   3.532, -8.115,   0.332 ],
	    " O3'": [   4.638, -6.523,   0.500 ]},
    "G": {
	    " P  ": [   0.288, -9.220,  -1.848 ],
	    " C4'": [   3.212, -6.864,  -1.355 ],
	    " H4'": [   4.069, -6.713,  -1.848 ],
	    " O4'": [   2.387, -5.664,  -1.352 ],
	    " C1'": [   2.281, -5.198,  -0.016 ],
	    " H1'": [   3.030, -4.578,   0.218 ],
	    " C2'": [   2.304, -6.454,   0.850 ],
	    " H2'": [   1.451, -6.973,   0.796 ],
	    "H2''": [   2.514, -6.256,   1.807 ],
	    " O1P": [   0.421,-10.485,  -2.605 ],
	    " O2P": [  -0.692, -9.226,  -0.740 ],
	    " O5'": [   1.721, -8.770,  -1.295 ],
	    " C5'": [   2.585, -7.995,  -2.146 ],
	    " H5'": [   2.038, -7.632,  -2.900 ],
	    "H5''": [   3.293, -8.601,  -2.509 ],
	    " N9 ": [   1.021, -4.412,   0.101 ],
	    " C4 ": [   0.923, -3.037,   0.035 ],
	    " N2 ": [   2.478,  0.061,  -0.344 ],
	    " H21": [   2.180,  1.015,  -0.361 ],
	    " H22": [   3.446, -0.163,  -0.456 ],
	    " N3 ": [   1.969, -2.189,  -0.149 ],
	    " C2 ": [   1.584, -0.914,  -0.172 ],
	    " N1 ": [   0.269, -0.516,  -0.025 ],
	    " H1 ": [   0.079,  0.465,  -0.057 ],
	    " C6 ": [  -0.820, -1.378,   0.165 ],
	    " O6 ": [  -1.951, -0.913,   0.284 ],
	    " C5 ": [  -0.409, -2.735,   0.190 ],
	    " N7 ": [  -1.161, -3.897,   0.354 ],
	    " C8 ": [  -0.268, -4.856,   0.294 ],
	    " H8 ": [  -0.503, -5.824,   0.381 ],
	    " C3'": [   3.458, -7.146,   0.127 ],
	    " H3'": [   3.568, -8.115,   0.350 ],
	    " O3'": [   4.677, -6.520,   0.518 ]},
    "C": {
	    " P  ": [   0.249, -9.221,  -1.866 ],
	    " C4'": [   3.176, -6.866,  -1.373 ],
	    " H4'": [   4.034, -6.717,  -1.864 ],
	    " O4'": [   2.348, -5.666,  -1.370 ],
	    " C1'": [   2.243, -5.199,  -0.034 ],
	    " H1'": [   2.994, -4.581,   0.200 ],
	    " C2'": [   2.265, -6.453,   0.832 ],
	    " H2'": [   1.411, -6.971,   0.778 ],
	    "H2''": [   2.474, -6.254,   1.789 ],
	    " O1P": [   0.383,-10.486,  -2.623 ],
	    " O2P": [  -0.730, -9.227,  -0.758 ],
	    " O5'": [   1.683, -8.771,  -1.313 ],
	    " C5'": [   2.547, -7.996,  -2.164 ],
	    " H5'": [   2.000, -7.632,  -2.918 ],
	    "H5''": [   3.254, -8.602,  -2.528 ],
	    " N1 ": [   0.982, -4.412,   0.083 ],
	    " C6 ": [  -0.215, -5.030,   0.274 ],
	    " H6 ": [  -0.250, -6.027,   0.336 ],
	    " C5 ": [  -1.361, -4.314,   0.381 ],
	    " H5 ": [  -2.245, -4.760,   0.521 ],
	    " C2 ": [   1.072, -3.026,  -0.008 ],
	    " O2 ": [   2.183, -2.509,  -0.181 ],
	    " N3 ": [  -0.065, -2.290,   0.097 ],
	    " C4 ": [  -1.249, -2.888,   0.284 ],
	    " N4 ": [  -2.329, -2.129,   0.380 ],
	    " H41": [  -3.225, -2.550,   0.520 ],
	    " H42": [  -2.249, -1.135,   0.312 ],
	    " C3'": [   3.421, -7.147,   0.109 ],
	    " H3'": [   3.532, -8.115,   0.332 ],
	    " O3'": [   4.638, -6.523,   0.500 ]},
    "A3": {
            " P  ": [   0.288, -9.220,  -1.848 ],
	    " C4'": [   3.212, -6.864,  -1.355 ],
	    " H4'": [   4.069, -6.713,  -1.848 ],
	    " O4'": [   2.387, -5.664,  -1.352 ],
	    " C1'": [   2.281, -5.198,  -0.016 ],
	    " H1'": [   3.030, -4.578,   0.218 ],
	    " C2'": [   2.304, -6.454,   0.850 ],
	    " H2'": [   1.451, -6.973,   0.796 ],
	    "H2''": [   2.514, -6.256,   1.807 ],
	    " C3'": [   3.458, -7.146,   0.127 ],
	    " H3'": [   3.568, -8.115,   0.350 ],
	    " O3'": [   4.677, -6.520,   0.518 ],
	    " O1P": [   0.421,-10.485,  -2.605 ],
	    " O2P": [  -0.692, -9.226,  -0.740 ],
	    " O5'": [   1.721, -8.770,  -1.295 ],
	    " C5'": [   2.585, -7.995,  -2.146 ],
	    " H5'": [   2.038, -7.632,  -2.900 ],
	    "H5''": [   3.293, -8.601,  -2.509 ],
	    " N9 ": [   1.021, -4.412,   0.101 ],
	    " C5 ": [  -0.409, -2.735,   0.190 ],
	    " N7 ": [  -1.161, -3.897,   0.354 ],
	    " C8 ": [  -0.268, -4.856,   0.294 ],
	    " H8 ": [  -0.503, -5.824,   0.381 ],
	    " N1 ": [   0.269, -0.516,  -0.025 ],
	    " C2 ": [   1.584, -0.914,  -0.172 ],
	    " H2 ": [   2.282, -0.210,  -0.302 ],
	    " N3 ": [   1.969, -2.189,  -0.149 ],
	    " C4 ": [   0.923, -3.037,   0.035 ],
	    " C6 ": [  -0.820, -1.378,   0.165 ],
	    " N6 ": [  -1.951, -0.913,   0.284 ],
	    " H61": [  -2.731, -1.524,   0.420 ],
	    " H62": [  -2.091,  0.076,   0.245 ]},
    "T3": { 
	    " P  ": [   0.249, -9.221,  -1.866 ],
	    " C4'": [   3.176, -6.866,  -1.373 ],
	    " H4'": [   4.034, -6.717,  -1.864 ],
	    " O4'": [   2.348, -5.666,  -1.370 ],
	    " C1'": [   2.243, -5.199,  -0.034 ],
	    " H1'": [   2.994, -4.581,   0.200 ],
	    " C2'": [   2.265, -6.453,   0.832 ],
	    " H2'": [   1.411, -6.971,   0.778 ],
	    "H2''": [   2.474, -6.254,   1.789 ],
	    " C3'": [   3.421, -7.147,   0.109 ],
	    " H3'": [   3.532, -8.115,   0.332 ],
	    " O3'": [   4.638, -6.523,   0.500 ],
	    " O1P": [   0.383,-10.486,  -2.623 ],
	    " O2P": [  -0.730, -9.227,  -0.758 ],
	    " O5'": [   1.683, -8.771,  -1.313 ],
	    " C5'": [   2.547, -7.996,  -2.164 ],
	    " H5'": [   2.000, -7.632,  -2.918 ],
	    "H5''": [   3.254, -8.602,  -2.528 ],
	    " N1 ": [   0.982, -4.412,   0.083 ],
	    " C6 ": [  -0.215, -5.030,   0.274 ],
	    " H6 ": [  -0.250, -6.027,   0.336 ],
	    " C2 ": [   1.072, -3.026,  -0.008 ],
	    " O2 ": [   2.183, -2.509,  -0.181 ],
	    " N3 ": [  -0.065, -2.290,   0.097 ],
	    " H3 ": [  -0.018, -1.293,   0.034 ],
	    " C4 ": [  -1.249, -2.888,   0.284 ],
	    " O4 ": [  -2.329, -2.129,   0.380 ],
	    " C5 ": [  -1.361, -4.314,   0.381 ],
	    " C5M": [  -2.245, -4.760,   0.521 ],
	    " H51": [  -2.964, -4.067,   0.574 ],
	    " H52": [  -2.432, -5.376,  -0.244 ],
	    " H53": [  -2.223, -5.283,   1.373 ]},
    "G3": {
	    " P  ": [   0.288, -9.220,  -1.848 ],
	    " C4'": [   3.212, -6.864,  -1.355 ],
	    " H4'": [   4.069, -6.713,  -1.848 ],
	    " O4'": [   2.387, -5.664,  -1.352 ],
	    " C1'": [   2.281, -5.198,  -0.016 ],
	    " H1'": [   3.030, -4.578,   0.218 ],
	    " C2'": [   2.304, -6.454,   0.850 ],
	    " H2'": [   1.451, -6.973,   0.796 ],
	    "H2''": [   2.514, -6.256,   1.807 ],
	    " C3'": [   3.458, -7.146,   0.127 ],
	    " H3'": [   3.568, -8.115,   0.350 ],
	    " O3'": [   4.677, -6.520,   0.518 ],
	    " O1P": [   0.421,-10.485,  -2.605 ],
	    " O2P": [  -0.692, -9.226,  -0.740 ],
	    " O5'": [   1.721, -8.770,  -1.295 ],
	    " C5'": [   2.585, -7.995,  -2.146 ],
	    " H5'": [   2.038, -7.632,  -2.900 ],
	    "H5''": [   3.293, -8.601,  -2.509 ],
	    " N9 ": [   1.021, -4.412,   0.101 ],
	    " C4 ": [   0.923, -3.037,   0.035 ],
	    " N2 ": [   2.478,  0.061,  -0.344 ],
	    " H21": [   2.180,  1.015,  -0.361 ],
	    " H22": [   3.446, -0.163,  -0.456 ],
	    " N3 ": [   1.969, -2.189,  -0.149 ],
	    " C2 ": [   1.584, -0.914,  -0.172 ],
	    " N1 ": [   0.269, -0.516,  -0.025 ],
	    " H1 ": [   0.079,  0.465,  -0.057 ],
	    " C6 ": [  -0.820, -1.378,   0.165 ],
	    " O6 ": [  -1.951, -0.913,   0.284 ],
	    " C5 ": [  -0.409, -2.735,   0.190 ],
	    " N7 ": [  -1.161, -3.897,   0.354 ],
	    " C8 ": [  -0.268, -4.856,   0.294 ],
	    " H8 ": [  -0.503, -5.824,   0.381 ]},
    "C3": {
	    " P  ": [   0.249, -9.221,  -1.866 ],
	    " C4'": [   3.176, -6.866,  -1.373 ],
	    " H4'": [   4.034, -6.717,  -1.864 ],
	    " O4'": [   2.348, -5.666,  -1.370 ],
	    " C1'": [   2.243, -5.199,  -0.034 ],
	    " H1'": [   2.994, -4.581,   0.200 ],
	    " C2'": [   2.265, -6.453,   0.832 ],
	    " H2'": [   1.411, -6.971,   0.778 ],
	    "H2''": [   2.474, -6.254,   1.789 ],
	    " C3'": [   3.421, -7.147,   0.109 ],
	    " H3'": [   3.532, -8.115,   0.332 ],
	    " O3'": [   4.638, -6.523,   0.500 ],
	    " O1P": [   0.383,-10.486,  -2.623 ],
	    " O2P": [  -0.730, -9.227,  -0.758 ],
	    " O5'": [   1.683, -8.771,  -1.313 ],
	    " C5'": [   2.547, -7.996,  -2.164 ],
	    " H5'": [   2.000, -7.632,  -2.918 ],
	    "H5''": [   3.254, -8.602,  -2.528 ],
	    " N1 ": [   0.982, -4.412,   0.083 ],
	    " C6 ": [  -0.215, -5.030,   0.274 ],
	    " H6 ": [  -0.250, -6.027,   0.336 ],
	    " C5 ": [  -1.361, -4.314,   0.381 ],
	    " H5 ": [  -2.245, -4.760,   0.521 ],
	    " C2 ": [   1.072, -3.026,  -0.008 ],
	    " O2 ": [   2.183, -2.509,  -0.181 ],
	    " N3 ": [  -0.065, -2.290,   0.097 ],
	    " C4 ": [  -1.249, -2.888,   0.284 ],
	    " N4 ": [  -2.329, -2.129,   0.380 ],
	    " H41": [  -3.225, -2.550,   0.520 ],
	    " H42": [  -2.249, -1.135,   0.312 ]}
    }

    return pdb
def init_pdb_amber():

    pdb = {
    "A": {
            " P  ": [   0.288, -9.220,  -1.848 ],
	    " O1P": [   0.421,-10.485,  -2.605 ],
	    " O2P": [  -0.692, -9.226,  -0.740 ],
	    " O5'": [   1.721, -8.770,  -1.295 ],
	    " C5'": [   2.585, -7.995,  -2.146 ],
	    "H5'1": [   2.038, -7.632,  -2.900 ],
	    "H5'2": [   3.293, -8.601,  -2.509 ],
	    " C4'": [   3.212, -6.864,  -1.355 ],
	    " H4'": [   4.069, -6.713,  -1.848 ],
	    " O4'": [   2.387, -5.664,  -1.352 ],
	    " C1'": [   2.281, -5.198,  -0.016 ],
	    " H1'": [   3.030, -4.578,   0.218 ],
	    " N9 ": [   1.021, -4.412,   0.101 ],
	    " C8 ": [  -0.268, -4.856,   0.294 ],
	    " H8 ": [  -0.503, -5.824,   0.381 ],
	    " N7 ": [  -1.161, -3.897,   0.354 ],
	    " C5 ": [  -0.409, -2.735,   0.190 ],
	    " C6 ": [  -0.820, -1.378,   0.165 ],
	    " N6 ": [  -1.951, -0.913,   0.284 ],
	    " H61": [  -2.731, -1.524,   0.420 ],
	    " H62": [  -2.091,  0.076,   0.245 ],
	    " N1 ": [   0.269, -0.516,  -0.025 ],
	    " C2 ": [   1.584, -0.914,  -0.172 ],
	    " H2 ": [   2.282, -0.210,  -0.302 ],
	    " N3 ": [   1.969, -2.189,  -0.149 ],
	    " C4 ": [   0.923, -3.037,   0.035 ],
	    " C3'": [   3.458, -7.146,   0.127 ],
	    " H3'": [   3.568, -8.115,   0.350 ],
	    " C2'": [   2.304, -6.454,   0.850 ],
	    "H2'1": [   1.451, -6.973,   0.796 ],
	    "H2'2": [   2.514, -6.256,   1.807 ],
	    " O3'": [   4.677, -6.520,   0.518 ]},
    "T": {
            " P  ": [   0.249, -9.221,  -1.866 ],
	    " O1P": [   0.383,-10.486,  -2.623 ],
	    " O2P": [  -0.730, -9.227,  -0.758 ],
	    " O5'": [   1.683, -8.771,  -1.313 ],
	    " C5'": [   2.547, -7.996,  -2.164 ],
	    "H5'1": [   2.000, -7.632,  -2.918 ],
	    "H5'2": [   3.254, -8.602,  -2.528 ],
	    " C4'": [   3.176, -6.866,  -1.373 ],
	    " H4'": [   4.034, -6.717,  -1.864 ],
	    " O4'": [   2.348, -5.666,  -1.370 ],
	    " C1'": [   2.243, -5.199,  -0.034 ],
	    " H1'": [   2.994, -4.581,   0.200 ],
	    " N1 ": [   0.982, -4.412,   0.083 ],
	    " C6 ": [  -0.215, -5.030,   0.274 ],
	    " H6 ": [  -0.250, -6.027,   0.336 ],
	    " C5 ": [  -1.361, -4.314,   0.381 ],
	    " C7 ": [  -2.245, -4.760,   0.521 ],
	    " H71": [  -2.964, -4.067,   0.574 ],
	    " H72": [  -2.432, -5.376,  -0.244 ],
	    " H73": [  -2.223, -5.283,   1.373 ],
	    " C4 ": [  -1.249, -2.888,   0.284 ],
	    " O4 ": [  -2.329, -2.129,   0.380 ],
	    " N3 ": [  -0.065, -2.290,   0.097 ],
	    " H3 ": [  -0.018, -1.293,   0.034 ],
	    " C2 ": [   1.072, -3.026,  -0.008 ],
	    " O2 ": [   2.183, -2.509,  -0.181 ],
	    " C3'": [   3.421, -7.147,   0.109 ],
	    " H3'": [   3.532, -8.115,   0.332 ],
	    " C2'": [   2.265, -6.453,   0.832 ],
	    "H2'1": [   1.411, -6.971,   0.778 ],
	    "H2'2": [   2.474, -6.254,   1.789 ],
	    " O3'": [   4.638, -6.523,   0.500 ]},
    "G": {
            " P  ": [   0.288, -9.220,  -1.848 ],
	    " O1P": [   0.421,-10.485,  -2.605 ],
	    " O2P": [  -0.692, -9.226,  -0.740 ],
	    " O5'": [   1.721, -8.770,  -1.295 ],
	    " C5'": [   2.585, -7.995,  -2.146 ],
	    "H5'1": [   2.038, -7.632,  -2.900 ],
	    "H5'2": [   3.293, -8.601,  -2.509 ],
	    " C4'": [   3.212, -6.864,  -1.355 ],
	    " H4'": [   4.069, -6.713,  -1.848 ],
	    " O4'": [   2.387, -5.664,  -1.352 ],
	    " C1'": [   2.281, -5.198,  -0.016 ],
	    " H1'": [   3.030, -4.578,   0.218 ],
	    " N9 ": [   1.021, -4.412,   0.101 ],
	    " C8 ": [  -0.268, -4.856,   0.294 ],
	    " H8 ": [  -0.503, -5.824,   0.381 ],
	    " N7 ": [  -1.161, -3.897,   0.354 ],
	    " C5 ": [  -0.409, -2.735,   0.190 ],
	    " C6 ": [  -0.820, -1.378,   0.165 ],
	    " O6 ": [  -1.951, -0.913,   0.284 ],
	    " N1 ": [   0.269, -0.516,  -0.025 ],
	    " H1 ": [   0.079,  0.465,  -0.057 ],
	    " C2 ": [   1.584, -0.914,  -0.172 ],
	    " N2 ": [   2.478,  0.061,  -0.344 ],
	    " H21": [   2.180,  1.015,  -0.361 ],
	    " H22": [   3.446, -0.163,  -0.456 ],
	    " N3 ": [   1.969, -2.189,  -0.149 ],
	    " C4 ": [   0.923, -3.037,   0.035 ],
	    " C3'": [   3.458, -7.146,   0.127 ],
	    " H3'": [   3.568, -8.115,   0.350 ],
	    " C2'": [   2.304, -6.454,   0.850 ],
	    "H2'1": [   1.451, -6.973,   0.796 ],
	    "H2'2": [   2.514, -6.256,   1.807 ],
	    " O3'": [   4.677, -6.520,   0.518 ]},
    "C": {
            " P  ": [   0.249, -9.221,  -1.866 ],
	    " O1P": [   0.383,-10.486,  -2.623 ],
	    " O2P": [  -0.730, -9.227,  -0.758 ],
	    " O5'": [   1.683, -8.771,  -1.313 ],
	    " C5'": [   2.547, -7.996,  -2.164 ],
	    "H5'1": [   2.000, -7.632,  -2.918 ],
	    "H5'2": [   3.254, -8.602,  -2.528 ],
	    " C4'": [   3.176, -6.866,  -1.373 ],
	    " H4'": [   4.034, -6.717,  -1.864 ],
	    " O4'": [   2.348, -5.666,  -1.370 ],
	    " C1'": [   2.243, -5.199,  -0.034 ],
	    " H1'": [   2.994, -4.581,   0.200 ],
	    " N1 ": [   0.982, -4.412,   0.083 ],
	    " C6 ": [  -0.215, -5.030,   0.274 ],
	    " H6 ": [  -0.250, -6.027,   0.336 ],
	    " C5 ": [  -1.361, -4.314,   0.381 ],
	    " H5 ": [  -2.245, -4.760,   0.521 ],
	    " C4 ": [  -1.249, -2.888,   0.284 ],
	    " N4 ": [  -2.329, -2.129,   0.380 ],
	    " H41": [  -3.225, -2.550,   0.520 ],
	    " H42": [  -2.249, -1.135,   0.312 ],
	    " N3 ": [  -0.065, -2.290,   0.097 ],
	    " C2 ": [   1.072, -3.026,  -0.008 ],
	    " O2 ": [   2.183, -2.509,  -0.181 ],
	    " C3'": [   3.421, -7.147,   0.109 ],
	    " H3'": [   3.532, -8.115,   0.332 ],
	    " C2'": [   2.265, -6.453,   0.832 ],
	    "H2'1": [   1.411, -6.971,   0.778 ],
	    "H2'2": [   2.474, -6.254,   1.789 ],
	    " O3'": [   4.638, -6.523,   0.500 ]}
    }


    return pdb


## parse command-line arguments
usage="python3 cadnano2pdb.py --ff=(CHARMM|AMBER) --lattice=(honeycomb|square) [--scaf=SCAFFOLD_SEQUENCE] [--stap=STAPLE_CSV_FILE] [--help] JSON_FILE"
description='''
The cadnano2pdb utility converts a JSON file from caDNAno program to a PDB file and write the PDB file to STDOUT. 
JSON_FILE is mandatory, and the other options are optional. The options not specified in the commandline will be determined interactively.
Scaffold and staple sequence information can be provided using --scaf and --staple options, respectively. 
When only one of those two are given, the sequence of the other will be determined complementarily. 
'''

parser = optparse.OptionParser(usage=usage,
description=description)
parser.add_option('--ff', default='CHARMM', action="store", type="string", dest="ff", help='(CHARMM|AMBER): determine either CHARMM- and AMBER-compatible PDB outputs. For example, nucleotide and atoms names and the order of atoms in a nucleotide will be consistent with the chosen force field.')
parser.add_option('--lattice', default='square', action="store", type="string", dest="lattice", help='(honeycomb|square): specify the lattice type in which JSON file was created.')
parser.add_option('--noloop', default=None, action="store_true", dest="noloop")
parser.add_option('--noskip', default=None, action="store_true", dest="noskip")
parser.add_option('--stap', default=None, action="store", type="string", dest="stap", help='STAPLE_CSV_FILE: CSV file from caDNAno.')
parser.add_option('--scaf', default=None, action="store", type="string", dest="scaf", help='SCAFFOLD_SEQUENCE: the sequence of scaffold. The sequence data is taken from the caDNAno program. When not given, staple sequences must be given using --stap option.')

## get the options and arguments
(options, args) = parser.parse_args()

forcefield=options.ff.upper()
lattice=options.lattice.upper()
noloop=options.noloop
noskip=options.noskip
csvfile=options.stap
scafseqname=options.scaf


## get the scaffold sequence from the sequence file
pwd=os.getcwd()
seq_files = glob.glob("%s/*.txt" % (pwd))

scafseq = {}
for f in seq_files:
    m = re.search("*.txt",f)
    name = m.group(1)
    with open(f, 'r') as fseq:
        tmp = []
        for line in fseq:

            if re.search("^;(.*)", line): 
                continue
            l1 = line.replace(' ','')
            l2 = l1.strip()
            l3 = l2.replace('\n','')
            l4 = l3.upper()
            tmp.extend(list(l4))

    scafseq[name] = tmp

try: 
    scafseq
except ValueError:
    sys.stderr.write("Error: no sequence files are found. Place sequence files in the same directory as the perl script.\n")


if len(args) < 1:
    parser.print_help()
    sys.stderr.write("\n\t\tSCAFFOLD_SEQUENCE=\n")
    for s in scafseq.keys():
        sys.stderr.write("\t\t\t %-15s (%5d nts)\n" % (s, len(scafseq[s])))
    raise Exception("Insufficient arguments")


## global variables
scaf_resi=0    # global variable for counting multiple scaffold strands.
dd = 23.0  # interaxial distance in A


sys.stderr.write("json=%s\n" % args[0])

### select lattice type.
while not (lattice == "HONEYCOMB" or lattice == "SQUARE"):  
    sys.stderr.write("Warning: --lattice=%s is undefined. Choose again!\n" % lattice)
    sys.stderr.write("Select lattice type (honeycomb or square): ")
    lattice = input()
    lattice = lattice.upper()


### select forcefield type.
while not (forcefield == "CHARMM" or forcefield == "AMBER"):  
    sys.stderr.write("Warning: --forcefield=%s is undefined. Choose again!\n" % forcefield)
    sys.stderr.write("Select forcefield type (charmm or amber): ")
    forcefield = input()
    forcefield = forcefield.upper()


####################################
# Read PDB of A,T,G,C
# - AT and GC pairs by rotating one NT by 180 wrt x axis.
# - AT and GC pairs are parallel to x axis.

if (forcefield == "CHARMM"):
    sys.stderr.write("CHARMM36 force field\n")
    pdb = init_pdb_charmm()
    
else:
    sys.stderr.write("AMBER99 force field\n")
    pdb = init_pdb_amber()


# determine sequence
while scafseqname and not scafseqname in scafseq: 
    sys.stderr.write("Warning: --scaf=%s is not defined.\n" % scafseqname)
    sys.stderr.write("\t\tSCAFFOLD_SEQUENCE=\n")
    for s in scafseq.keys():
        sys.stderr.write("\t\t\t %-15s (%5d nts)\n" % (s, len(scafseq[s])+1))
    sys.stderr.write("Select SCAFFOLD_SEQUENCE: ")
    scafseqname = input()
    
if not scafseqname and not csvfile: 
    raise Exception("at least one of --scaf and --stap options must be given")


####################################
# Read sequence (CSV)
csv = []
if csvfile:
    with open(csvfile, 'r') as FH:
        sys.stderr.write("Reading CSV file %s\n" % csvfile)
        for line in FH:
            m = re.search("([0-9]+)\[([0-9]+)\],([0-9]+)\[([0-9]+)\],([ATGC]+),[0-9]+", line)
            if m:
                tmp = {
                    '5helix' : int(m.group(1)),
                    '5base'  : int(m.group(2)),
                    '3helix' : int(m.group(3)),
                    '3base'  : int(m.group(4)),
                    'seq'    : list(m.group(5))
                }
                csv.append(tmp)

    sys.stderr.write("#staples in %s = %d\n" % (csvfile, len(csv)))

#print(csv)

####################################
# Read json
# Indices for helix and base starts from 0.
# @vstrands array has a hash of these:
#       num : helix index from 0. NOTE that @vstrands order and $num don't match! 
#             Use ih2id() function.
#             Order of helix id doesn't matter because helix positions determined by (row,col)
#       stap    : connectivity array of bases [5' helix, 5' base, 3' helix, 3' base]
#       scaf    : connectivity array of bases
#       row     
#       col     
#       stap_colors     
#       skip        ; deletion
#       loop        ; insertion
#       scafLoop    ; deprecated?
#       stapLoop    ; deprecated?
### defined by me #############################
#       scaf_nt    : sequence
#       stap_nt    : sequence
#       stap_ri    : staple residue id (one-based)
#       scaf_ri    : scaffold residue id (one-based)
#       theta       : rotation considering insertion/deletion.
#       z           : z coordinate considering insertion/deletion.
#       stap_atoms : hash of atom name & atom index
#       scaf_atoms : hash of atom name & atom index
#       hb_major        : two atom indices of H-bond in major groove. 
#       hb_minor        : two atom indices of H-bond in minor groove. DON'T USE IT!
#       hb_center       : two atom indices of H-bond in center. DON'T USE IT!
#       c1p_scaf                : atom index of C1' atom DON'T USE IT!
#       c1p_stap                : atom index of C1' atom DON'T USE IT!
#       ( N1 | C2 | N3 | C4 | C5 | C6 | N7 | C8 | N9 ) : base atom indices
####################################

with open(args[0], 'r') as inFile:
    data = json.load(inFile)

vstrands = data['vstrands']
#print(len(vstrands[0]['scaf']))

ih2id={}
for i in range(len(vstrands)):
    ih = vstrands[i]['num']
    ih2id[ih] = i

####################### 
# renumber $ih so that we can forget about it.
for ih in range(len(vstrands)): 
    # scaffold
    for ib in range(len(vstrands[ih]['scaf'])): 
        pih, pib, nih, nib = vstrands[ih]['scaf'][ib]
        if (pih > -1): vstrands[ih]['scaf'][ib][0] = ih2id[pih] 
        if (nih > -1): vstrands[ih]['scaf'][ib][2] = ih2id[nih] 
    
    # staples
    for ib in range(len(vstrands[ih]['stap'])): 
        pih, pib, nih, nib = vstrands[ih]['stap'][ib]
        if (pih > -1): vstrands[ih]['stap'][ib][0] = ih2id[pih] 
        if (nih > -1): vstrands[ih]['stap'][ib][2] = ih2id[nih] 
    
for i in range(len(csv)): 
    csv[i]['5helix'] = ih2id[int(csv[i]['5helix'])]

global_ri = 0


## Ignore skip or loop if needed. Shit. csv should be also modified accordingly.
if noloop: 
    for ih in range(len(vstrands)): 
        for ib in range(len(vstrands[ih]['loop'])): 
            vstrands[ih]['loop'][ib] = 0

if noskip: 
    for ih in range(len(vstrands)): 
        for ib in range(len(vstrands[ih]['skip'])): 
            vstrands[ih]['skip'][ib] = 0
        

####################################################
# determine theta and z considering insertion/deletion
# For honeycomb, 2 turns (720 degrees) per 21 bps.
# For square, 3 turns (1080 degrees) per 32 bps.
nblock = 7 if (lattice == "HONEYCOMB") else 8
theta_per_bp = 360.0*2.0/21 if (lattice == "HONEYCOMB") else 360.0*3.0/32
z_per_bp = 3.4 

for ih in range(len(vstrands)):

    vstrands[ih]['theta'] = {}
    vstrands[ih]['z'] = {}

for ih in range(len(vstrands)):

    ib = 0
    while ib < int(len(vstrands[ih]['stap']) / nblock) * nblock:
        dn = 0
        #sys.stderr.write("%d, %d, %d\n" % (ih, ib, len(vstrands[ih]['stap'])))
        for i in range(ib, ib + nblock): 
            dn += vstrands[ih]['loop'][i] + vstrands[ih]['skip'][i]

        resi=0
        theta0 = (int(float(ib)/nblock)) * nblock * theta_per_bp
        z0 = (int(float(ib)/nblock)) * nblock * z_per_bp

        for i in range(ib, ib + nblock):
            if vstrands[ih]['skip'][i] == -1: continue

            #########################################
            # important convention for loops:
            # For upstream (+z direction), theta & z go like 0, 1, 2 ...
            # For downstream (-z direction), theta & z go in opposite way.
            vstrands[ih]['theta'][i] = {}
            vstrands[ih]['z'][i] = {}
            for il in range(0, vstrands[ih]['loop'][i] + 1): 
                theta = theta0 + theta_per_bp * resi * float(nblock)/(nblock + dn)
                z = z0 + z_per_bp * resi * float(nblock)/(nblock + dn)

                vstrands[ih]['theta'][i][il] = theta
                vstrands[ih]['z'][i][il] = z

                resi += 1


        ib += nblock            

#print(len(vstrands[0]['z']))
ch = 65    # chain starting from A
nseg = 1

####################################################
# Assign sequences to scaffolds & staples 
sys.stderr.write("Start assigning sequences............\n")

for ih in range(len(vstrands)): 
    vstrands[ih]['stap_nt'] = {}
    vstrands[ih]['scaf_nt'] = {}

for ih in range(len(vstrands)): 

    for ib in range(len(vstrands[ih]['stap'])): 
        if csvfile: 
            pih, pib, nih, nib = vstrands[ih]['stap'][ib]
            if (pib == -1 and nib != -1):  # found 5' end
                assign_seq_staple(ih, ib)
        
        if scafseqname: 
            pih, pib, nih, nib = vstrands[ih]['scaf'][ib]
            if (pib == -1 and nib != -1):  # found 5' end
                assign_seq_scaffold(ih, ib)
#print(vstrands[3]['scaf_nt'])            
#print(csv)
#print(get_direction_stap(3))

####################################################
# Check sequence complementarity
sys.stderr.write("Checking sequence complementarity............\n")
#print(vstrands[0]['scaf_nt'])
for ih in range(len(vstrands)): 
    for ib in range(len(vstrands[ih]['stap'])): 
        list_il = list(range(vstrands[ih]['loop'][ib] + 1))
        for il in list_il: 
            
            if ib in vstrands[ih]['stap_nt']: 
                #sys.stderr.write("%s\t%s\n" % (vstrands[ih]['stap_nt'], vstrands[ih]['scaf_nt']))
                if ib in vstrands[ih]['scaf_nt']: 
                    if (wc(vstrands[ih]['stap_nt'][ib][il]) != vstrands[ih]['scaf_nt'][ib][il]): 
                        sys.stderr.write("Warning: scaf=%s stap=%s at ih=%d and ib=%d\n" % (vstrands[ih]['scaf_nt'][ib][il], vstrands[ih]['stap_nt'][ib][il], ih, ib))
                
                else:
                    vstrands[ih]['scaf_nt'][ib][il] = wc(vstrands[ih]['stap_nt'][ib][il])
                
            
            elif ib in vstrands[ih]['scaf_nt']:
                vstrands[ih]['stap_nt'][ib][il] = wc(vstrands[ih]['scaf_nt'][ib][il])
            

#################################################
# Print scaffold and staple sequence for checking
sys.stderr.write("ib      =    ")
for ib in range(len(vstrands[0]['scaf'])): 
    if (ib % nblock) == 0: sys.stderr.write("|") 
    sys.stderr.write(str(ib % nblock))

sys.stderr.write("\n")
for ih in range(len(vstrands)): 
    sys.stderr.write("ih_scaf =%3d " % ih)
    for ib in range(len(vstrands[ih]['scaf'])): 
        if (ib % nblock) == 0: sys.stderr.write("|") 
        list_il = list(range(vstrands[ih]['loop'][ib] + 1))
        for il in list_il: 
            sys.stderr.write(vstrands[ih]['scaf_nt'][ib][il]) if ib in vstrands[ih]['scaf_nt'] else sys.stderr.write(".")
    sys.stderr.write("\n")

    sys.stderr.write("ih_stap =%3d " % ih)
    for ib in range(len(vstrands[ih]['stap'])): 
        if (ib % nblock) == 0: sys.stderr.write("|") 
        list_il = list(range(vstrands[ih]['loop'][ib] + 1))
        for il in list_il: 
            sys.stderr.write(vstrands[ih]['stap_nt'][ib][il]) if ib in vstrands[ih]['stap_nt'] else sys.stderr.write(".")
    
    sys.stderr.write("\n")


####################################################
# Print scaffold
# now supports multiple scaffolds.
sys.stderr.write("Printing scaffold PDB............\n")
global_aid = 0 # atom index

for ih in range(len(vstrands)): 
    vstrands[ih]['stap_ri'] = {}
    vstrands[ih]['scaf_ri'] = {}
    vstrands[ih]['stap_atoms'] = {}
    vstrands[ih]['scaf_atoms'] = {}
    for dir in [1, -1]:
        vstrands[ih]['x%s' % dir] = {}
        vstrands[ih]['y%s' % dir] = {}
        vstrands[ih]['z%s' % dir] = {}


for ih in range(len(vstrands)): 
    for ib in range(len(vstrands[ih]['scaf'])): 
        pih, pib, nih, nib = vstrands[ih]['scaf'][ib]
        if (pib == -1 and nib != -1):  # found 5' end
            trace_scaffold(ih,ib)
            nseg += 1
        

####################################################
# Trace staples 
sys.stderr.write("Printing staple PDB............\n")
for ih in range(len(vstrands)): 
    for ib in range(len(vstrands[ih]['stap'])): 

        pih, pib, nih, nib = vstrands[ih]['stap'][ib]
        if (pib == -1 and nib != -1):  # found 5' end
            trace_staple(ih, ib)
            nseg += 1
    

####################################################
# inter-helical push 
# usage:
# > cat push.helix.for.make_ndx  | make_ndx  -f hextube.pdb -o tmp.ndx
# > awk '/^[0-9]/ && NF == 2 {printf "bond %8d %8d %8d %8d\n", $1-1, $2-1, 1, 30}' tmp.ndx > push.extrabonds
with open('push.helix.for.make_ndx', 'w') as FH:
    FH.write("del 0 - 30\n")
    nnih = {} # nearest neighbor list of helices
    for ih in range(len(vstrands)):
        dir = is_upstream(ih) # +1 for upstream, -1 for downstream.

        if ih not in nnih:
            nnih[ih] = {}

        for ib in range(len(vstrands[ih]['stap'])):
            if (vstrands[ih]['stap'][ib][1] == -1 and vstrands[ih]['stap'][ib][3] == -1):
                continue
            if vstrands[ih]['skip'][ib] == -1:
                continue
    
            stap_pih, stap_pib, stap_nih, stap_nib = vstrands[ih]['stap'][ib]
            scaf_pih, scaf_pib, scaf_nih, scaf_nib = vstrands[ih]['scaf'][ib]
 
            if stap_nih not in nnih[ih]:
                nnih[ih][stap_nih] = {}
            if stap_nih not in nnih:
                nnih[stap_nih] = {}
            if ih not in nnih[stap_nih]:
                nnih[stap_nih][ih] = {}    
            if scaf_nih not in nnih[ih]:
                nnih[ih][scaf_nih] = {}
            if scaf_nih not in nnih:
                nnih[scaf_nih] = {}
            if ih not in nnih[scaf_nih]:
                nnih[scaf_nih][ih] = {}    

            # Number of base-pairs within HJ which would be excluded from extrabonds
            cutOff = 8 if (lattice == "HONEYCOMB") else 11

            if (ih >= 0 and stap_nih >= 0 and ih != stap_nih): 
                #print STDERR "stap cross from H$ih to H$stap_nih at $ib\n";
                for off in range(cutOff): 
                    nnih[ih][stap_nih][ib+off] = 1
                    nnih[stap_nih][ih][ib+off] = 1
                    nnih[ih][stap_nih][ib-off] = 1
                    nnih[stap_nih][ih][ib-off] = 1
                
            
            if (ih >= 0 and scaf_nih >= 0 and ih != scaf_nih): 
                #print STDERR "scaf cross from H$ih to H$scaf_nih at $ib\n";
                for off in range(cutOff): 
                    nnih[ih][scaf_nih][ib+off] = 1
                    nnih[scaf_nih][ih][ib+off] = 1
                    nnih[ih][scaf_nih][ib-off] = 1
                    nnih[scaf_nih][ih][ib-off] = 1
                
    
    for ih in range(len(vstrands)-1): 
        for jh in range(ih+1, len(vstrands)): 
            if not (ih in nnih and jh in nnih[ih]):
                continue
            sys.stderr.write("H%d and H%d are in contact\n" % (ih, jh))
            for ib in range(len(vstrands[ih]['stap'])): 
                ## NOTE: $scaf_ri or $stap_ri can be undefined.
                if not ('scaf_ri' in vstrands[ih] and ib in vstrands[ih]['scaf_ri'] and vstrands[ih]['scaf_ri'][ib] != ''):
                    continue
                if not ('stap_ri' in vstrands[ih] and ib in vstrands[ih]['stap_ri'] and vstrands[ih]['stap_ri'][ib] != ''):
                    continue
                if not ('scaf_ri' in vstrands[jh] and ib in vstrands[jh]['scaf_ri'] and vstrands[jh]['scaf_ri'][ib] != ''):
                    continue
                if not ('stap_ri' in vstrands[jh] and ib in vstrands[jh]['stap_ri'] and vstrands[jh]['stap_ri'][ib] != ''):
                    continue
    
                # skip if close to HJ.
                if (ih in nnih and jh in nnih[ih] and ib in nnih[ih][jh]):
                    continue
    
                scaf_ri = vstrands[ih]['scaf_ri'][ib][0]
                scaf_rj = vstrands[jh]['scaf_ri'][ib][0]
                stap_ri = vstrands[ih]['stap_ri'][ib][0]
                stap_rj = vstrands[jh]['stap_ri'][ib][0]
                FH.write("ri %d %d & a P\n" % (scaf_ri, stap_rj))
                FH.write("ri %d %d & a P\n" % (scaf_rj, stap_ri))
        
        
    
    FH.write("q\n")

    
####################################################
# chickenwire (CW)
 
with open('chickenwire.for.make_ndx', 'w') as FH:

    FH.write("del 0 - 30\n")
    FH.write("!a C1' C2'  C3' C4' O4' C5'  P  O1P  O2P  O3'  O5' H*\n")
    FH.write("name 0 BASE\n")
    # hbond pair
    with open('hbonds.for.make_ndx', 'w') as FHB:
        FHB.write("del 0 - 30\n")
        FHB.write("r DT* DC* & a N3\n")
        FHB.write("r DA* DG* & a N1\n")
        FHB.write("0 | 1\n")
        FHB.write("del 0\n")
        FHB.write("del 0\n")
        FHB.write("name 0 N1N3\n")

        CWatom_num   = 0
        CWatoms      = []  
        CWres        = []
        CWseg        = []
        CWatoms_hash = {}
        CWbonds      = [] 
        CWjuncs      = [] 
        CW_num_cross = {}   # crossover number at ib

        for ih in range(len(vstrands)): 
            dir = is_upstream(ih) # +1 for upstream, -1 for downstream.
            for ib in range(len(vstrands[ih]['stap'])): 
                if ((vstrands[ih]['stap'][ib][1] == -1 and vstrands[ih]['stap'][ib][3] == -1) and
(vstrands[ih]['scaf'][ib][1] == -1 and vstrands[ih]['scaf'][ib][3] == -1)):
                    continue
                if vstrands[ih]['skip'][ib] == -1:
                    continue
                
                (stap_pih,stap_pib,stap_nih,stap_nib) = vstrands[ih]['stap'][ib]
                (scaf_pih,scaf_pib,scaf_nih,scaf_nib) = vstrands[ih]['scaf'][ib]
                
                # Atoms
                loops = list(range(vstrands[ih]['loop'][ib] + 1))
                if dir == -1:
                    loops = loops[::-1]

                #my $inside_loop;
                for il in loops: 
                    CWatom_num += 1
                    ## NOTE: $scaf_ri or $stap_ri can be undefined.
                    scaf_ri = vstrands[ih]['scaf_ri'][ib][il] if ('scaf_ri' in vstrands[ih] and ib in vstrands[ih]['scaf_ri'] and il in vstrands[ih]['scaf_ri'][ib]) else ''
                    stap_ri = vstrands[ih]['stap_ri'][ib][il] if ('stap_ri' in vstrands[ih] and ib in vstrands[ih]['stap_ri'] and il in vstrands[ih]['stap_ri'][ib]) else ''
        
                    FH.write("ri %s %s & 0\n" % (scaf_ri, stap_ri) )
                    sname = "H%d" % ih
                    rname = "B%d" % ib
                    aname = "L%d" % il
                    CWatoms.append(aname)
                    CWres.append(rname)
                    CWseg.append(sname)
                    key = "%s%s%s" % (sname, rname, aname)
                    CWatoms_hash[key] = CWatom_num
                    FH.write("name %d %s\n" % (CWatom_num, key))

                    #################
                    # hbonds
                    if ((scaf_ri != '') and
                        (stap_ri != '')): 
                        #x1 = vstrands[ih][x1][ib][il]
                        #y1 = vstrands[ih][y1][ib][il]
                        #z1 = vstrands[ih][z1][ib][il]
                        #x2 = vstrands[ih]["x-1"][ib][il]
                        #y2 = vstrands[ih]["y-1"][ib][il]
                        #z2 = vstrands[ih]["z-1"][ib][il]
                        #dx = x1-x2
                        #dy = y1-y2
                        #dz = z1-z2
                        #dd = sqrt(dx*dx + dy*dy + dz*dz)
                        #if (dd > 4): 
                        #    sys.stderr.write("%d-%d too far! (%d,%d)-(%d,%d)\n" % (scaf_ri, stap_ri, x1, y1, x2, y2))
                        
                        FHB.write("ri %d & 0\n" % (scaf_ri))
                        FHB.write("ri %d & 0\n" % (stap_ri))
 
                    #################
                # ib-to-ib bonds.
                # only with next one to avoid double counting.
                if (stap_nih != -1 and stap_nib != -1): 
                    (stap_nih,stap_nib) = next_nonskip_ih_ib("stap",stap_nih,stap_nib)
                    txt1 = "H%dB%dL%d" % (ih, ib, loops[-1])
                    txt2 = "H%dB%dL0" % (stap_nih, stap_nib)
                    CWbonds.append([txt1, txt2])
                    #print STDERR "CW: stap cross H$ih B$ib L$loops[$#loops] to H$stap_nih B$stap_nib L0\n";
                    if ib in CW_num_cross:
                        CW_num_cross[ib] += 1
                    else:
                        CW_num_cross[ib] = 0
                
                if (scaf_nih != -1 and scaf_nib != -1): 
                    (scaf_nih, scaf_nib) = next_nonskip_ih_ib("scaf", scaf_nih, scaf_nib)
                    txt1 = "H%dB%dL%d" % (ih, ib, loops[-1])
                    txt2 = "H%dB%dL0" % (scaf_nih, scaf_nib)
                    CWbonds.append([txt1, txt2])
                    #print STDERR "CW: stap cross H$ih B$ib L$loops[$#loops] to H$scaf_nih B$scaf_nib L0\n";
                    if ib in CW_num_cross:
                        CW_num_cross[ib] += 1
                    else:
                        CW_num_cross[ib] = 0

        FH.write("quit")
        FHB.write("quit")

## write CW_num_cross
with open('chickenwire.num.cross', 'w') as FH:
    k = list(CW_num_cross.keys())
    for i in range(max(k)): 
        FH.write("%d %d\n" % (i+1, CW_num_cross[i] if (i in CW_num_cross) else 0))
    

## write psf
with open('chickenwire.psf', 'w') as FPSF: 
    FPSF.write("PSF EXT\n")
    FPSF.write("%10d !NTITLE\n" % 1)
    FPSF.write("* Chickenwire representation made by cadnano2pdb.py EXT\n\n")
    
    FPSF.write("%10d !NATOM\n" % (len(CWatoms)) )
    #fmt02='(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6)'

    index=0
    resid=1
    for a in range(len(CWatoms)): 
        index += 1
        FPSF.write("%10d " % index)
        FPSF.write("%-8s " % CWseg[a])
        FPSF.write("%-8d " % index)
        FPSF.write("%-8s " % CWres[a])
        FPSF.write("%-6s " % CWatoms[a])
        FPSF.write("%6d "  % 6)
        FPSF.write("  0.000000       1.00800")
        FPSF.write("           0   0.00000     -0.301140E-02\n")

    ## replace CWbonds using hash
    CWbonds_num = []
    for i in range(len(CWbonds)): 
        #print STDERR "key=$CWbonds[$i][0] $CWatoms_hash{$CWbonds[$i][0]}\n";
        #print STDERR "key=$CWbonds[$i][1] $CWatoms_hash{$CWbonds[$i][1]}\n";
        if (CWbonds[i][0] not in CWatoms_hash): 
            sys.stderr.write("Undefined %s\n" % CWbonds[i][0])
            continue
        
        if (CWbonds[i][1] not in CWatoms_hash): 
            sys.stderr.write("Undefined %s\n" % CWbonds[i][1])
            continue
        
        bi = CWatoms_hash[ CWbonds[i][0] ]
        bj = CWatoms_hash[ CWbonds[i][1] ]
        if (bi < 1 or bi > len(CWatoms)+1):
            sys.stderr.write("Bond(%s,%s) out of range\n" % (bi,bj)) 
        if (bj < 1 or bj > len(CWatoms)+1):
            sys.stderr.write("Bond(%s,%s) out of range\n" % (bi,bj)) 
        CWbonds_num.append([bi,bj])

    FPSF.write("\n%10d !NBOND: bonds\n" % len(CWbonds_num))
    count = 0
    for i in CWbonds_num:
        for j in i:
            FPSF.write("%10d" % j)
            count+=1
            if count % 8 == 0:
                FPSF.write("\n")
    
    FPSF.write("\n")
    #while (len(CWbonds_num) >= 7): 
    #    #printf STDERR "%10d%10d%10d%10d%10d%10d%10d%10d\n", splice(@CWbonds_num, 0, 8);
    #    FPSF.write("%10d%10d%10d%10d%10d%10d%10d%10d\n", splice(@CWbonds_num, 0, 8))
    #
    #while ($#CWbonds_num >= 0) {
    #    printf FPSF "%10d", shift(@CWbonds_num);
    #}
    #print FPSF "\n";
    
    FPSF.write("\n%10d !NTHETA: angles\n" % 0)
    FPSF.write("\n%10d !NPHI: dihedrals\n" % 0)
    FPSF.write("\n%10d !NIMPHI: impropers\n" % 0)
    FPSF.write("\n%10d !NDON: donors\n" % 0)
    FPSF.write("\n%10d !NACC: acceptors\n" % 0)








