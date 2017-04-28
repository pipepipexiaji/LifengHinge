import os
import base
import readers
import numpy as np


l = readers.LorenzoReader("prova.conf", "prova.top")
s = l.get_system()
s.map_nucleotides_to_strands()
system = base.System(s._box, s._time)
init = []

break_index = 0

for strand in s._strands:
    break_index = 0
    for i in range(len(strand._nucleotides)-1):
        dr = strand._nucleotides[i].distance(strand._nucleotides[i+1],box=system._box)
        d = np.sqrt(np.dot(dr,dr))
        if d > 0.7525+0.25:
            init.append([strand._nucleotides[i].index,strand._nucleotides[i+1].index,d])
            #if strand.get_slice(break_index, i).N != 0:
            system.add_strand(strand.get_slice(break_index, i+1), check_overlap=False)
            break_index = i+1
        #if strand.get_slice(break_index, strand._last-strand._first).N != 0:
    system.add_strand(strand.get_slice(break_index, strand._last-strand._first+1),check_overlap=False)
    if break_index == 0 and strand._circular==True:
        system._strands[-1]._circular = True

r = 0
while True:
    ini_ = []
    flag = 1
    for i in range(len(init)):
        if init[i][2] == 0.1:
            ini_.append(init[i])
            continue
        elif init[i][2] < 6.0:
            init[i][2] = 0.1
            ini_.append(init[i])
            flag = 0
        else:
            init[i][2] -= 6.0
            ini_.append(init[i])
            flag = 0
    if flag:
        break

    try:
        os.mkdir("ini%i" % r)
    except:
        pass

    f = open("ini%i/init" % r,"w")
    for i in range(len(ini_)):
        f.write("{\ntype = mutual_trap\nparticle = " +
        str(ini_[i][0]) +
        "\nref_particle = " +
        str(ini_[i][1]) +
        "\nstiff = 1.\nr0 = " +
        str(ini_[i][2]) + "\nPBC = 1\n}\n" +
        "{\ntype = mutual_trap\nref_particle = "+
        str(ini_[i][0]) + "\nparticle = " + str(ini_[i][1]) +
        "\nstiff = 1.\nr0 = " +
        str(ini_[i][2]) + "\nPBC = 1\n}\n")
    r += 1

os.mkdir("ini%i" % r)
f = open("ini%i/init" % r, "w")
for i in range(len(ini_)):
    f.write("{\ntype = mutual_trap\nparticle = " +
        str(ini_[i][0]) +
        "\nref_particle = " +
        str(ini_[i][1]) +
        "\nstiff = 20.\nr0 = " +
       "0.1\nPBC = 1\n}\n" +
        "{\ntype = mutual_trap\nref_particle = "+
        str(ini_[i][0]) + "\nparticle = " + str(ini_[i][1]) +
        "\nstiff = 20.\nr0 = " +
        "0.1\nPBC = 1\n}\n")

system._prepare(None)
system.print_lorenzo_output("ini.conf", "cadnano.top")
