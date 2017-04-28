# Creating mdin files and pbs file for equilibration process according to Beveridge's paper:
# initial minimization --> 100ps heating --> five 50ps-equilibration with 1000 steps of energy minimization

import os
from contextlib import contextmanager

def working_directory(path):
    current_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(current_dir)
        path.__exit__()


def equr_in(filename,num,res):
   """write equilibration mdin file"""
   f = open(filename + ".in", 'w')

   f.write('equilibration with constrain\n&cntrl\n'
           'imin=0, irest=0, ntx=1, \n'
           'ntb=2, pres0=1.0, ntp=1, \n'
           'ntt=1, temp0=300.0, tautp=0.2, \n'
           'ntc=2, ntf=2, cut=10.0, ntr=1, \n'
           'nstlim=25000, dt=0.002, \n'
           'ntpr=5000, ntwx=5000, ntwr=25000 \n/\n'
           'DNA restraint ')
   f.write(str(num) + 'kcal/mol/A\n' + 
           str(num) + '\nRES 1 ' + str(res) + 
           '\nEND\nEND')

def minr_in():
   """write the 1000 steps minimization process mdin file"""
   f = open("min.in", 'w')
   f.write('1000 step minimization\n'
           ' &cntrl\n'
           'imin=1, irest=1, ntx=7, \n'
           'maxcyc=1000, ncyc=500, \n'
           'ntb=1, ntr=0, cut=10.0 \n/')

def heating_in(res):
   """write the heating process mdin file"""
   f = open("heating.in", 'w')

   f.write('100ps heating process with res on DNA\n'
           ' &cntrl\n'
           'imin=0, irest=0, ntx=1, ntb=1, \n'
           'cut=10.0, ntr=1, ntc=2, ntf=2, \n'
           'tempi=0.0, temp0=300.0, ntt=3, gamma_ln=0.1, \n'
           'nstlim=5000, ntwx=5000, ntwr=25000\n/\n'
           'Keep DNA fixed with weak restraints\n'
           '25.0\nRes 1 ')
   f.write(str(res) + '\nEND\nEND')

def inimin_in():
   """write initial minimization mdin file"""
   f = open("inimin.in", 'w')

   f.write('initial minimization process\n &cntrl\n'
           'imin=1, maxcyc=2500, \n'
           'ncyc=1000, ntb=1, ntr=0, cut=10.0\n/')    

def md5ns_in(res):
   """write the 5ns md simulation input file"""
   f = open("5ns.in", 'w')
   f.write('5ns md simulation\n&cntrl\n'
           'imin=0, irest=1, ntx=7, \n'
           'ntb=2, pres0=1.0, ntp=1, \n'
           'ntt=1, temp0=300.0, tautp=5.0, \n'
           'ntc=2, ntf=2,cut=9.0, \n'
           'nstlim=2500000, dt=0.002, '
           'ntpr=20000, ntwx=10000, ntwr=500000\n/')       

def write_pbs():
   """write pbs file for the simulation process"""
   # task name
   taskname = "2x2-bugfix"

   # prmtop file
   prmtop = "2x2.prmtop"

   # inpcrd file
   inpcrd = "2x2.inpcrd"

   # DNA base pair number
   res = 274

   # constraint list
   cons_list = [5.0,4.0,3.0,2.0,1.0,0.5,0.0]

   # directory on server
   ser_dir = "/raid/aslovin/DNA-hinge/2x2-bugfix"

   # directory of amber
   soft_dir = "mpirun -np 12 /home/aslovin/amber14/exe/sander.MPI"

   f = open("p.pbs",'w')

   f.write('#!/bin/sh\n#PBS -N ' + taskname + 
           '\n#PBS -m n\n#PBS -d ' + ser_dir +
           '\n#PBS -l nodes=1:ppn=12\n\n')

   os.mkdir("inimin")
   current_dir = os.getcwd()
   os.chdir('inimin')
   inimin_in()
   os.chdir(current_dir)
   f.write('cd ' + ser_dir + '/inimin\n' + soft_dir + 
           ' -O -i inimin.in -o inimin.out -p ../' + prmtop +
           ' -c ../' + inpcrd + ' -r inimin.rst\n')

   os.mkdir("heating")
   os.chdir('heating')
   heating_in(res)
   os.chdir(current_dir)
   f.write('cd ../heating\n' + soft_dir + 
           ' -O -i heating.in -o heating.out -p ../' + prmtop + 
           ' -c ../inimin/inimin.rst' + ' -ref ../inimin/inimin.rst -r heating.rst -x heating.mdcrd\n')
   
   minr_in()
   for i,r in enumerate(cons_list):
      f.write('cd ../minr' + str(i) + '\n' +  soft_dir +
              ' -O -i ../min.in -o min.out -p ../' + prmtop + ' -r min.rst ')
      if i == 0:
         f.write('-c ../heating/heating.rst\n')
      else:
         f.write('-c ../equr' + str(i-1) + '/equr.rst\n')

      os.mkdir("minr" + str(i))
      os.mkdir("equr" + str(i))
      os.chdir('equr' + str(i))
      equr_in('equr' + str(i),r,res)
      os.chdir(current_dir)

      f.write('cd ../equr' + str(i) + '\n' + soft_dir +
              ' -O -i equr' + str(i) + '.in -o equr' + str(i) + '.out -p ../' + 
              prmtop + ' -c ../minr' + str(i) + '/min.rst -ref ../minr' + str(i) + 
              '/min.rst -r equr.rst -x equr' + str(i) + '.mdcrd\n')

   os.mkdir("5ns")
   os.chdir('5ns')
   md5ns_in(res)
   os.chdir(current_dir)
   f.write('cd ../5ns\n' + soft_dir + ' -p ' + prmtop + ' -O -i 5ns.in -o 5ns.out -c ../equr6/equr.rst -ref ../equr6/equ.rst -r 5ns.rst -x 5ns.mdcrd')
   f.close()



write_pbs()

