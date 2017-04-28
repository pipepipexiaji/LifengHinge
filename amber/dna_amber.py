# This code simply ignores "skipped" bases when generating the origami model. The resulting structure may contain highly stretched bonds and extreme local twist in regions where bases have been skipped. It is hoped that minimization and equilibration will resolve these tensions and converge to a stable geometry. Long minimization and equilibration runs may therefore be necessary when a design contains many skipped bases.

from subprocess import call
from contextlib import contextmanager
from shutil import copy
import sys
import os
import json

@contextmanager
def working_directory(path):
    current_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(current_dir)
        

def complement(base):
   """Returns the complement of a base"""
   if base == 'A' or base == 'a':
      cbase = 'T'
   elif base == 'T' or base == 't':
      cbase = 'A'
   elif base == 'C' or base == 'c':
      cbase = 'G'
   elif base == 'G' or base == 'g':
      cbase = 'C'
   else:
      # Change this to change the default identity of 'free' staple strand regions.
      cbase = 'A'
   return cbase


class Design:
   def __init__(self,cfilename,sfilename):
      # First import the caDNAno design
      cfile = open(cfilename)
      raw_data = json.load(cfile)
      cfile.close()

      # Sort the relevant data from the .json file into a nested dictionary structure indexed by helix #
      data = {}
      for i,s in enumerate(raw_data["vstrands"]):
         h = s["num"]
         data[h] = {'row': s["row"], 'col': s["col"], 'scaf': s["scaf"], 'stap': s["stap"], 'skip': s["skip"], 'base': [0 for _ in range(len(s["scaf"]))]}
      
      self.data = data
      
      # Then import the base sequence for the scaffold strand
      self.sfile = open(sfilename)
   
   def find_strands(self,h_input,z_input):
      # Create and populate the scaffold strand
      self.scaf = Scaffold(self.data,h_input,z_input,self.sfile)
      self.scaf.find_start()
      self.scaf.populate()
      
      # Create an empty list of staple strands
      self.stap_list = []
      
      # Scan through the data, creating and populating each staple strand in turn
      for h in self.data.keys():
         for z,node in enumerate(self.data[h]["stap"]):
            if (node != [-1]*4) and (self.data[h]["skip"][z] == 0):
               # The Staple.populate method sets the "skip" field to -1 for every node belonging to that staple strand...
               self.stap_list.append( Staple(self.data,h,z) )
               self.stap_list[-1].find_start()
               self.stap_list[-1].populate()
      
      # We no longer need the scaffold sequence file for anything; close it
      self.sfile.close()
               
   def fix_segments(self):
      # Figures out if each strand segment should run in the positive or negative z-direction, and also recalculates their positions relative to (x_min,y_min,z_min)
      scaf_seg_list = [seg for _,seg in enumerate(self.scaf.seg_list)]
      stap_seg_list = [seg for _,strand in enumerate(self.stap_list) for _,seg in enumerate(strand.seg_list)]
      
      # Find the minimum x,y and z values for any strand segment
      x_min,y_min,z_min = ( min(r) for r in zip(*[(seg.x, seg.y, seg.z) for seg in scaf_seg_list + stap_seg_list]) )
      
      # Scaffold strand segments run in the positive-z (sense) direction if the sum of their x and y-coordinates is even
      for seg in scaf_seg_list:
         seg.sense = (seg.x + seg.y) % 2 == 0
         seg.x -= x_min
         seg.y -= y_min
         seg.z -= z_min
      
      # Staple strand segments run in the positive-z (sense) direction if the sum of their x and y-coordinates is odd
      for seg in stap_seg_list:
         seg.sense = (seg.x + seg.y) % 2 != 0
         seg.x -= x_min
         seg.y -= y_min
         seg.z -= z_min
   
   def export(self):
      print "\nScaffold strand " + ("(linear)", "(loop)")[self.scaf.loop] + ":"
      for i,seg in enumerate(self.scaf.seg_list):
         print "  Segment", i, ": x =", seg.x, ", y =", seg.y, ", z =", seg.z, ", seq =", seg.seq, ", sense =", seg.sense
      
      for j,stap in enumerate(self.stap_list):
         print "\nStaple strand #", j, ("(linear)", "(loop)")[stap.loop] + ":"
         for i,seg in enumerate(stap.seg_list):
            print "  Segment", i, ": x =", seg.x, ", y =", seg.y, ", z =", seg.z, ", seq =", seg.seq, ", sense =", seg.sense
   
   def nab(self,temp_dir):
      # Create NAB files in the temporary directory and run them to obtain initial PDB files for each strand
      self.charge_count = 0
      with working_directory(temp_dir):
         self.scaf.export_nab('scaf')
         self.charge_count += self.scaf.base_count
         call("$AMBERHOME/bin/nab scaf.nab", shell=True)
         call("./a.out", shell=True)
         for i,stap in enumerate(self.stap_list):
            strand_name = 'stap' + str(i)
            stap.export_nab(strand_name)
            self.charge_count += stap.base_count
            call("nab " + strand_name + ".nab", shell=True)
            call("./a.out", shell=True)
         call("rm a.out", shell=True)
   
   def bare_amberparm(self,top_file,temp_dir):
      # Create initial amber parameter file for the bare (unsolvated, no salt ions) DNA origami
      with working_directory(temp_dir):
         self.export_amber_bare('make_amber_bare',top_file)
         call("$AMBERHOME/bin/tleap -f make_amber_bare.script", shell=True)
      
   def solvate_amber(self,top_file,temp_dir):
      #Creat prmtop and inpcrd files for the neutralized DNA origami
      with working_directory(temp_dir):
         self.export_amber_solvate('make_amber_solvate',top_file)
         call("$AMBERHOME/bin/tleap -f make_amber_solvate.script", shell=True)
      
   def export_amber_bare(self,filename,top_file):
      # top_file specify the parm file used in tleap, like $AMBERHOME/dat/leap/cmd/oldff/leaprc.ff99bsc0
      f = open(filename + ".script", 'w')
      
      # Load scaffold
      f.write('source ' + top_file+
              '\nmols = loadpdb scaf.pdb\n')
      # Creat a string for the command combining scaffold and staples together
      self.combine_list='mols'
   
      # Load staples
      for i,strand in enumerate(self.stap_list):
         f.write('mol' + str(i) + ' = loadpdb stap' + str(i) + '.pdb\n')
         self.combine_list += ' mol' + str(i)
      
      # Combine the strands together
      f.write('mol = combine { ' + self.combine_list + ' }\n' +
              'savePdb mol bare_dna.pdb' +
              '\nsaveamberparm mol bare_dna.prmtop bare_dna.inpcrd\n' +
              'quit')
      
      f.close()
      
   def export_amber_solvate(self,filename,top_file):
      f = open(filename + ".script", 'w')
   
      f.write('source ' + top_file+
              '\nmols = loadpdb scaf.pdb\n')
      # Creat a string for the command combining scaffold and staples together
      self.combine_list='mols'
   
      # Load staples
      for i,strand in enumerate(self.stap_list):
         f.write('mol' + str(i) + ' = loadpdb stap' + str(i) + '.pdb\n')
         self.combine_list += ' mol' + str(i)

      # Combine the strands together and solvate it after neutralization
      f.write('mol = combine { ' + self.combine_list + ' }\n'+
              'addions mol Na+ 0\n'+
              'solvatebox mol TIP3PBOX 10.0\n'+
              'savePdb mol solvate_dna.pdb\n'+
              'saveamberparm mol solvate_dna.prmtop solvate_dna.inpcrd\n'+
              'quit')
   
   def autoionize(self,temp_dir,top_file):
      # add extra ions
      with working_directory(temp_dir):
         f = open("autoionize.script", 'w')
         f.write('source ' + top_file +
                 '\nmol = loadpdb solvate_dna.pdb\n' +
                 'addionsRand mol Na+ ' + str(self.getions(0.2)) + ' Cl- ' + str(self.getions(0.2)) +
                 '\nsaveamberparm mol dna_ionized.prmtop dna_ionized.inpcrd\n' +
                 'quit')
         f.close()
         call("$AMBERHOME/bin/tleap -f autoionize.script", shell=True)
      
   def getions(self,concentration): 
      # calculate the number of ions to get expected concentration(M)
      f = open("solvate_dna.pdb", 'r')
      self.wat = 0
      for line in f:
         if 'WAT' in line:
            self.wat += 1
      print self.wat
      # Calculate the number of ions (according to the autoionize.tcl script in vmd)
      self.ions = int(0.5 + 0.0187 * concentration * self.wat)
      
      # Return the number of ions
      return self.ions     
 
class Segment:
   """(One strand of) a double helix of DNA with a defined sequence and location in space"""
   def __init__(self,x,y,z,first_base):
      # Initialize a new helix by specifying the coordinates and identity of its first base
      # Helices are oriented along the z-axis, so the (x,y)-coordinates of the helix will never change
      self.x = x
      self.y = y
      self.z = z
      self.seq = first_base
      # TO_DO Check that the coordinates are positive numbers and that input bases are 'A', 'C', 'G' or 'T'. Use assertions.
   
   def add_base(self,z,base):
      # Adds a base to the helix at the specified z-position
      # This may change the z-position of the helix, which is defined as the smallest z-position of any of its constituent bases
      self.seq += base
      self.z = min(self.z, z)
      


class Strand:
   """A continuous covalently bonded sequence of nucleotides arranged in space, represented as an ordered list of linked segments"""
   def __init__(self,data,h_init,z_init):
      self.data = data
      self.h_start = h_init
      self.z_start = z_init
      self.seg_list = []
   
   def find_start(self):
      # Finds the position of the first nucleotide in the strand by walking backwards from an initial guess
      h = self.h_start
      z = self.z_start
      while True:
         # If the previous base along the strand is non-existent, we have reached the start!
         if -1 in self.data[h][self.type][z][0:2]:
            self.loop = False
            self.h_start = h
            self.z_start = z
            break
         # If we wind up back where we started, the strand must be a loop, so we can use the original starting point...
         elif self.data[h][self.type][z][0:2] == [self.h_start,self.z_start]:
            self.loop = True
            break
         # Otherwise update the coordinates to that of the previous base and iterate
         else:
            h, z = self.data[h][self.type][z][0], self.data[h][self.type][z][1]
   
   def populate(self):
      # Start building the first strand segment at the specified starting position
      h = self.h_start
      z = self.z_start
      new_seg = True
      # For each base along the strand,
      while True:
   
         # If the base should not be skipped,
         if self.data[h]['skip'][z] == 0:

            # Read in the identity of the base
            base = self.get_base(h,z)
            
            # Start a new segment if necessary (containing the current base)
            if new_seg == True:
               self.seg_list.append( Segment(self.data[h]['col'], self.data[h]['row'], z, base) )
            # Otherwise add to the current segment
            else:
               self.seg_list[-1].add_base(z, base)
            
            # If the next base along the strand is non-existent (in the case of a linear strand), or if the next base is the original one (the strand is a loop), we are done!
            if (-1 in self.data[h][self.type][z][2:4]) or (self.data[h][self.type][z][2:4] == [self.h_start,self.z_start]):
               break
            
            # If the next base involves crossing to another helix, indicate that a new segment should be started
            if self.data[h][self.type][z][2] == h:
               new_seg = False
            else:
               new_seg = True
      
         # If the base should be skipped, ignore it and search for the next base along the strand
         else:
            
            # If the next base along the strand is non-existent (in the case of a linear strand), or if the next base is the original one (the strand is a loop), we are done!
            if (-1 in self.data[h][self.type][z][2:4]) or (self.data[h][self.type][z][2:4] == [self.h_start,self.z_start]):
               break
            
            # Always start a new segment after skipping one or more bases
            new_seg = True
            
         # Update the coordinates to that of the next base and iterate
         h, z = self.data[h][self.type][z][2], self.data[h][self.type][z][3]
   
   def export_nab(self,filename):
      f = open(filename + ".nab", 'w')
      
      # Declarations
      f.write('molecule m, m_seg;\n' +
              'matrix rotate, move;\n' +
              'string seq, cseq;\n')
      
      # Create main strand molecule
      f.write('\n// Create molecule to hold the main strand\n' +
              'm = newmolecule();\n' +
              'addstrand(m, "str");\n')
      
      self.base_count = 0
      for i,seg in enumerate(self.seg_list):
         # Create each segment in turn
         if seg.sense == False:
            seg.seq = seg.seq[::-1]
         f.write('\n// Build segment ' + str(i) + '\n' +
                 ('cseq','seq')[seg.sense] + ' = "' + seg.seq.lower() + '";\n' +
                 ('seq','cseq')[seg.sense] + ' = wc_complement( ' + ('cseq','seq')[seg.sense] + ', "", "dna" );\n')
         # The "endstr" string specifies whether each segment should be capped at either end (wherever the strand begins/ends)
         endstr = ''
         if not self.loop:
            if i == 0:
               endstr += ('a5','s5')[seg.sense]
            if i == len(self.seg_list) - 1:
               endstr += ('a3','s3')[seg.sense]
         # This line actually creates the segment, based on the values of various helical parameters
         f.write('m_seg = wc_helix( seq, "", "dna", cseq, "", "dna", 2.25, 0.0, 33.75, 3.38, "' + endstr + '" );\n')
         # TO_DO above and below, change the helical parameters to variables that can be set centrally
         
         # Now rotate and move the strand segment so it is positioned correctly
         f.write('\n// Position segment ' + str(i) + '\n' +
                 'rotate = newtransform( 0.0, 0.0, 0.0, 0.0, 0.0, ' + str(55. + 33.75*(seg.z-1) ) + ' );\n' +
                 'move = newtransform( ' + ', '.join( [str(21.*seg.x), str(21.*seg.y), str(3.38*seg.z)] ) + ', 0.0, 0.0, 0.0 );\n' +
                 'transformmol( rotate, m_seg, NULL );\n' +
                 'transformmol( move, m_seg, NULL );\n')
         
         # Add the relevent strand to the main molecule, and then delete the segment
         f.write('\n// Add segment ' + str(i) + ' to the main strand\n' +
                 'mergestr(m, "str", "last", m_seg, "' + ('anti','sense')[seg.sense] + '", "first");\n')
         if i > 0:
            f.write('connectres(m, "str", ' + str(self.base_count) + ', "O3\'", ' + str(self.base_count+1) + ', "P");\n')
         f.write('freemolecule( m_seg );\n')
         
         # Increase the length of the strand by the length of the segment added
         self.base_count += len(seg.seq)
      
      # Close the loop if the strand is indeed a loop
      if self.loop:
         f.write('\n// Close the loop of the main strand\n' +
                 'connectres(m, "str", ' + str(self.base_count) + ', "O3\'", 1, "P");\n')
      # Finally, output the strand to a pdb file
      f.write('\n// Write the completed main strand to a PDB file\n' +
              'putpdb( "' + filename + '.pdb", m, "-wwpdb" );')
      
      f.close()


class Scaffold(Strand):
   def __init__(self,data,h_init,z_init,seq_file):
      Strand.__init__(self,data,h_init,z_init)
      self.seq_file = seq_file
      self.type = 'scaf'
      self.check_start_coords()
   
   def check_start_coords(self):
      # First check that the specified starting helix actually exists
      if self.h_start not in self.data.keys():
         print "\nSpecified starting helix not found! Exiting...\n"
         sys.exit()

      # Then check that the specified initial z-position is valid
      if (self.data[self.h_start][self.type][self.z_start] == [-1]*4) or (self.data[self.h_start]['skip'][self.z_start] == -1):
         print "\nInvalid scaffold starting location! Exiting...\n"
         sys.exit()
      
   def get_base(self,h,z):
      # Read in the identity of the next base in the scaffold sequence file and assign it to position (h,z)
      base = self.seq_file.read(1)
      self.data[h]['base'][z] = base
      return base


class Staple(Strand):
   def __init__(self,data,h_init,z_init):
      Strand.__init__(self,data,h_init,z_init)
      self.type = 'stap'
   
   def get_base(self,h,z):
      # Mark that we have visited the position (h,z)
      self.data[h]['skip'][z] = -1
      # Return the identity of the base complementary to the one present at position (h,z) in the scaffold strand
      return complement(self.data[h]['base'][z])


def convert_origami():
   """Top-level function to convert caDNAno design to model"""
   
   # cfile is the name of your caDNAno .json file
   cfilename = '2x2_b_32bp.json'
   
   # sfile is the name of the text file that contains the scaffold sequence
   sfilename = 'p7560.txt'
   
   # temp_dir is the name of the temporary directory (within the current active directory) in which intermediate NAB, PDB and other files will be built
   temp_dir = 'beam32'

   # At what location should the scaffold sequence start? [h_in, z_in] is the [helix #, z-position] double that appears in the bottom-left corner of your caDNAno window when you hover the mouse over the desired base pair. For instance, if you see "0[20]", h_in = 0 and z_in = 20
   h_in = 1
   z_in = 26
   
   # Topology file
   top_file = "leaprc.ff99bsc0"
   
   
   d = Design(cfilename,sfilename)
   d.find_strands(h_in,z_in)
   d.fix_segments()
   d.export()
   
   os.mkdir(temp_dir)
   d.nab(temp_dir)
   d.bare_amberparm(top_file,temp_dir)
   #d.solvate_amber(top_file,temp_dir)
   #d.autoionize(temp_dir,top_file)


# Now we use the gathered sequence and position data to generate the scaffold 

# Figure out how this helix is oriented. If (row + col) is even, then the scaffold strand runs in the positive z-direction ("sense"), while staple strands run in the negative z-direction ("anti"). If (row + col) is odd, the directions and identities of the scaffold and staple strands are reversed.


# Set the amber environment variables
keys = os.environ.keys()
from re import search
for key in keys:
   if not search("AMBERHOME", key):
      os.environ["AMBERHOME"]="/home/hy/amber14"

if __name__ == "__main__":
    convert_origami()
