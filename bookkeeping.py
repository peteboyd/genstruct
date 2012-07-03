#!/usr/bin/env python

"""
reads and writes files, sets up logging
"""

import logging
import sys
import textwrap
import copy
from datetime import date
from time import time
from atoms import Atoms
from pyspglib import spglib
from operations import *

xyzbondfmt = "%s%12.5f%12.5f%12.5f " +\
             "atom_vector%12.5f%12.5f%12.5f " +\
             "atom_vector%12.5f%12.5f%12.5f\n"
xyzcellfmt1 = "%s%12.5f%12.5f%12.5f " +\
             "atom_vector%12.5f%12.5f%12.5f\n"
xyzcellfmt2 = "%s%12.5f%12.5f%12.5f " +\
             "atom_vector%12.5f%12.5f%12.5f " +\
             "atom_vector%12.5f%12.5f%12.5f " +\
             "atom_vector%12.5f%12.5f%12.5f\n"
xyzatomfmt = "%s%12.5f%12.5f%12.5f\n"

#pdb format determined from description in pdb document
pdbatmfmt = "%-6s%5i%5s%1s%3s%2s%4i%1s%11.3f%8.3f%8.3f%6.2f%6.2f%12s%2s\n"
pdbcellfmt = "%-6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P1\n"

class Log:

    def __init__(self, file=None):
        if file is None:
            self.file = "genstruct.out"
        else:
            self.file = file
        self.quiet = False
        self.verbose = True 
        self._init_logging()
        # set up writing to file and terminal

    def _init_logging(self):
        if self.quiet:
            stdout_level = logging.ERROR
            file_level = logging.INFO
        elif self.verbose:
            stdout_level = logging.DEBUG
            file_level = logging.DEBUG
        else:
            stdout_level = logging.INFO
            file_level = logging.INFO

        logging.basicConfig(level=file_level,
                            format='[%(asctime)s] %(levelname)s %(message)s',
                            datefmt='%Y%m%d %H:%M:%S',
                            filename=self.file,
                            filemode='a')
       
        logging.addLevelName(10, '--')
        logging.addLevelName(20, '>>')
        logging.addLevelName(30, '**')
        logging.addLevelName(40, '!!')
        logging.addLevelName(50, 'XX')

        console = ColouredConsoleHandler(sys.stdout)
        console.setLevel(stdout_level)
        console.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
        logging.getLogger('').addHandler(console)

class ColouredConsoleHandler(logging.StreamHandler):
    """Makes colourised and wrapped output for the console."""
    def emit(self, record):
        """Colourise and emit a record."""
        # Need to make a actual copy of the record
        # to prevent altering the message for other loggers
        myrecord = copy.copy(record)
        levelno = myrecord.levelno
        if levelno >= 50:  # CRITICAL / FATAL
            front = '\033[30;41m'  # black/red
            text = '\033[30;41m'  # black/red
        elif levelno >= 40:  # ERROR
            front = '\033[30;41m'  # black/red
            text = '\033[1;31m'  # bright red
        elif levelno >= 30:  # WARNING
            front = '\033[30;43m'  # black/yellow
            text = '\033[1;33m'  # bright yellow
        elif levelno >= 20:  # INFO
            front = '\033[30;42m'  # black/green
            text = '\033[1m'  # bright
        elif levelno >= 10:  # DEBUG
            front = '\033[30;46m'  # black/cyan
            text = '\033[0m'  # normal
        else:  # NOTSET and anything else
            front = '\033[0m'  # normal
            text = '\033[0m'  # normal

        myrecord.levelname = '%s%s\033[0m' % (front, myrecord.levelname)
        myrecord.msg = textwrap.fill(
            myrecord.msg, initial_indent=text, width=76,
            subsequent_indent='\033[0m   %s' % text) + '\033[0m'
        logging.StreamHandler.emit(self, myrecord)

class CIF(object):
    """
    Write cif files
    """
    def __init__(self, atoms, connect_table=None):
        #self.tol = 1.267526
        self.tol = 0.01
        self.atoms = atoms
        self.symbols = []
        self.symdata = None
        self.symdata = spglib.get_symmetry_dataset(self.atoms, 
                                                    symprec=self.tol)
        if connect_table is not None:
            self.connect_table = connect_table
        else:
            self.connect_table = None
        
    def add_labels(self, symbols):
        """
        Include numbered labels for each atom.
        """
        atomdic = {}

        labels = []
        for atom in symbols:
            atomdic.setdefault(atom,0)
            atomdic[atom] += 1
            labels.append(atom + str(atomdic[atom]))

        return labels

    def write_cif(self, name=None):

        method = "Genstruct - created by Aliens."

        if name is None:
            filename = "default.cif"
        else:
            filename = name + ".cif"

        # get refined data from symmetry finding
        cell, frac_coords, numbers = spglib.refine_cell(self.atoms,symprec=self.tol)
        
        # determine cell parameters
        cellparams = self.get_cell_params(cell)

        lines = "data_" + name.split('/')[1] + "\n"
        today = date.today()
        prefix = "_audit"
        # creation date (under _audit)
        lines += "%--34s"%(prefix + "_creation_date") + \
                today.strftime("%A %d %B %Y") + "\n"
        # creation method (under _audit)
        lines += "%-34s"%(prefix + "_creation_method") + \
                method + "\n"

        prefix = "_symmetry"
        # space group name (under _symmetry)
        lines += "%-34s"%(prefix + "_space_group_name_H-M") + \
                self.symdata["international"] + "\n"
        # space group number (under _symmetry)
        lines += "%-34s"%(prefix + "_Int_Tables_number") + \
                str(self.symdata["number"]) + "\n"
        # cell setting (under _symmetry)

        lines += "%-34s"%(prefix + "_cell_setting") + \
                cell_setting[self.symdata['number']] + "\n"

        prefix = "_cell"
        # cell parameters (under _cell)
        lines += "%-34s"%(prefix + "_length_a") + \
                "%(a)-7.4f\n"%(cellparams)
        lines += "%-34s"%(prefix + "_length_b") + \
                "%(b)-7.4f\n"%(cellparams)
        lines += "%-34s"%(prefix + "_length_c") + \
                "%(c)-7.4f\n"%(cellparams)
        lines += "%-34s"%(prefix + "_angle_alpha") + \
                "%(alpha)-7.4f\n"%(cellparams)
        lines += "%-34s"%(prefix + "_angle_beta") + \
                "%(beta)-7.4f\n"%(cellparams)
        lines += "%-34s"%(prefix + "_angle_gamma") + \
                "%(gamma)-7.4f\n"%(cellparams)
        # fractional coordinates
        lines += "loop_\n"

        lines += "_atom_site_label\n"
        lines += "_atom_site_type_symbol\n"
        lines += "_atom_site_fract_x\n"
        lines += "_atom_site_fract_y\n"
        lines += "_atom_site_fract_z\n"

        
        equiv_atoms = self.symdata['equivalent_atoms']
        unique_atoms = list(set(equiv_atoms))
        unique_coords = [list(frac_coords[i]) for i in unique_atoms]
        unique_symbols = [self.atoms.symbols[i] for i in unique_atoms]
        unique_labels = self.add_labels(unique_symbols)

        for atom in range(len(unique_atoms)):
            line = [unique_labels[atom], unique_symbols[atom]] + \
                    unique_coords[atom]
            lines += "%-7s%-3s%10.5f%10.5f%10.5f\n"%(tuple(line))
        
        if self.connect_table is not None:
            lines += "loop_\n"
            lines += "_geom_bond_atom_site_label_1\n"
            lines += "_geom_bond_atom_site_label_2\n"
            lines += "_geom_bond_distance\n"
            lines += "_ccdc_geom_bond_type\n"
            keys = self.connect_table.keys()
            keys.sort()
            for bond in keys:
                lines += "%-7s%-7s%10.5f%5s\n"%(unique_labels[bond[0]],
                        unique_labels[bond[1]], self.connect_table[bond], "S")

        ciffile = open(filename, 'w')
        ciffile.writelines(lines)
        ciffile.close()
        return

    def get_cell_params(self, cell):
        """ return alen, blen, clen, alpha, beta, gamma """

        cellparams = {}
        cellparams['a'] = length(cell[0])
        cellparams['b'] = length(cell[1])
        cellparams['c'] = length(cell[2])
        cellparams['alpha'] = calc_angle(cell[1], cell[2])*RAD2DEG
        cellparams['beta'] = calc_angle(cell[0], cell[2])*RAD2DEG
        cellparams['gamma'] = calc_angle(cell[0], cell[1])*RAD2DEG
        return cellparams


class Time:
    """
    Class to time executions
    """
    def __init__(self):
        self.timer = 0.
        self.currtime = time() 

    def timestamp(self):
        currtime = time()
        self.timer = currtime - self.currtime 
        self.currtime = currtime

cell_setting = {
    1   :   "triclinic",       
    2   :   "triclinic",       
    3   :   "monoclinic",      
    4   :   "monoclinic",      
    5   :   "monoclinic",      
    6   :   "monoclinic",      
    7   :   "monoclinic",      
    8   :   "monoclinic",      
    9   :   "monoclinic",      
   10   :   "monoclinic",      
   11   :   "monoclinic",      
   12   :   "monoclinic",      
   13   :   "monoclinic",      
   14   :   "monoclinic",      
   15   :   "monoclinic",      
   16   :   "orthorhombic",      
   17   :   "orthorhombic",      
   18   :   "orthorhombic",      
   19   :   "orthorhombic",      
   20   :   "orthorhombic",      
   21   :   "orthorhombic",      
   22   :   "orthorhombic",      
   23   :   "orthorhombic",      
   24   :   "orthorhombic",      
   25   :   "orthorhombic",      
   26   :   "orthorhombic",      
   27   :   "orthorhombic",      
   28   :   "orthorhombic",      
   29   :   "orthorhombic",      
   30   :   "orthorhombic",      
   31   :   "orthorhombic",      
   32   :   "orthorhombic",      
   33   :   "orthorhombic",      
   34   :   "orthorhombic",      
   35   :   "orthorhombic",      
   36   :   "orthorhombic",      
   37   :   "orthorhombic",      
   38   :   "orthorhombic",      
   39   :   "orthorhombic",      
   40   :   "orthorhombic",      
   41   :   "orthorhombic",      
   42   :   "orthorhombic",      
   43   :   "orthorhombic",      
   44   :   "orthorhombic",      
   45   :   "orthorhombic",      
   46   :   "orthorhombic",      
   47   :   "orthorhombic",      
   48   :   "orthorhombic",      
   49   :   "orthorhombic",      
   50   :   "orthorhombic",      
   51   :   "orthorhombic",      
   52   :   "orthorhombic",      
   53   :   "orthorhombic",      
   54   :   "orthorhombic",      
   55   :   "orthorhombic",      
   56   :   "orthorhombic",      
   57   :   "orthorhombic",      
   58   :   "orthorhombic",      
   59   :   "orthorhombic",      
   60   :   "orthorhombic",      
   61   :   "orthorhombic",      
   62   :   "orthorhombic",      
   63   :   "orthorhombic",      
   64   :   "orthorhombic",      
   65   :   "orthorhombic",      
   66   :   "orthorhombic",      
   67   :   "orthorhombic",      
   68   :   "orthorhombic",      
   69   :   "orthorhombic",      
   70   :   "orthorhombic",      
   71   :   "orthorhombic",      
   72   :   "orthorhombic",      
   73   :   "orthorhombic",      
   74   :   "orthorhombic",      
   75   :   "tetragonal",        
   76   :   "tetragonal",        
   77   :   "tetragonal",        
   78   :   "tetragonal",        
   79   :   "tetragonal",        
   80   :   "tetragonal",        
   81   :   "tetragonal",        
   82   :   "tetragonal",        
   83   :   "tetragonal",        
   84   :   "tetragonal",        
   85   :   "tetragonal",        
   86   :   "tetragonal",        
   86   :   "tetragonal",        
   87   :   "tetragonal",        
   88   :   "tetragonal",        
   89   :   "tetragonal",        
   90   :   "tetragonal",        
   91   :   "tetragonal",        
   92   :   "tetragonal",        
   93   :   "tetragonal",        
   94   :   "tetragonal",        
   95   :   "tetragonal",        
   96   :   "tetragonal",        
   97   :   "tetragonal",        
   98   :   "tetragonal",        
   99   :   "tetragonal",        
  100   :   "tetragonal",        
  101   :   "tetragonal",        
  102   :   "tetragonal",        
  103   :   "tetragonal",        
  104   :   "tetragonal",        
  105   :   "tetragonal",        
  106   :   "tetragonal",        
  107   :   "tetragonal",        
  108   :   "tetragonal",        
  109   :   "tetragonal",        
  110   :   "tetragonal",        
  111   :   "tetragonal",        
  112   :   "tetragonal",        
  113   :   "tetragonal",        
  114   :   "tetragonal",        
  115   :   "tetragonal",        
  116   :   "tetragonal",        
  117   :   "tetragonal",        
  118   :   "tetragonal",        
  119   :   "tetragonal",        
  120   :   "tetragonal",        
  121   :   "tetragonal",        
  122   :   "tetragonal",        
  123   :   "tetragonal",        
  124   :   "tetragonal",        
  125   :   "tetragonal",        
  126   :   "tetragonal",        
  127   :   "tetragonal",        
  128   :   "tetragonal",        
  129   :   "tetragonal",        
  130   :   "tetragonal",        
  131   :   "tetragonal",        
  132   :   "tetragonal",        
  133   :   "tetragonal",        
  134   :   "tetragonal",        
  135   :   "tetragonal",        
  136   :   "tetragonal",        
  137   :   "tetragonal",        
  138   :   "tetragonal",        
  139   :   "tetragonal",        
  140   :   "tetragonal",        
  141   :   "tetragonal",        
  142   :   "tetragonal",        
  143   :   "trigonal",          
  144   :   "trigonal",          
  145   :   "trigonal",          
  146   :   "rhombohedral",   
  147   :   "trigonal",       
  148   :   "rhombohedral",   
  149   :   "trigonal",       
  150   :   "trigonal",       
  151   :   "trigonal",       
  152   :   "trigonal",       
  153   :   "trigonal",       
  154   :   "trigonal",       
  155   :   "rhombohedral",   
  156   :   "trigonal",       
  157   :   "trigonal",       
  158   :   "trigonal",       
  159   :   "trigonal",       
  160   :   "rhombohedral",   
  161   :   "rhombohedral",   
  162   :   "trigonal",       
  163   :   "trigonal",       
  164   :   "trigonal",       
  165   :   "trigonal",       
  166   :   "rhombohedral",   
  167   :   "rhombohedral",   
  168   :   "hexagonal",      
  169   :   "hexagonal",      
  170   :   "hexagonal",      
  171   :   "hexagonal",      
  172   :   "hexagonal",      
  173   :   "hexagonal",      
  174   :   "hexagonal",      
  175   :   "hexagonal",      
  176   :   "hexagonal",      
  177   :   "hexagonal",      
  178   :   "hexagonal",      
  179   :   "hexagonal",      
  180   :   "hexagonal",      
  181   :   "hexagonal",      
  182   :   "hexagonal",      
  183   :   "hexagonal",      
  184   :   "hexagonal",      
  185   :   "hexagonal",      
  186   :   "hexagonal",      
  187   :   "hexagonal",      
  188   :   "hexagonal",      
  189   :   "hexagonal",      
  190   :   "hexagonal",      
  191   :   "hexagonal",      
  192   :   "hexagonal",      
  193   :   "hexagonal",      
  194   :   "hexagonal",      
  195   :   "cubic",          
  196   :   "cubic",          
  197   :   "cubic",          
  198   :   "cubic",          
  199   :   "cubic",          
  200   :   "cubic",          
  201   :   "cubic",          
  202   :   "cubic",          
  203   :   "cubic",          
  204   :   "cubic",          
  205   :   "cubic",          
  206   :   "cubic",          
  207   :   "cubic",          
  208   :   "cubic",          
  209   :   "cubic",          
  210   :   "cubic",          
  211   :   "cubic",          
  212   :   "cubic",          
  213   :   "cubic",          
  214   :   "cubic",          
  215   :   "cubic",          
  216   :   "cubic",          
  217   :   "cubic",          
  218   :   "cubic",          
  219   :   "cubic",          
  220   :   "cubic",          
  221   :   "cubic",          
  222   :   "cubic",          
  223   :   "cubic",          
  224   :   "cubic",          
  225   :   "cubic",          
  226   :   "cubic",          
  227   :   "cubic",          
  228   :   "cubic",          
  229   :   "cubic",          
  230   :   "cubic",          
  } 

# main group values taken from J.Phys.Chem.A, 113, 5811, 2009
# transition metal radii taken from Bondi JPC '64
#Radii = {
#  "H"   :   0.10,
#  "He"  :   0.40,
#  "Li"  :   0.81,
#  "Be"  :   0.53,
#  "B"   :   0.92,
#  "C"   :   0.70,
#  "N"   :   0.55,
#  "O"   :   0.52,
#  "F"   :   0.47,
#  "Ne"  :   0.54,
#  "Na"  :   0.27,
#  "Mg"  :   0.73,
#  "Al"  :   0.84,
#  "Si"  :   0.10,
#  "P"   :   0.80,
#  "S"   :   0.80,
#  "Cl"  :   0.75,
#  "Ar"  :   0.88,
#  "K"   :   0.75,
#  "Ca"  :   0.31,
#  "Ga"  :   0.87,
#  "Ge"  :   0.11,
#  "As"  :   0.85,
#  "Se"  :   0.90,
#  "Br"  :   0.83,
#  "Kr"  :   0.02,
#  "Rb"  :   0.03,
#  "Sr"  :   0.49,
#  "In"  :   0.93,
#  "Sn"  :   0.17,
#  "Sb"  :   0.06,
#  "Te"  :   0.06,
#  "I"   :   0.98,
#  "Xe"  :   0.16,
#  "Cs"  :   0.43,
#  "Ba"  :   0.68,
#  "Tl"  :   0.96,
#  "Pb"  :   0.02,
#  "Bi"  :   0.07,
#  "Po"  :   0.97,
#  "At"  :   0.02,
#  "Rn"  :   0.20,
#  "Fr"  :   0.48,
#  "Ra"  :   0.83,
## Transition metals.  4s means the radius of the 4s shell was used,
## B means taken from Bondi DOI: 10.1246/bcsj.20100166
#  "Sc"  :   0.08,   # 4s
#  "Ti"  :   0.99,   # 4s
#  "V"   :   0.91,   # 4s
#  "Cr"  :   0.92,   # 4s
#  "Mn"  :   0.77,   # 4s
#  "Fe"  :   0.71,   # 4s
#  "Co"  :   0.65,   # 4s
#  "Ni"  :   0.63,   # B
#  "Cu"  :   0.40,   # B
#  "Zn"  :   0.39,   # B
#  "Y"   :   0.23,   # 5s
#  "Zr"  :   0.12,   # 5s
#  "Nb"  :   0.03,   # 5s
#  "Mo"  :   0.95,   # 5s
#  "Tc"  :   0.89,   # 5s
#  "Ru"  :   0.89,   # 5s
#  "Rh"  :   0.86,   # 5s
#  "Pd"  :   0.63,   # B
#  "Ag"  :   0.72,   # B
#  "Cd"  :   0.58    # B
# }

Radii = {
 "H"   :   1.10,
 "He"  :   1.40,
 "Li"  :   1.81,
 "Be"  :   1.53,
 "B"   :   1.92,
 "C"   :   1.70,
 "N"   :   1.55,
 "O"   :   1.52,
 "F"   :   1.47,
 "Ne"  :   1.54,
 "Na"  :   2.27,
 "Mg"  :   1.73,
 "Al"  :   1.84,
 "Si"  :   2.10,
 "P"   :   1.80,
 "S"   :   1.80,
 "Cl"  :   1.75,
 "Ar"  :   1.88,
 "K"   :   2.75,
 "Ca"  :   2.31,
 "Ga"  :   1.87,
 "Ge"  :   2.11,
 "As"  :   1.85,
 "Se"  :   1.90,
 "Br"  :   1.83,
 "Kr"  :   2.02,
 "Rb"  :   3.03,
 "Sr"  :   2.49,
 "In"  :   1.93,
 "Sn"  :   2.17,
 "Sb"  :   2.06,
 "Te"  :   2.06,
 "I"   :   1.98,
 "Xe"  :   2.16,
 "Cs"  :   3.43,
 "Ba"  :   2.68,
 "Tl"  :   1.96,
 "Pb"  :   2.02,
 "Bi"  :   2.07,
 "Po"  :   1.97,
 "At"  :   2.02,
 "Rn"  :   2.20,
 "Fr"  :   3.48,
 "Ra"  :   2.83,
# Transition metals.  4s means the radius of the 4s shell was used,
# B means taken from Bondi DOI: 10.1246/bcsj.20100166
 "Sc"  :   2.08,   # 4s
 "Ti"  :   1.99,   # 4s
 "V"   :   1.91,   # 4s
 "Cr"  :   1.92,   # 4s
 "Mn"  :   1.77,   # 4s
 "Fe"  :   1.71,   # 4s
 "Co"  :   1.65,   # 4s
 "Ni"  :   1.63,   # B
 "Cu"  :   1.40,   # B
 "Zn"  :   1.39,   # B
 "Y"   :   2.23,   # 5s
 "Zr"  :   2.12,   # 5s
 "Nb"  :   2.03,   # 5s
 "Mo"  :   1.95,   # 5s
 "Tc"  :   1.89,   # 5s
 "Ru"  :   1.89,   # 5s
 "Rh"  :   1.86,   # 5s
 "Pd"  :   1.63,   # B
 "Ag"  :   1.72,   # B
 "Cd"  :   1.58    # B
 }
