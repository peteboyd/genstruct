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
    def __init__(self, atoms=None):
        #self.tol = 1.267526
        self.tol = 1.3
        self.atoms = atoms
        self.symbols = []
        self.symdata = None
        if atoms:
            # populate all parameters
            self.symdata = spglib.get_symmetry_dataset(self.atoms, 
                                                    symprec=self.tol)
        else:
            print "Need atoms"

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
        cell, frac_coords, numbers = spglib.refine_cell(self.atoms,
                                                    symprec=self.tol)
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
