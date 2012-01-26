#!/usr/bin/env python

"""
GenStruct -- Generation of Structures for fapping.
"""

__version__ = "$Revision$"

import subprocess
import textwrap
import sys
import math
import re
import numpy as np
import copy
from numpy import array
from elements import WEIGHT, ATOMIC_NUMBER
import sample
from bookkeeping import Log, Time
from random import random, uniform, randrange, choice
from datetime import date
from logging import warning, debug, error, info, critical

# Constants

RAD2DEG = 180./np.pi

class Generate(object):
    """
    This will contain all the wrappers for iterative and non iterative
    structure generations.
    """
   
    def __init__(self):
        
        # moflib is the list containing MOF structures for extensive
        # branching
        self.moflib = []
        # string history stores how the MOF was built
        self.stringhist = []
        # max number of sbus in a structure
        self.nsbumax = 99

    def random_insert(self, isbu, dataset, struct):
        """
        Randomly places a SBU anywhere
        """
        struct.mof.append(copy.deepcopy(dataset[isbu]))
        struct.connectivity.append([None]*len(dataset[isbu].connectpoints))
        # TODO(pboyd): randomly perturb the coordinates

    def branched_generation(self):
        dataset = [SBU("linear Fe", xsdfile="linear_Fe"),
                   SBU("tetrahedral", xsdfile="tetrahedral", ismetal=False)]
        structcount = 0
        # initialize timing
        stopwatch = Time()
        # Start the string off at "0-0-0-0-0"
        string = "0-0-0-0-0"
        self.stringhist = [string]
        maxstring = self.setmax(dataset)
        # start the timer
        stopwatch.timestamp()
        # only one structure generated (a "leaf")
        struct = Structure(dataset)
        self.random_insert(0, dataset, struct)
        done = False
        while not done:
            stringtype = self.valid_string(string, struct, dataset)
            ints = [int(i) for i in string.split("-")]
            if stringtype == "True":
                info("Applying %s"%(string))
                debug("Number of SBUs in structure: %i"%
                        (len(struct.mof)))
                sbu2 = len(struct.mof)
                struct.apply_string(string)
                if struct.bad_addition(string, struct):
                    warning("SBU overlap occured")
                    struct.mof.pop()
                else:
                    struct.join_sbus(ints[0], ints[1], sbu2, ints[3], False)
                    struct.sbu_check(sbu2)
                    self.stringhist.append(string)
                    
                struct.pbc.complete_box()
                if struct.saturated():
                    struct.final_coords()
                    struct.mof_reorient()
                    struct.xyz_debug()
                    dump = struct.coordinate_dump()
                    write_pdb("final", dump[0], dump[1], struct.pbc.acell)
                    write_xyz("final", dump[0], dump[1], 
                          struct.pbc.cell, struct.pbc.origins)
                    structcount += 1
                
                if structcount == 1:
                    done = True
                    # stop the timer
                    stopwatch.timestamp()
                    info("Structure generated! Timing reports "+
                         "%6.4f seconds."%(stopwatch.timer))
            elif stringtype == "Backup":
                # back to the previous string
                string = self.stringhist[ints[0]]
                self.readjust(struct, string, dataset)

            string = self.iterate_string(string, maxstring)

    def readjust(self, struct, string, dataset):
        """
        Restart the structure generation up to the point where
        the string was successful
        """
        ints = [int(i) for i in string.split("-")]

        sbu0 = copy.deepcopy(struct.mof[0])
        # Reset struct to a blank slate
        struct.reset() 
        struct.mof.append(sbu0)
        for oldstring in self.stringhist[:ints[0]]:
            newsbu = len(struct.mof)
            oldints = oldstring.split("-")
            # Rebuild with the good strings
            struct.apply_string(oldstring)
            struct.join_sbus(oldints[0], oldints[1], 
                             newsbu, oldints[3], False)
            struct.sbu_check(newsbu)


    def valid_string(self, string, struct, database):
        """
        No fucking idea how this will work
        """
        ints = [int(i) for i in string.split("-")]
        # all possibilities are exhausted without joining two SBUs
        # then go back to the previous SBU and change it.
        # TODO(pboyd): this needs to be hashed out.
        if ints[0] >= len(struct.mof):
            return "Backup"

        if ints[1] >= len(struct.mof[ints[0]].connectpoints):
            bondcount = [i for i in struct.connectivity[ints[0]] if
                         i is not None]
            if len(bondcount) <= 1:
                return "Backup"

            return "False"

        if ints[2] >= len(database):
            return "False"

        if ints[3] >= len(database[ints[2]].connectpoints):
            return "False"
        
        if struct.connectivity[ints[0]][ints[1]] is not None:
            return "False"

        if struct.mof[ints[0]].ismetal == database[ints[2]].ismetal:
            return "False"

        if struct.mof[ints[0]].ismetal:
            if ints[4] >= len(database[ints[2]].connectangles[ints[3]]):
                return "False"
        else:
            if ints[4] >= len(struct.mof[ints[0]].connectangles[ints[1]]):
                return "False"

        return "True"

    def iterate_string(self, string, maxstring):
        """
        iterate the string according to the max
        """
        liststring = [int(i) for i in string.split("-")]

        for i, j in enumerate(reversed(liststring)):
            if j == maxstring[-i-1]:
                liststring[-i-1] = 0
            else:
                liststring[-i-1] += 1
                return "%i-%i-%i-%i-%i"%(tuple(liststring))

    def setmax(self, dataset):
        """
        Set maximum string
        """
        maxangle, maxbond = 0, 0
        maxsbu = len(dataset)
        for sbu in dataset:
            anglelen = max([len(i) for i in sbu.connectangles])
            numbonds = len(sbu.connectpoints)
            maxangle = max(maxangle, anglelen)
            maxbond = max(maxbond, numbonds)

        return([self.nsbumax, maxbond, maxsbu, maxbond, maxangle])

    def exhaustive_generation(self):
        """
        generates MOFs via string iteration
        """
        dataset = [SBU("Copper paddlewheel", xsdfile="Cu_paddlewheel"),
                   SBU("tetrahedral", xsdfile="tetrahedral", ismetal=False)]

        # call the Structure class with a set of SBU's and rules to
        # attempt to build the MOF
        self.moflib.append(Structure(dataset))
        # Start the string off at "0-0-0-0-0"
        self.moflib[0].string = "0-0-0-0-0"
        self.random_insert(0, dataset, self.moflib[0])
        done = False
        stringcount = 0
        while not done:
            stringcount+=1
            newlist = []
            for ind, struct in enumerate(self.moflib):
                newstruct = self.pythonsuxcopy(struct)
                if struct.try_iterate_string(struct.string):
                    dir = [int(i) for i in struct.string.split("-")]
                    info("Applying %s"%(struct.string))
                    debug("Number of SBUs in structure: %i"%(len(struct.mof)))
                    debug("Size of MOF array: %i"%(len(self.moflib)))
                    sbu2 = len(newstruct.mof)
                    newstruct.apply_string(struct.string)
                    if newstruct.bad_addition(struct.string, struct):
                        warning("SBU overlap occured")
                    else:
                        newstruct.join_sbus(dir[0], dir[1], sbu2, dir[3], False)
                        newstruct.sbu_check(sbu2)
                        newlist.append(newstruct)
                
                newstruct.pbc.complete_box()
                if newstruct.saturated():
                    newstruct.final_coords()
                    newstruct.mof_reorient()
                    newstruct.xyz_debug()
                    dump = newstruct.coordinate_dump()
                    write_pdb("final", dump[0], dump[1], newstruct.pbc.acell)
                    write_xyz("final", dump[0], dump[1], 
                          newstruct.pbc.cell, newstruct.pbc.origins)
                    done = True
            [self.moflib.append(i) for i in newlist]
    def pythonsuxcopy(self, struct):
        """
        module for copying Structure class which contains SBU class
        and the Cell class.  copy.deepcopy fails for these class - 
        within - class copies
        """
        dataset = copy.deepcopy(struct.sbu_array)
        new = Structure(dataset)
        new.connectivity = copy.deepcopy(struct.connectivity)
        new.pbc = copy.deepcopy(struct.pbc)
        new.mof = copy.deepcopy(struct.mof)
        new.string = copy.deepcopy(struct.string)
        new.sbu_array = copy.deepcopy(struct.sbu_array)
        return new

class SBU(object):
    """
    """

    def __init__(self, type=None, atomlabel=None, coordinates=None,
                 pdbfile=None, xsdfile=None, ismetal=None):
        # index is used to keep track of the SBU's order when building
        # new structures.
        self.index = None
        if type is not None:
            self.type = type
        else: 
            self.type = "Null"
        if atomlabel is not None and coordinates is not None:
            self.atomlabel = atomlabel
            self.coordinates = array(coordinates)
        else:
            self.atomlabel = []
            self.coordinates = []
        if xsdfile is not None:
            self.from_xsd(xsdfile)
        if pdbfile is not None:
            self.from_pdb(pdbfile)
        
        self.mass = 0. 
        # TODO(pboyd): change self.connect to something else
        # this is just the intramolecular connectivity matrix
        self.connect = []
        self.nbonds = [] 
        self.bonding = []
        self.unsat = []
        # vectors of connection points for SBU
        self.connect_vector = []
        # terminal points of connecting vectors for SBU
        self.connectpoints = []
        # angle vectors used to align SBUs once bonded together
        self.anglevect = []
        # angles for possible rotations
        self.connectangles = []

        # Default metal
        if ismetal is not None:
            self.ismetal = ismetal
        else:
            self.ismetal = True

        # populate arrays
        # structcalc is a tuple of important connectivity information
        # of the linker
        structcalc = coord_matrix(self.atomlabel, self.coordinates)
        self.connect = structcalc[0]
        self.bonding = structcalc[1]
        self.nbonds = structcalc[2]

        # calculate centre of mass      
        self.COM = self.centre_of_mass()

        # shift all coordinates back by the centre of mass
        self.COM_shift()

        # check unsaturated bonds for possible connect_vector
        self.check_unsat()

        # add linking points to the SBU.
        # TODO(pboyd): add possible linking types
        self.add_connect_vector()

    def from_pdb(self, filename):
        """Reads atom coordinates of an SBU from a pdb file"""

        filename = filename + ".pdb"
        pdbfile = open(filename)
        pdblines = pdbfile.readlines()
        for line in pdblines:
            if line.lower().startswith('atom'):
                self.atomlabel.append(line[12:16].strip())
                self.coordinates.append(array([float(line[30:38]), \
                    float(line[38:46]), float(line[47:54])]))
        pdbfile.close()

    def from_xsd(self, filename):
        """
        Reads atom coordinates of an SBU from a Materials Studio .xsd
        file.
        """
        filename = filename + ".xsd"
        xsdfile = open(filename)
        xsdlines = xsdfile.readlines()
        xsdfile.close()

        for line in xsdlines:
            line = line.strip()
            line = re.split("[<>=\"]", line)
            line = [x for x in line if x]
            if line[0].strip().lower() == "atom3d id":
                for index, entry in enumerate(line):
                    if entry.strip().lower() == "name":
                        self.atomlabel.append(line[index+1].strip())
                    elif entry.strip().lower() == "xyz":
                        coord = line[index+1].strip().split(",")
                        self.coordinates.append([float(x) for x in coord])

        self.coordinates = array(self.coordinates)
    def check_unsat(self):
        """
        checking for unsaturated bonds
        """
        unsat = []
        # need to improve robustness of this code
        for i in range(len(self.nbonds)):
            if ((self.nbonds[i] <= 2)and(self.atomlabel[i] != "H")and
                (self.atomlabel[i] != "X")and(self.atomlabel[i] != "Y")
                and(self.atomlabel[i]!="Z")):
                unsat.append(i)
        self.unsat = unsat

    def add_connect_vector(self):
        """
        checks for certain atoms and adds topological connecting
        points.  X = the connection point, Y = the vector designed
        to orient the bond, Z = the bond vector(Z-X)
        """

        purgeatoms = []
        for indx, atom in enumerate(self.atomlabel):
            # test for "X" atoms (considered a linking point)
            if atom == "X":
                # there should only be a "Y" and a "Z" bonding to an "X"
                # TODO(pboyd): add error checking if Y and Z not on X
                # also, if no X atom
                for i in self.bonding[indx]:
                    if self.atomlabel[i] == "Y":
                        # angle vector required to align SBUs
                        self.anglevect.append(self.coordinates[i] - 
                                self.coordinates[indx])
                        purgeatoms.append(i)
                    elif self.atomlabel[i] == "Z":
                        vect = (self.coordinates[i] - 
                                self.coordinates[indx])
                        self.connect_vector.append(vect)
                        purgeatoms.append(i)
                pt = self.coordinates[indx]
                self.connectpoints.append(pt)
                self.connectangles.append([0.])
                purgeatoms.append(indx)

        # purge any coordinates in self.coordinates belonging to "X"'s
        # and "Y"'s as they are now part of self.connect_vector 
        # and self.connectpoints and self.anglevect
        self.purge(purgeatoms)

        if len(self.anglevect) == 0:
            sys.exit(0)
        self.connect_vector = array(self.connect_vector)
        self.connectpoints = array(self.connectpoints)
        self.anglevect = array(self.anglevect)

    def purge(self, atoms):
        """Remove entries for X and Y atoms in the atomic coordinates"""
        atoms.sort()
        for patm in reversed(atoms):
            self.coordinates = np.delete(self.coordinates, patm, 0) 
            self.atomlabel.pop(patm)
            self.nbonds.pop(patm)
            self.bonding.pop(patm)

        self.connect = coord_matrix(self.atomlabel, self.coordinates)[0]

    def carboxylate_test(self, atom):
        """Simple test for a carboxylate group"""
        oxcount = 0
        bondvector =[]

        for iatm in self.bonding[atom]:
            if (self.atomlabel[iatm] == "O") and (iatm in self.unsat):
                oxcount += 1
                bondvector.append(self.coordinates[iatm])

        if oxcount == 2:
            return True, bondvector
        else:
            return False, None


    def centre_of_mass(self):
        """
        calculates the centre of mass:
        sum(mass*coordinate) / sum(masses)
        
        or top / bottom

        In this calculation the bottom value is also stored
        as the self.mass
        """
        top = 0.
        bottom = 0.

        for i in range(len(self.atomlabel)):
            top += self.coordinates[i] * WEIGHT[self.atomlabel[i]]
            bottom += WEIGHT[self.atomlabel[i]]
        self.mass = bottom
        return top / bottom

    def COM_shift(self):
        """
        Shifts all the coordinates such that the centre of mass (COM)
        is at the origin.
        This will make it easier for geometry shifts and rotations
        """
        for icoord, xyz in enumerate(self.coordinates):
            self.coordinates[icoord] = xyz - self.COM
class Cell(object):
    """
    Handles all operations related to developing and maintaining 
    periodic boundaries
    """

    def __init__(self):
        self.cell = np.zeros((3,3))
        self.origins = np.zeros((3,3)) 
        # cell parameters
        self.acell = np.zeros(6)
        # cell with unit vectors
        self.ncell = np.zeros((3,3))
        # inverted cell
        self.icell = np.zeros((3,3))
        # index to keep track of number of vectors added
        self.index = -1

    def valid_vector(self, vect):
        """
        checks a vector against existing cell vectors
        """
        for i in range(self.index+1):
            if self.line_test(i, vect):
                print "parallel boundary vector exists!"
                return False
        if self.index == 1:
            if self.plane_test(vect):
                print "vector lies in the same plane!"
                return False
        if self.index == 2:
            # This may change in the future, if changeable 
            # periodic boundaries are considered.
            return False
        return True

    def add_vector(self, vect):
        self.index +=1
        self.cell[self.index] = vect
        self.ncell[self.index] = normalize(vect)

    def plane_test(self, vector, tol=None):
        """
        Check to see if the vector is planar with the other
        cell vectors.
        """
        if tol is not None:
            tol = tol
        else:
            tol = 1e-2

        nvect = normalize(vector)
        for i in range(self.index+1):
            for j in range(i+1, self.index+1):
                test = dot(nvect, cross(self.ncell[i],self.ncell[j]))
                if np.allclose(test, 0., atol=tol):
                    return True
        return False

    def line_test(self, icell, vector, tol=None):
        """
        Check to see if the vector is on the same line as other
        cell vectors.
        """
        if tol is not None:
            tol = tol
        else:
            tol = 1e-2
        nvect = normalize(vector)
        ncell = normalize(self.cell[icell])
        xprod = cross(nvect,ncell)
        return np.allclose(xprod, np.zeros(3), atol=tol)

    def invert_cell(self):
        """ get the inverted cell for pbc corrections"""
        det = determinant(self.cell)

        a = self.cell[1][1] * self.cell[2][2] - \
            self.cell[1][2] * self.cell[2][1]
        b = self.cell[1][2] * self.cell[2][0] - \
            self.cell[1][0] * self.cell[2][2]
        c = self.cell[1][0] * self.cell[2][1] - \
            self.cell[1][1] * self.cell[2][0]
        d = self.cell[0][2] * self.cell[2][1] - \
            self.cell[0][1] * self.cell[2][2]
        e = self.cell[0][0] * self.cell[2][2] - \
            self.cell[0][2] * self.cell[2][0]
        f = self.cell[2][0] * self.cell[0][1] - \
            self.cell[0][0] * self.cell[2][1]
        g = self.cell[0][1] * self.cell[1][2] - \
            self.cell[0][2] * self.cell[1][1]
        h = self.cell[0][2] * self.cell[1][0] - \
            self.cell[0][0] * self.cell[1][2]
        i = self.cell[0][0] * self.cell[1][1] - \
            self.cell[0][1] * self.cell[1][0]

        self.icell = array([[a,b,c],[d,e,f],[g,h,i]]) / det

        #TODO(pboyd): this is the transpose of the np.linalg.inv() routine
        # figure out which one is correct for the cell vectors

    def cellparams(self):
        """ get cell parameters based on the cell vectors"""

        self.acell[:3] = [length(x) for x in self.cell]
        self.acell[3] = calc_angle(self.cell[1], self.cell[2])*RAD2DEG
        self.acell[4] = calc_angle(self.cell[0], self.cell[2])*RAD2DEG
        self.acell[5] = calc_angle(self.cell[0], self.cell[1])*RAD2DEG

    def complete_box(self):
        """Test to see if the full set of vectors are in place"""
        tol = 1.e-3
        volume = np.dot(np.cross(self.cell[0],self.cell[1]),
                        self.cell[2])
        if np.allclose(volume, 0., atol=tol):
            return False
        else:
            self.cellparams()
            self.invert_cell()
            return True


class Structure(object):
    """
    Structure contains all the machinery necessary to build a MOF from
    simple organic and metallic building units

    """
 
    def __init__(self, sbu_array):
        self.pbc = Cell() 
        # the list mof keeps track of all the SBUs
        self.mof = []
        # keep track of the most recent string
        self.string = ""
        # list to keep track of bonding between SBU's
        self.connectivity = []
        # maxsbu is the max number of sbus before termination..
        self.maxsbu = 99
        # list of SBU's
        self.sbu_array = sbu_array
        for i, j in enumerate(self.sbu_array):
            j.index = i

    def reset(self):
        sbuarray = self.sbu_array
        self.pbc = None
        self.__init__(sbuarray)

    def try_iterate_string(self, string):
        """check whether a string will form an appropriate bond"""
        values = string.split("-")
        xarray = [int(i) for i in values]
        isbu = xarray[0]
        ibond = xarray[1]
        nsbu = xarray[2]
        nbond = xarray[3]
        nangle = xarray[4]
        xarray.reverse()

        maxisbu = self.maxsbu # isbu will be continually growing.
        if isbu >= len(self.mof):
            return False
        maxibond = len(self.mof[isbu].connectpoints)
        maxnsbu = len(self.sbu_array)
        maxnbond = len(self.sbu_array[nsbu].connectpoints)
        if self.mof[isbu].ismetal:
            maxnangle = len(self.sbu_array[nsbu].connectangles[nbond])
        else:
            maxnangle = len(self.mof[isbu].connectangles[ibond])
        # reversed from string order for iteration purposes
        maxarray = [maxnangle, maxnbond, maxnsbu, maxibond, maxisbu]
        # assume only one type of bond (correct for both BTC and Cu-paddlewheel)
        maxarray = [maxnangle, 1, maxnsbu, maxibond, 99]

        # iterate the string
        for ival, value in enumerate(xarray):
            if value == (maxarray[ival]-1):
                xarray[ival] = 0
            else: 
                xarray[ival] += 1
                self.string =  "%i-%i-%i-%i-%i"%(tuple(xarray[::-1]))
                break
        if isbu == maxisbu:
            # TODO(pboyd): Raise flag for termination.
            return False
        xarray.reverse()
        if xarray[0] >= len(self.mof):
            return False
        # TODO(pboyd): further checks?
        # not a metal to organic bond?
        falsetest2 = self.sbu_array[xarray[2]].ismetal == \
                        self.mof[xarray[0]].ismetal
        # existing SBU is already bonded?
        falsetest3 = self.connectivity[xarray[0]][xarray[1]] is not None

        if (falsetest2) or (falsetest3):
            return False
        else:
            return True 

    def apply_string(self, string):
        """applies the directions from a string"""

        self.string = string
        strlist = string.split("-")
        sbu1 = int(strlist[0])
        bond1 = int(strlist[1])
        sbutype2 = int(strlist[2])
        bond2 = int(strlist[3])
        # if sbu1 is a metal, the angle should be from the organic
        # index (sbu2)
        iangle = int(strlist[4])

        # add sbu1 to sbu2
        self.add_new_sbu(sbu1, bond1, sbutype2, bond2, iangle)

    def bad_addition(self, string, struct):
        """apply checks to see if the new addition is valid"""
        strlist = string.split("-")

        # assumes sbu2 is the newest addition to the MOF
        sbu1 = int(strlist[0])
        bond1 = int(strlist[1])
        sbu2 = len(self.mof) - 1
        bond2 = int(strlist[3])
        # FIXME(pboyd): overlap not returning true in some cases.
        if self.overlap(sbu2):
            return True
        return False

    def sbu_check(self, sbu_ind1, tol=None):
        """
        Checks an SBU for 
        1. bonding locally,
        2. bonding by existing pbc
        """
        # FIXME(pboyd): instead of del(bond[i]) change to removing the 
        # value from the list of sbus     {1:[[1,2],[2,3],[4,5]]}
        # so instead of removing the key---^         ^
        # remove the list entry ---------------------|

        # tolerance for bond distance
        if tol is not None:
            tol = tol
        else:
            tol = 5e-1
        # return dictionary of anti-parallel, unsaturated bonds within MOF
        # note, may be a problem if bonds are close but not exactly
        # anti-parallel
        bonds = self.antiparallel_bonds(sbu_ind1)
        # check to see if any bonds are formed
        if not bonds:
            info("No new bonds formed with SBU %i"%(sbu_ind1))
            return

        # TODO(pboyd):check for alignment (angle_vect + angles)??

        # generate vectors for all the eligible bonds 
        bondvects = self.gen_bond_vectors(sbu_ind1, bonds)

        bondlist = [(i,tuple(j)) for i in bonds for j in bonds[i]]
        # !!!problems if connection site forms multiple possible bonds!!!
        for bond in bondlist:
            if length(bondvects[bond]) <= tol:
                # bonds are close
                info("Local bond found between "+
                        "SBU %i, bond %i, "%(sbu_ind1, bond[0])+
                        "and SBU %i, bond %i."%(bond[1][0],bond[1][1]))
                self.join_sbus(sbu_ind1, bond[0], bond[1][0], 
                               bond[1][1], False)
                try:
                    del bonds[bond[0]]
                except:
                    info("Tried to bond SBU %i bond %i to SBU %i "
                         %(sbu_ind1, bond[0], bond[1][0]) +
                         "bond %i, but it has already been bonded."
                         %(bond[1][1]))

        bondlist = [(i,tuple(j)) for i in bonds for j in bonds[i]]
        # check for possible periodic boundary
        for bond in bondlist:
            if self.pbc.valid_vector(bondvects[bond]):
                info("New periodic boundary formed between "+
                    "SBU %i, bond %i, "%(sbu_ind1, bond[0])+
                    "and SBU %i, bond %i "%(bond[1][0],bond[1][1])+
                    " with the vector (%f, %f, %f)."%(tuple(bondvects[bond])))
                self.join_sbus(sbu_ind1, bond[0], bond[1][0], 
                              bond[1][1], True)
                self.pbc.add_vector(bondvects[bond])
                self.pbc.origins[self.pbc.index] = \
                        self.mof[sbu_ind1].connectpoints[bond[0]][:]
                try:
                    del bonds[bond[0]]
                except:
                    info("Tried to bond SBU %i bond %i to SBU %i "
                         %(sbu_ind1, bond[0], bond[1][0]) +
                         "bond %i, but it has already been bonded."
                         %(bond[1][1]))

        bondlist = [(i,tuple(j)) for i in bonds for j in bonds[i]]

        # check for "local" bonds via periodic boundaries
        for bond in bondlist:
            pbc_shifted = self.apply_pbc(bondvects[bond])
            if length(pbc_shifted) <= tol:
                info("Local bond found between " +
                    "SBU %i, bond %i, "%(sbu_ind1, bond[0])+
                    "and SBU %i, bond %i."%(bond[1][0],bond[1][1])+
                    " Found with the application of pbc's.")
                self.join_sbus(sbu_ind1, bond[0], bond[1][0], 
                              bond[1][1], True)
                try:
                    del bonds[bond[0]]
                except:
                    info("Tried to bond SBU %i bond %i to SBU %i "
                         %(sbu_ind1, bond[0], bond[1][0]) +
                         "bond %i, but it has already been bonded."
                         %(bond[1][1]))

        # announce connections which did not form bonds
        bondlist = [(i,tuple(j)) for i in bonds for j in bonds[i]]
        if bondlist is not None:
            for bond in bondlist:
                pbc_shifted = self.apply_pbc(bondvects[bond])
                if self.remote_attach_sbus(pbc_shifted):
                    info("Remote bond found between " +
                        "SBU %i, bond %i, "%(sbu_ind1, bond[0])+
                        "and SBU %i, bond %i."%(bond[1][0],bond[1][1])+
                        " Found between existing pbc's.")
                    self.join_sbus(sbu_ind1, bond[0], bond[1][0], 
                              bond[1][1], True)
                else:
                    info("Nothing was done between "+
                        "SBU %i, bond %i, "%(sbu_ind1, bond[0]) +
                        "and SBU %i, bond %i. "%(bond[1][0], bond[1][1])+
                        "Even though they had parallel alignment.")
                    info("(%f, %f, %f)"%(tuple(bondvects[bond])))

    def remove_like_sbus(self, sbu_ind1, bonds):
        """
        Eliminate SBUs which are the same type.  This may be changed
        in the future to allow for metal-metal bonds?
        """
        copybonds = bonds.copy()
        for i in copybonds.iteritems():
            for ind, sbu in enumerate(i[1]):
                if self.mof[sbu_ind1].ismetal == self.mof[sbu[0]].ismetal:
                    bonds[i[0]].pop(ind)

            if not i[1]:
                del bonds[i[0]]
        return bonds

    def remote_attach_sbus(self, bondvector):
        """
        Join two SBUs if they have anti-parallel bonding vectors,
        and existing pbc vectors can join them.
        """
        # tolerance for matching vector lengths.
        tol = 0.5
        # check if a local bond can be applied after applying pbc
        # conditions.  NOTE: this is still considered a "remote"
        # attachment since pbc conditions are used.

        if self.pbc.index >= 0:
            proj_a = project(bondvector, self.pbc.cell[0])
            test_a = np.allclose(length(self.pbc.cell[0]), 
                                 length(proj_a), atol=tol)
            null_a = np.allclose(proj_a, np.zeros(3), atol=tol)
            # case 1: the bond between two sbus can already be
            # established by a single boundary condition
            for i in range(self.pbc.index+1):
                if self.pbc.line_test(i, bondvector):
                    if np.allclose(length(self.pbc.cell[i]), 
                                   length(bondvector), atol=tol):
                        return True

        if self.pbc.index >= 1:
            proj_b = project(bondvector, self.pbc.cell[1])
            test_b = np.allclose(length(self.pbc.cell[1]), 
                                 length(proj_b), atol=tol)
            null_b = np.allclose(proj_b, np.zeros(3), atol=tol)
            # case 2: the bond between two sbus is at the corner of
            # two periodic boundaries.
            if self.pbc.index == 1:
                if test_a and test_b:
                    return True

        if self.pbc.index == 2:
            proj_c = project(bondvector, self.pbc.cell[2])
            test_c = np.allclose(length(self.pbc.cell[2]), 
                                 length(proj_c), atol=tol)
            null_c = np.allclose(proj_c, np.zeros(3), atol=tol)
            # case 2: the bond between two sbus is at the corner of
            # two periodic boundaries.
            if test_a and test_b and null_c:
                return True
            if test_a and test_c and null_b:
                return True
            if test_c and test_b and null_a:
                return True

            # case 3: the bond between two sbus is at the corner of 
            # three periodic boundaries.
            if test_a and test_b and test_c:
                return True

        return False

    def gen_bond_vectors(self, sbu1, bond_dictionary):
        """
        returns a list of bond vectors 
        """
        sorted_list = [x for x in bond_dictionary.iteritems()]
        sorted_list.sort(key = lambda x: x[0])

        vectors = {}
        for bond1 in sorted_list:
            for sbu2 in bond1[1]:
                connecting_point1 = self.mof[sbu1].connectpoints[bond1[0]]
                connecting_point2 = self.mof[sbu2[0]].connectpoints[sbu2[1]]
                vect = connecting_point2 - connecting_point1
                # key is (bond, [sbu, bond])
                vectors[(bond1[0],tuple(sbu2))] = vect
        return vectors

    def antiparallel_bonds(self, sbu1):
        """
        Returns a dictionary of parallel bonds with each unsaturated
        bond of sbu1.  The format is {bond :[sbu#, bond#]}
        """
        # exclude bonds in sbu1 which are already bonded
        bonds = [i for i in range(len(self.connectivity[sbu1])) 
                 if self.connectivity[sbu1][i] is None]
        # remove sbu1 from scan
        mofscan = [i for i in range(len(self.mof)) if i is not sbu1]
        bonding_dic = {}
        
        for sbu2 in mofscan:
            if self.mof[sbu1].type != self.mof[sbu2].type:
                bonds2 = [i for i in range(len(self.connectivity[sbu2]))
                      if self.connectivity[sbu2][i] is None]
                for bond2 in bonds2:
                    for bond1 in bonds:
                        vect1 = self.mof[sbu1].connect_vector[bond1]
                        vect2 = self.mof[sbu2].connect_vector[bond2]
                        # note the tolerance here will dictate how many
                        # structures will be found
                        if antiparallel_test(vect1, vect2, 1e-1):
                            bonding_dic.setdefault(bond1,[]).append([sbu2,bond2])
        return bonding_dic

    def overlap(self, sbu1, tol=None):
        """
        checks if one of the SBUs in structure (istruct) is
        overlapping with the other SBUs
        """
        # distance tolerance in angstroms
        if tol is not None:
            tol = tol
        else:
            tol = .5

        mofrange = [i for i in range(len(self.mof)) if i != sbu1]
        for sbu2 in mofrange:
            sbu2 = self.mof[sbu2]
            for atom2 in sbu2.coordinates:
                for atom1 in self.mof[sbu1].coordinates:
                    dist = self.apply_pbc(atom1 - atom2)
                    if length(dist) <= tol:
                        return True
        return False

    def coordinate_dump(self):
        """ 
        Dumps all the atom labels, xyz coordinates and what SBU index
        they belong to
        """
        coords, atoms, sbuind = [], [], []
        for sbu in range(len(self.mof)):
            coords.append([i for i in self.mof[sbu].coordinates])
            atoms.append([i for i in self.mof[sbu].atomlabel])
            sbuind.append([sbu for i in range(len(self.mof[sbu].coordinates))])

        coords = array(list(flatten(coords)))
        atoms = list(flatten(atoms))
        sbuind = list(flatten(sbuind))
        return (atoms, coords, sbuind)

    def add_new_sbu(self, sbu1, bond1, sbutype2, bond2, iangle):
        """adds a new sbutype2 to the growing MOF"""

        #self.xyz_debug()
        #dump = self.coordinate_dump()
        #write_xyz("history", dump[0], dump[1], self.pbc.cell, self.pbc.origins)

        # use angle listed on organic species bond type
        if self.sbu_array[sbutype2].ismetal:
            angle = self.mof[sbu1].connectangles[bond1][iangle]
        else:
            angle = self.sbu_array[sbutype2].connectangles[bond2][iangle]

        self.mof.append(copy.deepcopy(self.sbu_array[sbutype2]))
        sbu2 = len(self.mof) - 1
        self.connectivity.append([None]*len(self.mof[sbu2].connectpoints))
        info(
            "Added %s, SBU %i, bond %i to SBU %i, %s bond %i."
            %(self.mof[sbu2].type, sbu2, bond2, sbu1,
              self.mof[sbu1].type, bond1))
        #self.xyz_debug()
        #dump = self.coordinate_dump()
        #write_xyz("history", dump[0], dump[1], self.pbc.cell, self.pbc.origins)

        # align sbu's by Z vector
        self.sbu_align(sbu1, bond1, sbu2, bond2)

        self.xyz_debug()
        dump = self.coordinate_dump()
        write_xyz("history", dump[0], dump[1], self.pbc.cell, self.pbc.origins)

        # rotate by Y vector
        self.bond_align(sbu1, bond1, sbu2, bond2, angle) 
    
    def bond_align(self, sbu1, bond1, sbu2, bond2, rotangle):
        """
        Align two sbus by their orientation vectors
        """

        sbu1 = self.mof[sbu1]
        sbu2 = self.mof[sbu2]

        axis = normalize(sbu2.connect_vector[bond2])
        angle = calc_angle(sbu2.anglevect[bond2],
                           sbu1.anglevect[bond1])

        # origin of rotation
        rotpoint = sbu2.connectpoints[bond2] 

        # test to see if the rotation is in the right direction
        R = rotation_matrix(axis, angle)
        test_vect = array(sbu2.anglevect[bond2][:] * R)[0]
        test_angle = calc_angle(sbu1.anglevect[bond1],
                           test_vect)
        # FIXME(pboyd): the tolerance below is quite large,
        if not np.allclose(test_angle, 0., atol=5e-2):
            axis = -1.*axis

        # rotation of SBU
        self.sbu_rotation(sbu2, rotpoint, axis, angle + rotangle)
        angle2 = calc_angle(sbu1.anglevect[bond1],
                           sbu2.anglevect[bond2])
        if not np.allclose(angle2, 0., atol = 5e-2):
            print "PROBLEM", angle2, "rotational axis", axis, "rotational angle", angle, "test vector", test_vect, "test angle", test_angle, "test vector2", sbu1.anglevect[bond1]

    def sbu_align(self, sbu1, bond1, sbu2, bond2):
        """
        Align two sbus, were sbu1 is held fixed and
        sbu2 is rotated and shifted to bond with sbu1
        """

        # the order of the cross product is important to 
        # signify which vector remains stationary (statvect)
        # and which will rotate.  This axis will determine the
        # proper rotation when applying the rotation matrix

        sbu1 = self.mof[sbu1]
        sbu2 = self.mof[sbu2]

        axis = rotation_axis(sbu2.connect_vector[bond2], 
                            -1.*sbu1.connect_vector[bond1])
        angle = calc_angle(sbu2.connect_vector[bond2], 
                           -1.*sbu1.connect_vector[bond1])
        # origin of rotation
        rotpoint = sbu2.centre_of_mass()

        while (np.allclose(axis, np.zeros(3), atol=1e-4)):
            randaxis = normalize(array([random(), random(), random()]))
            randangle = uniform(0., np.pi/3.)
            self.sbu_rotation(sbu2, rotpoint, randaxis, randangle)
            axis = rotation_axis(sbu2.connect_vector[bond2], 
                             -1.*sbu1.connect_vector[bond1])
            angle = calc_angle(sbu2.connect_vector[bond2], 
                           -1.*sbu1.connect_vector[bond1])

        # rotate sbu2
        self.sbu_rotation(sbu2, rotpoint, axis, angle)

        # shift the coordinates by the shift vectors 
        shiftvector = (sbu1.connectpoints[bond1] - 
                       sbu2.connectpoints[bond2])
        # bond points are then shifted by the shiftvector
        sbu2.connectpoints += shiftvector
        sbu2.coordinates += shiftvector

    def sbu_rotation(self, sbu, rotpoint, axis, angle):
        """
        Takes a growing structure and rotates the coordinates
        of sbu about a selected bond.  
        This only applies to linkers with a single bond to the 
        SBU - bidentate bonds must have other considerations
        """

        R = rotation_matrix(axis, angle)
        # R rotates these points about the origin, hence a correction
        # must be applied if the position of rotation is not the 
        # origin, using rotpoint * (I - R) 
        C = array(rotpoint * (np.identity(3) - R))
        
        sbu.coordinates = array(sbu.coordinates * R) + C
        sbu.connectpoints = array(sbu.connectpoints * R) + C
        # connect_vector and anglevect are ALWAYS centered at the origin!!
        sbu.connect_vector = array(sbu.connect_vector * R)
        sbu.anglevect = array(sbu.anglevect * R)
    
    def final_coords(self):
        """
        convert all coordinates to the shifted coordinates derived
        from the periodic boundary conditions
        """
        for mof in self.mof:
            for ibond, bond in enumerate(mof.connectpoints):
                mof.connectpoints[ibond] = self.apply_pbc(bond)
            for icoord, coord in enumerate(mof.coordinates):
                mof.coordinates[icoord] = self.apply_pbc(coord)
        return
                
    def mof_reorient(self):
        """
        Re-orients the cell vectors and the coordinates such that the
        first cell vector points in the x cartesian axis and the 
        second cell vector points in the xy cartesian plane
        """
        xaxis = array([1.,0.,0.])
        yaxis = array([0.,1.,0.])
        # first rotation to the x-axis
        x_rotangle = calc_angle(self.pbc.cell[0], xaxis)
        x_rotaxis = rotation_axis(self.pbc.cell[0], xaxis)
        RX = rotation_matrix(x_rotaxis, x_rotangle)
        self.pbc.cell = array(self.pbc.cell * RX)
        for sbu in self.mof:
            sbu.coordinates = array(sbu.coordinates * RX)
            sbu.connect_vector = array(sbu.connect_vector * RX)
            sbu.anglevect = array(sbu.anglevect * RX)
            sbu.connectpoints = array(sbu.connectpoints * RX)

        # second rotation to the xy - plane
        projx_b = self.pbc.cell[1] - project(self.pbc.cell[1], xaxis)
        self.icell = np.zeros((3,3))
        xy_rotangle = calc_angle(projx_b, yaxis)
        xy_rotaxis = xaxis
        RXY = rotation_matrix(xy_rotaxis, xy_rotangle)
        # test to see if the rotation is in the right direction
        testvect = array(projx_b * RXY)[0]
        testangle = calc_angle(testvect, yaxis)
        if not np.allclose(testangle, 0., atol = 1e-3):
            RXY = rotation_matrix(-1.*xy_rotaxis, xy_rotangle)

        self.pbc.cell = array(self.pbc.cell * RXY)
        for sbu in self.mof:
            sbu.coordinates = array(sbu.coordinates * RXY)
            sbu.connect_vector = array(sbu.connect_vector * RXY)
            sbu.anglevect = array(sbu.anglevect * RXY)
            sbu.connectpoints = array(sbu.connectpoints * RXY)

        # third change: -z to +z direction
        if self.pbc.cell[2][2] < 0.:
            RZ = rotation_matrix(array([1.,0.,0]), np.pi) 
            self.pbc.cell = array(self.pbc.cell * RZ)
            for sbu in self.mof:
                sbu.coorinates = array(sbu.coordinates * RZ)
                sbu.connect_vector = array(sbu.connect_vector * RZ)
                sbu.anglevect = array(sbu.anglevect * RZ)
                sbu.connectpoints = array(sbu.connectpoints * RZ)

        return

    def saturated(self):
        """
        returns False if at least one bond is unsaturated
        """
        for imof, mof in enumerate(self.connectivity):
            if [no_bond for no_bond in mof if no_bond is None]:
                return False
        return True

    def join_sbus(self, sbu1, bond1, sbu2, bond2, islocal):
        """ Joins two sbus together"""
        if islocal:
            self.connectivity[sbu1][bond1] = sbu2
            self.connectivity[sbu2][bond2] = sbu1
        else:
            self.connectivity[sbu1][bond1] = -1*sbu2
            self.connectivity[sbu2][bond2] = -1*sbu1
        

    def apply_pbc(self, vector):
        """
        applies periodic boundary conditions to a vector
        """

        # FIXME(pboyd): the pbcs when the index = 1 appear broken
        A = self.pbc.cell[0]
        B = self.pbc.cell[1]
        C = self.pbc.cell[2]
        if np.allclose(vector, np.zeros(3)):
            return vector
        if self.pbc.index == 0:
            vector = vectorshift(vector, A)
            #vector = vector
        elif self.pbc.index == 1:
            proj_a = project(vector, self.pbc.cell[0])
            proj_b = project(vector, self.pbc.cell[1])
            a = length(proj_a) / length(self.pbc.cell[0])
            b = length(proj_b) / length(self.pbc.cell[1])
            if antiparallel_test(proj_a, self.pbc.cell[0]):
                a = -1.*a
            if antiparallel_test(proj_b, self.pbc.cell[1]):
                b = -1.*b
            # round a and b to 4 decimal places.  Precision can cause
            # problems if a and b are really close to 1, but under
            # so math.floor doesn't subtract by 1.  Same goes with c
            # below.
            a = round(a, 4)
            b = round(b, 4)
            vector = vector - (math.floor(a)*A) - (math.floor(b)*B)

        elif self.pbc.index == 2:
            bca = dot(B, cross(C,A))
            a = dot(C, cross(vector, B)) / bca
            b = dot(C, cross(vector, A)) / (-1.*bca)
            c = dot(B, cross(vector, A)) / bca
            a = a - math.floor(round(a,4))
            b = b - math.floor(round(b,4))
            c = c - math.floor(round(c,4))
            vector = (a)*A + (b)*B + (c)*C

           #vector = vector
        return vector

    def xyz_debug(self, filename=None):
        """
        Writes all the bonding and atom info to an xyz file
        """
        if filename is None:
            filename = "debug.xyz"
        else:
            filename = filename + ".xyz"
        xyzfile = open(filename, "a")

        bondformat = "%s%12.5f%12.5f%12.5f " +\
                     "atom_vector%12.5f%12.5f%12.5f " +\
                     "atom_vector%12.5f%12.5f%12.5f\n"
        #cellformat = "%s%12.5f%12.5f%12.5f " +\
        #             "atom_vector%12.5f%12.5f%12.5f\n" 
        cellformat = "%s%12.5f%12.5f%12.5f " +\
                     "atom_vector%12.5f%12.5f%12.5f " +\
                     "atom_vector%12.5f%12.5f%12.5f " +\
                     "atom_vector%12.5f%12.5f%12.5f\n" 
        atomformat = "%s%12.5f%12.5f%12.5f\n"
        
        line = []
        atomcount = 1
        line.append(cellformat%(tuple(["C"] + [0.,0.,0.] +
                    list(self.pbc.cell[0]) + list(self.pbc.cell[1])+ list(self.pbc.cell[2]))))
        #atomcount = 3
        #for icell, cell in enumerate(self.pbc.cell):
        #    line.append(cellformat%(tuple(["C"] + list(self.pbc.origins[icell])
        #                + list(cell))))
        for sbu in self.mof:
            for ibond, bondpt in enumerate(sbu.connectpoints):
                atomcount += 1
                bondpt = self.apply_pbc(bondpt)
                line.append(bondformat%(tuple(["F"]+list(bondpt)+
                                        list(sbu.connect_vector[ibond])+
                                        list(sbu.anglevect[ibond]))))
            for icoord, coord in enumerate(sbu.coordinates): 
                atomcount += 1
                coord = self.apply_pbc(coord)
                line.append(atomformat%(tuple([sbu.atomlabel[icoord]]+list(coord))))

        xyzfile.write("%i\nDebug\n"%atomcount)
        xyzfile.writelines(line)
        xyzfile.close()
    
def vectorshift(shift_vector, root_vector):
    """Shift a vector to within the bounds of a root vector"""

    if(np.allclose(np.zeros(3), root_vector)):
        return shift_vector
    else:
        # break function if the shift_vector is zero
        if np.allclose(shift_vector, np.zeros(3)):
            return shift_vector
        # project the vector to be shifted onto the root vector
        proj_vect = project(shift_vector, root_vector)

        # express the projected vector as a fraction of the root vector
        fraction = length(proj_vect) / length(root_vector)
        if antiparallel_test(proj_vect, root_vector):
            fraction = -1. * fraction
        # return the shift_vector which has been shifted such that
        # it's projection lies within the bounds of the root_vector
        return (shift_vector - math.floor(fraction)*root_vector)

def determinant(matrix):
    """ calculates the determinant of a 3x3 matrix"""
    return np.dot(matrix[0],np.cross(matrix[1],matrix[2]))


def coord_matrix(atom, coord):
    """
    generate a coordination matrix for a list of x,y,z coordinates
    """
    connect = np.matrix(np.zeros((len(coord),len(coord))))
    bonding = []
    numbonds = []
    for i in range(len(coord)):
        bonding.append([])
        numbonds.append(0)
        for j in range(len(coord)):
            dist = length(coord[i],coord[j])
            tol = bond_tolerance(atom[i], atom[j])
            if (i != j) and (dist <= tol):
                connect[i,j] = dist
                connect[j,i] = dist
                numbonds[i] += 1
                bonding[i].append(j)
            else:
                connect[i,j] = 0.
                connect[j,i] = 0.
    return (connect, bonding, numbonds)

def bond_tolerance(atom1, atom2):
    """
    returns a tolerance value for atom - atom distances
    case specific
    """
    # This function will need some elaboration
    # TODO(pboyd): these conditions are not robust and
    # should be changed when the code becomes bigger
    if("O" in [atom1,atom2])and("C" in [atom1,atom2]):
        return 1.6
    elif("X" in (atom1,atom2))and(("Y" not in (atom1,atom2))and("Z" not in (atom1,atom2))):
        return 0.
    elif("X" in (atom1,atom2))and(("Y" in (atom1,atom2))or("Z" in (atom1,atom2))):
        return 1.1
    elif("Y" in (atom1,atom2))and("X" not in (atom1,atom2)):
        return 0.
    elif("Y" in (atom1,atom2))and("X" in (atom1,atom2)):
        return 1.1
    elif("Z" in (atom1,atom2))and("X" not in (atom1,atom2)):
        return 0.
    elif("Z" in (atom1,atom2))and("X" in (atom1,atom2)):
        return 1.1
    else:
        return 1.6

def length(coord1, coord2=None):
    """ 
    Returns the length between two vectors.
    If only one vector is specified, it returns
    the length of that vector from the origin
    """
    if coord2 is not None:
        coord2 = coord2
    else:
        coord2 = np.zeros(3)
    if sum(coord1) == 0. and sum(coord2) == 0.:
        return 0.
    vector = coord2 - coord1
    return math.sqrt(dot(vector, vector))

def rotation_axis(vect1, vect2):
    """
    axis of rotation for two vectors (90 degrees from both)
    """
    vect1 = normalize(vect1)
    vect2 = normalize(vect2)
    return normalize(np.cross(vect1, vect2))

def calc_angle(vect1, vect2):
    """ determines angle between vector1 and vector2"""
    dot12 = dot(vect1, vect2) 
    dist11 = length(vect1)
    dist22 = length(vect2)
    # clamps the angle coefficient to a min or max of -1, 1 so no 
    # error is returned when calculating the acos.
    angle_coefficient =  min(max(dot12/(dist11*dist22), -1.0),1.0)
    return math.acos(angle_coefficient)

def rotation_matrix(axis, angle):
    """
    returns a (3,3) rotation matrix based on
    axis and angle.
    The rotation is counter-clockwise when looking
    down the axis so choice of sign is important here
    """
    ux = axis[0]
    uy = axis[1]
    uz = axis[2]

    matrix = np.matrix([[
        np.cos(angle) + ux * ux * (1 - np.cos(angle)),
        uy * ux * (1 - np.cos(angle)) + uz * np.sin(angle),
        uz * ux * (1 - np.cos(angle)) - uy * np.sin(angle)],
        [
        ux * uy * (1 - np.cos(angle)) - uz * np.sin(angle),
        np.cos(angle) + uy * uy * (1 - np.cos(angle)),
        uz * uy * (1 - np.cos(angle)) + ux * np.sin(angle)],
        [
        ux * uz * (1 - np.cos(angle)) + uy * np.sin(angle),
        uy * uz * (1 - np.cos(angle)) - ux * np.sin(angle),
        np.cos(angle) + uz * uz * (1 - np.cos(angle))]])
    return matrix

def planar_test(coordinates, tol=None):
    """
    Tests whether four or more points are within a plane.
    Returns False if only one of the points is co-planar
    """

    # FIXME(pboyd): this appears to be broken
    if tol is None:
        tol = 1e-2
    else:
        tol = tol
    # Define a plane by three points
    planevector1 = normalize(coordinates[1] - coordinates[0])

    index = 1
    while np.allclose(sum(planevector1), 0., atol=tol):
        index += 1
        planevector1 = normalize(coordinates[index] - coordinates[0])
    # test points for colinearity with the first vector
    for coord in coordinates[index+1:]:
        if(not sum(coord) == 0.):
            if(not linear_test(array([coordinates[0],coordinates[1],coord]))):
                planevector2 = normalize(coord - coordinates[0])
                break

    # test to see if planevector2 is allocated (if there is a plane
    # within the coordinates) otherwise co-linear
    try:
        planevector2
    except:
        return False

    for coord in coordinates[index+1:]:
       # If not colinear test for planarity
        if sum(coord) != 0.:
            newvect = normalize(coord - coordinates[0])

            # test = 0 if co-planar
            test = dot(newvect, cross(planevector1, planevector2))
            if (not np.allclose(0., test, atol=tol)):
                return False
    
    return True
        

def normalize(vector):
    """changes vector length to unity"""
    if sum(vector) == 0.:
        return vector
    return vector / math.sqrt(dot(vector, vector))

def linear_test(coordinates, tol=None):
    """
    Tests for three or more points lying within the same
    line.
    """
    if tol is None:
        tol = 1e-2
    else:
        tol = tol
 
    # test by cross-product.  If linear, should be the zero vector
    vector1 = normalize(coordinates[1] - coordinates[0])

    for point in coordinates[2:]:
        if sum(point) != 0.:
            vector2 = normalize(point - coordinates[0])
            crossprod = np.cross(vector2,vector1)
            if np.allclose(crossprod, np.zeros(3), atol=tol):
                return True

    return False

def points_close(point1, point2, tol=None):
    """
    Checks for the proximity of two points in 3d space with a given
    tolerance
    """
    if tol is not None:
        tol = tol
    else:
        tol = 1e-1

    dist = length(point1, point2)
    if dist <= tol:
        return True
    return False


def project(vector1, vector2):
    """ projects vector1 onto vector2, returns that vector"""
    angle = calc_angle(vector1, vector2)
    unit2 = normalize(vector2)
    return length(vector1) * np.cos(angle) * unit2

def degenerate_test(vector1, vector2, tol=None):
    """
    test for two vectors lying on the same line, parallel or 
    anti-parallel.
    """
    if tol is not None:
        tol = tol
    else:
        tol = 5e-2
    vector1 = normalize(vector1)
    vector2 = normalize(vector2)
    test1 = np.allclose(cross(vector1, vector2), 
                        np.zeros(3), atol=tol)
    test2 = np.allclose(abs(dot(vector1, vector2)), 1., atol=tol)
    if test1 and test2:
        return True
    else:
        return False

def antiparallel_test(vector1, vector2, tol=None):
    """test for two vectors being in anti-parallel orientations"""

    if tol is not None:
        tol = tol
    else:
        tol = 5e-2
    vector1 = normalize(vector1)
    vector2 = normalize(vector2)
    test1 = np.allclose(cross(vector1, vector2), 
                        np.zeros(3), atol=tol)
    test2 = np.allclose(dot(vector1, vector2), -1., atol=tol)
    if test1 and test2:
        return True
    else:
        return False

def parallel_test(vector1, vector2, tol=None):
    """test for two vectors being in parallel orientations"""
    if tol is not None:
        tol = tol
    else:
        tol = 5e-2 
    vector1 = normalize(vector1)
    vector2 = normalize(vector2)
    test1 = np.allclose(cross(vector1, vector2), 
                        np.zeros(3), atol=tol)
    test2 = np.allclose(dot(vector1, vector2), 1., atol=tol)
    if test1 and test2:
        return True
    else:
        return False

def write_pdb(label, atoms, coords, acell):
    pdbfile = open('%s.pdb' % label, 'w')
    today = date.today()
    lines = []
    atomformat1 = "%-6s%5i %-4s %3s %1s%4i%s   "
    atomformat2 = "%8.3f%8.3f%8.3f%6.2f%6.2f            %-3s%2s\n"
    lines.append("%6s   Created: %s \n%-6s" % ('REMARK', 
                  today.strftime("%A %d %B %Y"), 'CRYST1'))
    lines.append("%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P1\n"
                  % (tuple(acell)))
    for atom in range(len(atoms)):
        lines.append(atomformat1%('ATOM',atom+1,atoms[atom],"MOL","X",1," ") +
                     atomformat2%(tuple(list(coords[atom])+[1.,0.,atoms[atom],0])))
    lines.append("TER")
    pdbfile.writelines(lines)
    pdbfile.close()

def flatten(seq):
    """ flattens a nested list to a 1d list """
    for x in seq:
        if type(x) is list:
            for y in flatten(x):
                yield y
        else:
            yield x

def write_xyz(label, atoms, coords, cell, origin=None):
    if origin is not None:
        origin = origin
    else:
        origin = np.zeros((3,3))
    xyzfile = open('%s.xyz' % label, 'a')
    xyzfile.write('%i\ncoordinate dump\n'%(len(coords)+3))
    for j in range(len(cell)):
        xyzfile.write("%s%12.5f%12.5f%12.5f %s%12.5f%12.5f%12.5f\n"%(
            tuple(["C"] + list(origin[j]) +["atom_vector"] + list(cell[j]))))
    for i in range(len(coords)):
        xyzfile.write('%s%12.5f%12.5f%12.5f\n' % 
                (tuple([atoms[i]]+list(coords[i]))))

def dot(vector1, vector2):
    """
    returns the dot product of two three dimensional vectors
    """
    return sum(x*y for x,y in zip(vector1,vector2)) 

def cross(vector1, vector2):
    """
    returns the cross product of two three dimensional vectors
    """
    v0 = (vector1[1] * vector2[2]) - (vector1[2] * vector2[1])
    v1 = -1.*((vector1[0] * vector2[2]) - (vector1[2] * vector2[0]))
    v2 = (vector1[0] * vector2[1]) - (vector1[1] * vector2[0])
    return array([v0,v1,v2])


def main():
    """Default if run as an executable"""

    # initialize logging 
    Log()
    genstruct = Generate()

    genstruct.branched_generation()
if __name__ == '__main__':
    main()
