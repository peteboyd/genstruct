#!/usr/bin/env python

"""
GenStruct -- Generation of Structures for fapping.
"""

__version__ = "$Revision$"

import subprocess
import textwrap
import sys
import numpy as np
from numpy import array
from elements import WEIGHT, ATOMIC_NUMBER
import sample
import bookkeeping
from random import random, uniform, randrange, choice
from datetime import date

# Constants

RAD2DEG = 180./np.pi

class SBU(object):
    """
    """

    def __init__(self, type=None, atomlabel=None, coordinates=None,
                 pdbfile=None, ismetal=None):
        #self.log = bookkeeping.Log()
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
        if pdbfile is not None:
            self.from_pdb(pdbfile)
            if coordinates is not None:
                sys.exit(0)
                #self.log.error("conflicting coordinate reference")
            if atomlabel is not None:
                sys.exit(0)
                #self.log.error("conflicting atom reference")
        
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
        # keeping track of which connecting points are bonded
        self.connectivity = []
        # angle vectors used to align SBUs once bonded together
        self.anglevect = []

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

        pdbfile = open(filename)
        pdblines = pdbfile.readlines()
        atmindex = 0
        for line in pdblines:
            if line.lower().startswith('atom'):
                self.atomlabel.append(line[12:16].strip())
                self.coordinates.append(array([float(line[30:38]), \
                    float(line[38:46]), float(line[47:54])]))
        pdbfile.close()


    def check_unsat(self):
        """
        checking for unsaturated bonds
        """
        unsat = []
        # need to improve robustness of this code
        for i in range(len(self.nbonds)):
            if ((self.nbonds[i] <= 2)and(self.atomlabel[i] != "H")and
                (self.atomlabel[i] != "X")and(self.atomlabel[i] != "Y")):
                unsat.append(i)
        self.unsat = unsat

    def add_connect_vector(self):
        """
        checks for certain unsaturated atoms and adds connect_vector accordingly
        """

        purgeatoms = []
        for indx, atom in enumerate(self.atomlabel):
            # test for "X" atoms (considered a linking point)
            if atom == "X":
                # there should only be one bonding atom to an "X"
                bondatm=None
                for i in self.bonding[indx]:
                    if self.atomlabel[i] != "Y":
                        bondatm = i
                    if self.atomlabel[i] == "Y":
                        # angle vector required to align SBUs
                        self.anglevect.append(self.coordinates[i] - 
                                self.coordinates[indx])
                        purgeatoms.append(i)
                if bondatm is None:
                    self.log.error("X is not bound to an atom!!")
                vect = (self.coordinates[indx] - 
                        self.coordinates[bondatm]) 
                pt = self.coordinates[indx]
                self.connectpoints.append(pt)
                self.connect_vector.append(vect)
                self.connectivity.append(None)
                purgeatoms.append(indx)

            # TODO(pboyd): add other tests for different 
            # coordinating atoms

        # purge any coordinates in self.coordinates belonging to "X"'s
        # and "Y"'s as they are now part of self.connect_vector 
        # and self.connectpoints and self.anglevect
        self.purge(purgeatoms)

        if len(self.anglevect) == 0:
            sys.exit(0)
            #self.log.error("Error, no angle vectors read in input")
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

class Structure(object):
    """
    Structure contains all the machinery necessary to build a MOF from
    simple organic and metallic building units

    """
 
    def __init__(self):
        self.cell = np.zeros((3,3))
        self.acell = np.zeros(6)
        self.icell = np.zeros((3,3))
        # the list mof keeps track of all the SBUs
        self.mof = []
        # initialize logging for this MOF
        self.log = bookkeeping.Log()

    def build_one(self):
        """Adds one extra SBU to the growing MOF"""
        if len(self.mof) == 0:
            self.mof = [SBU("Fe", pdbfile="Oh.pdb")]
            dump = coordinate_dump(self.mof)
            write_xyz("history", dump[1], dump[0])
        else:
            self.add_new_sbus()

    def build(self):
        """Exhaustively builds MOF with metals and linkers"""
    
        self.mof = [SBU("Fe", pdbfile="Oh.pdb")]
        dump = coordinate_dump(self.mof)
        write_xyz("history", dump[1], dump[0])
   
        # progressive building of MOF, will complete when all
        # bonds are saturated
        done = False
        while not done:
            for num, itstruct in enumerate(self.mof):
                # adds bonds via proximity or periodic boundary
                self.sbu_check(num)

            # adds one set of new sbus to all unsaturated bonds
            self.add_new_sbus()
            # check full saturation condition
            # TODO(pboyd): add a routine to shift all atoms into the
            # box
            if self.saturated() and self.complete_box():
                done = True

        # dump framework into an xyz file
        dump = coordinate_dump(self.mof)
        write_xyz("RANDOM", dump[1], dump[0])
        write_pdb("RANDOM", dump[1], dump[0], self.acell)

    def saturated(self):
        """
        returns False if at least one bond is unsaturated
        """
        for imof in self.mof:
            for ibond in imof.connectivity:
                if ibond is None:
                    return False

        return True

    def sbu_check(self, isbu):
        """
        Checks an SBU for 
        1. if the SBU is outside the box 
        2. bonding locally,
        3. bonding by existing pbc
        4. if there are atomistic overlaps
        """
        for nbnd, btype in enumerate(self.mof[isbu].connectivity):
            if btype is None:
                # Now test for connecting with other SBUs
                for nsbu in range(len(self.mof)):

                    if self.local_bond(isbu, nbnd, nsbu):
                        #2. Local bond exists
                        self.mof[isbu].connectivity[nbnd] = nsbu
                    elif self.remote_bond(isbu, nbnd, nsbu):
                        #3. remote bond via PBC
                        self.mof[isbu].connectivity[nbnd] = nsbu

    def add_new_sbus(self):
        """adds a new sbus to the growing MOF"""
        
        sbus = len(self.mof)
        
        for sbu1 in range(sbus):
            for bond1, bondtype in enumerate(self.mof[sbu1].connectivity):
                if bondtype is None:
                    if self.mof[sbu1].ismetal:
                        # add organic linker to the metal link
                        self.mof.append(SBU("benzene", 
                            pdbfile="bdc.pdb",
                            ismetal=False))
                    else:
                        # add metal corner to the organic link
                        self.mof.append(SBU("Fe", 
                            pdbfile="Oh.pdb"))
    
                    sbu2 = len(self.mof) - 1
                    # randomly chose the linker bond to attach to the 
                    # metal
                    bond2 = randrange(
                            len(self.mof[sbu2].connect_vector))
                    self.log.info(
                            "added %s, sbu %i, bond %i to SBU %i, %s bond %i"
                            %(self.mof[sbu2].type, sbu2, bond2, sbu1,
                              self.mof[sbu1].type, bond1))
                    self.mof[sbu1].connectivity[bond1] = sbu2 
                    self.mof[sbu2].connectivity[bond2] = sbu1
                    self.sbu_align(sbu1, bond1, sbu2, bond2)
                    # rotate by Y vector
                    self.bond_align(sbu1, bond1, sbu2, bond2) 
                    # check to see if the new SBU will bond to
                    # exisiting SBUs
                    self.sbu_check(sbu2)

                    dump = coordinate_dump(self.mof)
                    write_xyz("history", dump[1], dump[0])

    def local_bond(self, sbu1, bond1, sbu2):
        """Checks if two SBUs can join by local bonding"""
        point = self.mof[sbu1].connectpoints[bond1]
        vector = self.mof[sbu1].connect_vector[bond1]

        # No local bond if both SBUs are the same type
        if self.mof[sbu1].ismetal == self.mof[sbu2].ismetal:
            return False

        for nbnd, type in enumerate(self.mof[sbu2].connectivity):
            if type is None:
                # Check if vectors are anti-parallel
                point2 = self.mof[sbu2].connectpoints[nbnd] - \
                        self.mof[sbu2].connect_vector[nbnd]
                vector2 = self.mof[sbu2].connect_vector[nbnd]
                test1 = anti_parallel_test(vector, vector2) and \
                        not self.pointing_away(sbu1, bond1, sbu2, nbnd)

                # Check if connectpoints overlaps with atom
                # with a given tolerance
                test2 = points_close(point, point2, 1.)
                if test1 and test2:
                    # add sbu1 to the sbu2 bond
                    self.mof[sbu2].connectivity[nbnd] = sbu1
                    return True

        return False

    def remote_bond(self, sbu1, bond1, sbu2):
        """Checks if two SBUs can join by a pbc"""
        vector = self.mof[sbu1].connect_vector[bond1]
        if self.mof[sbu1].ismetal == self.mof[sbu2].ismetal:
            return False

        for nbnd, type in enumerate(self.mof[sbu2].connectivity):
            if type is None:
                vector2 = self.mof[sbu2].connect_vector[nbnd]
                test = anti_parallel_test(vector, vector2) and \
                       self.pointing_away(sbu1, bond1, sbu2, nbnd)
                if (test):
                    # vectors are anti-parallel check if they
                    # can be connected by an existing periodic
                    # boundary or generate a new one
                    if self.box_gen(sbu1, bond1, sbu2, nbnd):
                        self.mof[sbu2].connectivity[nbnd] = sbu1
                        return True
                    
        return False
    
    def box_gen(self, sbu1, bond1, sbu2, bond2):
        """
        Routine to tackle the enormous problem of the 
        periodic boundary conditions
        """
        vector = self.mof[sbu2].connectpoints[bond2] - \
                 (self.mof[sbu1].connectpoints[bond1] - 
                 self.mof[sbu1].connect_vector[bond1])

        # check vector for parallel with existing cell vectors
        # larger, smaller than existing vectors
        indx = None
        for icell in range(3):
            if parallel_test(self.cell[icell], vector) or \
                    anti_parallel_test(
                            self.cell[icell], vector):
                if np.allclose(length(self.cell[icell]), 
                                 length(vector)):
                    # bonding!
                    return True
                elif length(vector) > length(self.cell[icell]):
                    # FIXME(pboyd): add some problem solving here
                    pass
                elif length(vector) < length(self.cell[icell]):
                    # No bonding
                    return False
            else:
                # create a new PBC
                if indx is None:
                    if sum(self.cell[icell]) == 0.:
                        indx = icell
        if indx is not None:
            self.cell[indx] = vector[:]
        return True

    def pointing_away(self, sbu1, bond1, sbu2, bond2):
        """
        Checks if two vectors are pointing away from each other
        """
        tail1 = self.mof[sbu1].connectpoints[bond1] -\
                self.mof[sbu1].connect_vector[bond1]
        head1 = self.mof[sbu1].connectpoints[bond1]

        tail2 = self.mof[sbu2].connectpoints[bond2] -\
                self.mof[sbu2].connect_vector[bond2]
        head2 = self.mof[sbu2].connectpoints[bond2]

        if (length(tail1, tail2) < length(head1, head2)):
            return True
        
        return False
    
    def overlap(self, istruct):
        """
        checks if one of the SBUs in structure (istruct) is
        overlapping with the other SBUs
        """
        # distance tolerance in angstroms
        disttol = 3.
    
        overlapping = False
        for nxyz, xyz in enumerate(self.mof[istruct].coordinates):
            exclrange = range(len(self.mof))
            exclrange.remove(istruct)
            for struct in exclrange:
                for ncoord,coords in enumerate(self.mof[struct].coordinates):
                    if (length(xyz, coords) <= disttol):
                        overlapping = True

    def bond_align(self, sbu1, bond1, sbu2, bond2):
        """
        Align two sbus by their orientation vectors
        """

        sbu1 = self.mof[sbu1]
        sbu2 = self.mof[sbu2]

        axis = rotation_axis(sbu1.anglevect[bond1],
                             sbu2.anglevect[bond2])
        angle = calc_angle(sbu1.anglevect[bond1],
                           sbu2.anglevect[bond2])
        # origin of rotation
        rotpoint = sbu2.centre_of_mass()
        # rotation
        self.sbu_rotation(sbu2, rotpoint, axis, angle)


    def sbu_align(self, sbu1, bond1, sbu2, bond2):
        """
        Align two sbus, were sbu1 is held fixed and
        sbu2 is rotated and shifted to bond with sbu1
        """

        # the order of the cross product is important to 
        # signify which vector remains stationary (statvect)
        # and which will rotate.  This axis will determine the
        # proper rotation when applying the rotation matrix

        # Change sbu1 and sbu2 from indices to reference to the
        # class
        sbu1 = self.mof[sbu1]
        sbu2 = self.mof[sbu2]

        axis = rotation_axis(sbu2.connect_vector[bond2], 
                             -1.*sbu1.connect_vector[bond1])
        angle = calc_angle(sbu1.connect_vector[bond1], 
                           sbu2.connect_vector[bond2])
        # origin of rotation
        rotpoint = sbu2.centre_of_mass()

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
        C = np.zeros(3)
        sbu.coordinates = array(sbu.coordinates * R) + C
        sbu.connect_vector = array(sbu.connect_vector * R) + C
        sbu.connectpoints = array(sbu.connectpoints * R) + C
        sbu.anglevect = array(sbu.anglevect * R) + C
        
        # recursive scan of all connecting SBUs, if they connect back to 
        # the stat_SBU, do not perform a rotation.
        # FIXME(pboyd): this doesn't work
    
#        dir = []
#        dir = self.recurse_bonds(rotstruct,[stat_SBU],stat_SBU)
#        dir.append(rotstruct)
#    
#        # apply rotation to the SBUs listed in dir
#        for i in dir:
#            self.mof[i].coordinates = \
#                    array(self.mof[i].coordinates * R) + C
#            self.mof[i].connect_vector = \
#                    array(self.mof[i].connect_vector * R) + C
#            self.mof[i].connectpoints = \
#                    array(self.mof[i].connectpoints * R) + C
#            self.mof[i].anglevect = \
#                    array(self.mof[i].anglevect * R) + C
    
    
    def recurse_bonds(self, bond, frombond, xbond):
        """
        Returns a list of indices of all bonded SBUs and their bonded
        SBUs and so on..
    
        entries are as follows:
        [branching point SBU, the SBU one level up the tree, 
        the termination SBU, the list of SBUs]
        """
        # FIXME(pboyd): not fully convinced this works properly
        if self.mof[bond].connectivity is None:
            return
    
        # frombond stores a list of SBUs already included in the bonding.
        # this is such that the recursion list doesn't repeat itself with
        # loops
        r = []
        frombond.append(bond)
        for i in self.mof[bond].connectivity:
            # TODO(pboyd): check for bonding to the xbond (if this happens,
            # kill the recursion and the rotation)
            if (i is not None)and(i not in frombond)and(i is not xbond):
                r.append(i)
                for tmp in self.recurse_bonds(i, frombond, xbond):
                    r.append(tmp)
        return r
    
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
    if("O" in [atom1,atom2])and("C" in [atom1,atom2]):
        return 1.6
    elif("Zn" in (atom1,atom2))and("X" in (atom1,atom2)):
        return 2.2
    elif("Fe" in (atom1,atom2))and("X" in (atom1,atom2)):
        return 2.2
    elif("X" in (atom1,atom2))and("O" in (atom1,atom2)):
        # temporary correction for carboxylate vector
        return 0.
    elif("X" in (atom1,atom2))and("C" in (atom1,atom2)):
        # TODO(pboyd): these conditions are not robust and
        # should be changed when the code becomes bigger
        return 2.6
    elif("Y" in (atom1,atom2))and("X" not in (atom1,atom2)):
        return 0.
    elif("Y" in (atom1,atom2))and("X" in (atom1,atom2)):
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
    vector = coord2 - coord1
    return np.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + 
                   vector[2] * vector[2])


def rotation_axis(vect1, vect2):
    """
    axis of rotation for two vectors (90 degrees from both)
    """
    vect1 = normalize(vect1)
    vect2 = normalize(vect2)
    return normalize(np.cross(vect1, vect2))

def calc_angle(vect1, vect2):
    """ determines angle between vector1 and vector2"""
    dot12 = np.dot(vect1, vect2)
    dist11 = np.sqrt(np.dot(vect1, vect1))
    dist22 = np.sqrt(np.dot(vect2, vect2))
    
    return np.arccos(dot12 / (dist11 * dist22))

def rotation_matrix(axis, angle):
    """
    returns a (3,3) rotation matrix based on
    axis and angle
    the rotation is counter-clockwise when looking
    down the axis so choice of sign is important here
    """
    ux = axis[0]
    uy = axis[1]
    uz = axis[2]

#    q = array([np.cos(angle / 2.),
#               np.sin(angle / 2.) * np.cos(ux),
#               np.sin(angle / 2.) * np.cos(uy),
#               np.sin(angle / 2.) * np.cos(uz)])
#    q = normalize(q)
# 
#    matrix = np.matrix([[q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3],
#                         2.*(q[1]*q[2] - q[0]*q[3]),
#                         2.*(q[0]*q[2] + q[1]*q[3])],
#                        [2.*(q[1]*q[2] + q[0]*q[3]),
#                         q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3],
#                         2.*(q[2]*q[3] - q[0]*q[1])],
#                        [2.*(q[1]*q[3] - q[0]*q[2]),
#                         2.*(q[0]*q[1] + q[2]*q[3]),
#                         q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3]]])
    matrix = np.matrix([[
        np.cos(angle) + ux * ux * (1 - np.cos(angle)),
        ux * uy * (1 - np.cos(angle)) - uz * np.sin(angle),
        ux * uz * (1 - np.cos(angle)) + uy * np.sin(angle)],
        [
        uy * ux * (1 - np.cos(angle)) + uz * np.sin(angle),
        np.cos(angle) + uy * uy * (1 - np.cos(angle)),
        uy * uz * (1 - np.cos(angle)) - ux * np.sin(angle)],
        [
        uz * ux * (1 - np.cos(angle)) - uy * np.sin(angle),
        uz * uy * (1 - np.cos(angle)) + ux * np.sin(angle),
        np.cos(angle) + uz * uz * (1 - np.cos(angle))]])

    return matrix

def planar_test(coordinates):
    """
    Tests whether four or more points are within a plane.
    You can run the test for fewer coordinates
    but you should note that you might be mentally
    handicapped if you do so.
    """

    # Define a plane by three points
    planevector1 = normalize(coordinates[1] - coordinates[0])

    if(np.allclose(sum(planevector1),0.)):
        # TODO(pboyd): raise error, two atoms are in the exact same position
        pass    
    # test points for colinearity with the first vector
    for coord in coordinates[2:]:
        if(not linear_test(array([coordinates[0],coordinates[1],coord]))):
            planevector2 = normalize(coord - coordinates[0])
            break

    # test to see if planevector2 is allocated (if there is a plane
    # within the coordinates) otherwise co-linear
    try:
        planevector2
    except:
        return False

    for coord in coordinates[2:]:
       # If not colinear test for planarity
        newvect = normalize(coord - coordinates[0])

        # test = 0 if co-planar
        test = np.inner(newvect, np.cross(planevector1, planevector2))
        # the tolerance is quite low, this has to do with the number
        # of sig. figs. of the training set.
        if (not np.allclose(0., test, atol=1e-2)):
            return False
    
    return True
        

def normalize(vector):
    """changes vector length to unity"""
    return vector / np.sqrt(np.dot(vector, vector))

def linear_test(coordinates):
    """
    Tests for three or more points lying within the same
    line.  You can perform this test for fewer points,
    but you should note that you might be mentally
    handicapped if you do so.
    """
    
    # test by cross-product.  If linear, should be the zero vector
    vector1 = normalize(coordinates[1] - coordinates[0])

    for point in coordinates[2:]:
        vector2 = normalize(point - coordinates[0])
        crossprod = np.cross(vector2,vector1)
        if np.allclose(crossprod, np.zeros(3), atol=1e-2):
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
        tol = 1e-3

    dist = length(point1, point2)
    if dist <= tol:
        return True
    return False

def anti_parallel_test(vector1, vector2):
    """test for two vectors being in anti-parallel orientations"""

    tol = 2.e-2
    vector1 = normalize(vector1)
    vector2 = normalize(vector2)
    test1 = np.allclose(np.cross(vector1, vector2), 
                        np.zeros(3), atol=tol)
    test2 = np.allclose(np.dot(vector1, vector2), -1., atol=tol)
    if test1 and test2:
        return True
    else:
        return False

def parallel_test(vector1, vector2):
    """test for two vectors being in parallel orientations"""

    tol = 2.e-2
    vector1 = normalize(vector1)
    vector2 = normalize(vector2)
    test1 = np.allclose(np.cross(vector1, vector2), 
                        np.zeros(3), atol=tol)
    test2 = np.allclose(np.dot(vector1, vector2), 1., atol=tol)
    if test1 and test2:
        return True
    else:
        return False

def write_pdb(label, atoms, coords, acell):
    pdbfile = open('%s.pdb' % label, 'w')
    today = date.today()
    lines = []
    lines.append("%6s   Created: %s \n%-6s" % ('REMARK', 
                  today.strftime("%A %d. %B %Y"), 'CRYST1'))
    lines.append("%8.3f %8.3f %8.3f %6.2f %6.2f %6.2f P1\n"
                  % (tuple(acell)))
    for atom in range(len(atoms)):
        lines.append("%-6s%4i %3s MOL%6i" %\
                     ('ATOM', atom, atoms[atom], 1,) +\
                     5*" " + "%8.3f %8.3f %8.3f  1.00  0.00 " 
                     % tuple(coords[atom]) + \
                     "%11s\n"%(atoms[atom]))
    pdbfile.writelines(lines)
    pdbfile.close()




def write_xyz(label, atoms, coords):
    xyzfile = open('%s.xyz' % label, 'a')
    xyzfile.write('%i\ncoordinate dump\n'% len(coords))
    for i in range(len(coords)):
        xyzfile.write('%s%12.5f%12.5f%12.5f\n' % (atoms[i],
                      coords[i][0], coords[i][1], coords[i][2]))

def coordinate_dump(structure):
    """ 
    Dumps all the atom labels and xyz coordinates from a list
    of classes"""

    coords = []
    atoms = []

    for itstruct in range(len(structure)):
        for itsite in range(len(structure[itstruct].coordinates)):
            coords.append(structure[itstruct].coordinates[itsite])
            atoms.append(structure[itstruct].atomlabel[itsite])

    return (array(coords), atoms)

def main():
    """Default if run as an executable"""

    struct = Structure()
    struct.build()
if __name__ == '__main__':
    main()
