#!/usr/bin/env python

"""
GenStruct -- Generation of Structures for fapping.
"""

__version__ = "$Revision$"

import subprocess
import sys
import numpy as np
from numpy import array
from elements import WEIGHT, ATOMIC_NUMBER
import sample
from random import random, uniform, randrange, choice

class SBU(object):
    """
    """

    def __init__(self, type, atomlabel, coordinates, ismetal=None):
        self.type = type
        self.atomlabel = atomlabel
        self.coordinates = array(coordinates)
        self.mass = 0. 
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


    def check_unsat(self):
        """
        checking for unsaturated bonds
        """
        unsat = []
        # need to improve robustness of this code
        for i in range(len(self.nbonds)):
            if ((self.nbonds[i] <= 2)and(self.atomlabel[i] != "H")and
                (self.atomlabel[i] != "X")):
                unsat.append(i)
        self.unsat = unsat

    def add_connect_vector(self):
        """
        checks for certain unsaturated atoms and adds connect_vector accordingly
        """

        purgeatoms = []
        for indx, atom in enumerate(self.atomlabel):
            # test for carboxylate
            if atom == "C":
                carb_test = self.carboxylate_test(indx)
                if carb_test[0]:
                    pt = carb_test[1][0] + carb_test[1][1] \
                         - self.coordinates[indx]
                    vect = carb_test[1][0] + carb_test[1][1] \
                           - 2 * self.coordinates[indx]

                    self.connectpoints.append(pt)
                    self.connect_vector.append(vect)
                    self.connectivity.append(None)
            # test for "X" atoms (considered a linking point)
            elif atom == "X":
                self.connectpoints.append(self.coordinates[indx])
                # there should only be one bonding atom to an "X"
                bondatm = self.bonding[indx][0]
                vect = self.coordinates[indx] - self.coordinates[bondatm]
                self.connect_vector.append(vect)
                self.connectivity.append(None)
                purgeatoms.append(indx)

            # TODO(pboyd): add other tests for different 
            # coordinating atoms

        # purge any coordinates in self.coordinates belonging to "X"'s
        # as they are now part of self.connect_vector and self.connectpoints
        self.purge(purgeatoms)

        self.connect_vector = array(self.connect_vector)
        self.connectpoints = array(self.connectpoints)

    def purge(self, atoms):
        """Remove entries for X atoms"""
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

    def align(self, ibond, statvect, statpt):
        """
        align two vectors and adjust the coordinates
        This is done by an axis-angle rotation then a shift
        """

        # the order of the cross product is important to 
        # signify which vector remains stationary (statvect)
        # and which will rotate.  This axis will determine the
        # proper rotation when applying the rotation matrix
        axis = rotation_axis(statvect, self.connect_vector[ibond])

        angle = calc_angle(statvect, self.connect_vector[ibond])

        # Transformation matrix for rotation
        transform = rotation_matrix(axis, angle)

        # shift bonding information by the rotation matrix and shift vectors
        self.connect_vector = array(self.connect_vector * transform)
        self.connectpoints = array(self.connectpoints * transform)
        # shiftvector is found after initial rotation of the bond vectors
        # and bond points
        shiftvector = (statpt) - (self.connectpoints[ibond] - self.connect_vector[ibond])
        # bond points are then shifted by the shiftvector
        self.connectpoints += shiftvector

        # coordinate rotation and translation
        self.coordinates = array(self.coordinates * transform)
        self.coordinates += shiftvector


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
    else:
        return 1.6

def length(coord1, coord2):
    """ returns the length between two vectors"""
    vector = coord2 - coord1
    return np.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + 
                   vector[2] * vector[2])


def rotation_axis(vect1, vect2):
    """
    axis of rotation for two vectors (90 degrees from both)
    """
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

def rotation(rotbond, rotstruct, structure):
    """
    Takes a growing structure and rotates the coordinates
    about a selected bond.  If other SBU's are attached
    to the structure, they are rotated too unless they 
    connect to the stationary SBU.
    This only applies to linkers with a single bond to the 
    SBU - bidentate bonds must have other considerations
    """

    axis = normalize(structure[rotstruct].connect_vector[rotbond])

    # random rotation for now...
    angle = uniform(0,2*np.pi)

    R = rotation_matrix(axis, angle)

    # find out what the stationary SBU is
    stat_SBU = structure[rotstruct].connectivity[rotbond]

    # correction for rotation about the origin
    # c * (I - R)
    c = structure[rotstruct].connectpoints[rotbond] + \
            structure[rotstruct].connect_vector[rotbond]

    C = array(c * (np.identity(3) - R))

    # recursive scan of all connecting SBUs, if they connect back to 
    # the stat_SBU, do not perform a rotation.

    dir = []
    dir = recurse_bonds(rotstruct,[stat_SBU],stat_SBU,structure)
    dir.append(rotstruct)

    # apply rotation to the SBUs listed in dir
    for i in dir:
        structure[i].coordinates = \
                array(structure[i].coordinates * R) + C
        structure[i].connect_vector = \
                array(structure[i].connect_vector * R) + C
        structure[i].connectpoints = \
                array(structure[i].connectpoints * R) + C


def recurse_bonds(bond, frombond, xbond, structure):
    """
    Returns a list of indices of all bonded SBUs and their bonded
    SBUs and so on..

    entries are as follows:
    [branching point SBU, the SBU one level up the tree, 
    the termination SBU, the list of SBUs]
    """
    # FIXME(pboyd): not fully convinced this works properly
    if not structure[bond].connectivity:
        return

    # frombond stores a list of SBUs already included in the bonding.
    # this is such that the recursion list doesn't repeat itself with
    # loops
    r = []
    frombond.append(bond)
    for i in structure[bond].connectivity:
        # TODO(pboyd): check for bonding to the xbond (if this happens,
        # kill the recursion and the rotation)
        if (i is not None)and(i not in frombond)and(i is not xbond):
            r.append(i)
            for tmp in recurse_bonds(i, frombond, xbond, structure):
                r.append(tmp)
    return r

def build():
    """Randomly builds MOF with metals and linkers"""

    organic = SBU("BTC", sample.atom_labels[0][:], 
            sample.atom_coordinates[0][:], False)

    benzene = SBU("Benzene", sample.atom_labels[2][:], 
            sample.atom_coordinates[2][:], False)

    anthracene = SBU("Anthracene", sample.atom_labels[3][:],
            sample.atom_coordinates[3][:], False)

    phosphate = SBU("Phosphate", sample.atom_labels[4][:],
            sample.atom_coordinates[4][:], False)

    # Create an initial instance of metal
    # more instances will be created as the system gets
    # bigger.

    structure = [SBU("Zn", sample.atom_labels[1][:], 
        sample.atom_coordinates[1][:])]

    count = 0
    done = False
    while not done:
        count += 1 
        # chose a random SBU
        randchoice = randrange(len(structure))

        # check for unsaturated connection points
        bondchoice = []
        for itbond, bonded in enumerate(
                structure[randchoice].connectivity):
            if bonded is None:
                bondchoice.append(itbond)

        # check if the SBU is fully saturated
        if bondchoice == []:
            # TODO(pboyd): include while loop to make sure we select
            # an unsaturated molecule
            pass
        # chose bond to mess with
        else:
            bondrand = choice(bondchoice)

        # the value of the connectivity of the nth bond is
        # the index of the added structure
        structure[randchoice].connectivity[bondrand] = len(structure) 

        link = structure[randchoice].connect_vector[bondrand] * -1.
        point = structure[randchoice].connectpoints[bondrand]

        if structure[randchoice].ismetal:
            # add organic linker to the metal link[randchoice]
            iorg = len(structure)

            structure.append(SBU("BTC", sample.atom_labels[0][:],
                             sample.atom_coordinates[0][:], False))

            # randomly chose the linker bond to attach to the metal
            jrand = randrange(len(structure[iorg].connect_vector))

            structure[iorg].connectivity[jrand] = randchoice 
            structure[iorg].align(jrand, link, point)

        else:
            # add metal corner to the organic link[randchoice]
            imet = len(structure)

            structure.append(SBU("Zn", sample.atom_labels[1][:],
                sample.atom_coordinates[1][:]))

            # randomly chose the linker bond to attach to the organic
            # linker 
            jrand = randrange(len(structure[imet].connect_vector))

            structure[imet].connectivity[jrand] = randchoice
            structure[imet].align(jrand, link, point)

        if count == 5:
            done = True

    dump = coordinate_dump(structure)

    write_xyz("RANDOM", dump[1], dump[0])

    # call a random rotation
    randrot = randrange(0, len(structure))

    s = []
    for k, ibond in enumerate(structure[randrot].connectivity):
        if ibond is not None:
            s.append(k)

    randbond = choice(s)
    rotation(randbond, randrot, structure)

    # TODO(pboyd): assign a bond length to the connection site, track
    # which metals are bonded to which organic linkers (for larger 
    # molecular rotations?)
    # TODO(pboyd): implement a periodic boundary test
    # TODO(pboyd): implement a test for a fully saturated system
    # TODO(pboyd): scan for atom - atom overlap

    dump = coordinate_dump(structure)

    write_xyz("RANDOM_ROT", dump[1], dump[0])

def write_xyz(label, atoms, coords):
    xyzfile = open('%s.xyz' % label, 'w')
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

    # TODO(pboyd): write all actions out to a log.
    build()

if __name__ == '__main__':
    main()

