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

class SBU(object):
    """
    """

    def __init__(self, type, sitelabel, coordinates):
        self.type = type
        self.sitelabel = sitelabel

        # TODO(pboyd): shift the coordinates by the COM such that
        # they are centred at the origin
        self.coordinates = array(coordinates)
        self.mass = 0. 
        self.numb = 0 
        self.connect = []
        self.nbonds = [] 
        self.bonding = []
        self.unsat = []
        self.links = []
        self.bondpoints = []
        # populate arrays
        # structcalc is a tuple of important connectivity information
        # of the linker
        structcalc = coord_matrix(self.sitelabel, self.coordinates)
        self.connect = structcalc[0]
        self.bonding = structcalc[1]
        self.nbonds = structcalc[2]

        # calculate centre of mass      
        self.COM = self.centre_of_mass()

        # shift all coordinates back by the centre of mass
        self.COM_shift()

        # check unsaturated bonds for possible links
        self.check_unsat()

        # add linking points to the SBU.
        # TODO(pboyd): add possible linking types
        self.add_links()


    def check_unsat(self):
        """
        checking for unsaturated bonds
        """
        unsat = []
        # need to improve robustness of this code
        for i in range(len(self.nbonds)):
            if ((self.nbonds[i] <= 2)and(self.sitelabel[i] != "H")and
                (self.sitelabel[i] != "X")):
                unsat.append(i)
        self.unsat = unsat

    def add_links(self):
        """
        checks for certain unsaturated atoms and adds links accordingly
        """

        purgeatoms = []
        for indx, atom in enumerate(self.sitelabel):
            # test for carboxylate
            if atom == "C":
                carb_test = self.carboxylate_test(indx)
                if carb_test[0]:
                    pt = carb_test[1][0] + carb_test[1][1] \
                         - self.coordinates[indx]
                    vect = carb_test[1][0] + carb_test[1][1] \
                           - 2 * self.coordinates[indx]

                    self.bondpoints.append(pt)
                    self.links.append(vect)
            # test for "X" atoms (considered a linking point)
            elif atom == "X":
                self.bondpoints.append(self.coordinates[indx])
                # there should only be one bonding atom to an "X"
                bondatm = self.bonding[indx][0]
                vect = self.coordinates[indx] - self.coordinates[bondatm]
                self.links.append(vect)
                purgeatoms.append(indx)

            # TODO(pboyd): add other tests for different 
            # coordinating atoms

        # purge any coordinates in self.coordinates belonging to "X"'s
        # as they are now part of self.links and self.bonding
        self.purge(purgeatoms)

        self.links = array(self.links)
        self.bondpoints = array(self.bondpoints)

    def purge(self, atoms):
        """Remove entries for X atoms"""
        for patm in reversed(atoms):
            self.coordinates = np.delete(self.coordinates, patm, 0) 
            self.sitelabel.pop(patm)
            self.nbonds.pop(patm)
            self.bonding.pop(patm)

        self.connect = coord_matrix(self.sitelabel, self.coordinates)[0]

    def carboxylate_test(self, atom):
        """Simple test for a carboxylate group"""
        oxcount = 0
        bondvector =[]

        for iatm in self.bonding[atom]:
            if (self.sitelabel[iatm] == "O") and (iatm in self.unsat):
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

        for i in range(len(self.sitelabel)):
            top += self.coordinates[i] * WEIGHT[self.sitelabel[i]]
            bottom += WEIGHT[self.sitelabel[i]]
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

def align(coords, alignvect, alignpt, statvect, statpt):
    """
    align two vectors and adjust the coordinates
    This is done by an axis-angle rotation then a shift
    """

    # the order of the cross product is important to 
    # signify which vector remains stationary (statvect)
    # and which will rotate.  This axis will determine the
    # proper rotation when applying the rotation matrix

    axis = np.cross(statvect, alignvect)
    axis = axis / np.sqrt(np.dot(axis, axis))

    angle = calc_angle(statvect, alignvect)

    # Transformation matrix for rotation
    transform = rotation_matrix(axis, angle)

    # array of shifted coordinates
    shiftcoords = np.zeros((len(coords),3))

    # shift the alignment vector by the rotation matrix
    alignvect = array(alignvect * transform)[0]
    alignpt = array(alignpt * transform)[0]
    shiftvector = (statpt) - (alignpt - alignvect)

    for ixyz, cart in enumerate(coords):
        shiftcoords[ixyz] = array(cart * transform)[0]
        shiftcoords[ixyz] = shiftcoords[ixyz] + shiftvector

    return shiftcoords


def rotation_axis(vect1, vect2):
    """
    axis of rotation for two vectors (90 degrees from both)
    """
    non_norm = np.cross(vect1, vect2)
    return non_norm / np.sqrt(np.dot(non_norm, non_norm))

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
    return vector / np.dot(vector, vector)

def linear_test(coordinates):
    """
    Tests for three or more points lying within the same
    line.  You can perform this test for fewer points,
    but you should note that you might be mentally
    handicapped if you do so.
    """
    
    # test by cross-product.  If linear should be the zero vector
    vector1 = normalize(coordinates[1] - coordinates[0])

    for point in coordinates[2:]:
        vector2 = normalize(point - coordinates[0])
        crossprod = np.cross(vector2,vector1)
        if np.allclose(crossprod, np.zeros(3), atol=1e-2):
            return True

    return False


def write_xyz(label, atoms, coords):
    xyzfile = open('%s.xyz' % label, 'w')
    xyzfile.write('%i\ncoordinate dump\n'% len(coords))
    for i in range(len(coords)):
        xyzfile.write('%s%12.5f%12.5f%12.5f\n' % (atoms[i],
                coords[i][0], coords[i][1], coords[i][2]) )

def main():
    """Default if run as an executable"""

    linkr = SBU("BTC", sample.atom_labels[0], sample.atom_coordinates[0])
    metal = SBU("Zn", sample.atom_labels[1], sample.atom_coordinates[1])

    benzene = SBU("Benzene", sample.atom_labels[2], 
                  sample.atom_coordinates[2])
    anthracene = SBU("Anthracene", sample.atom_labels[3],
                     sample.atom_coordinates[3])
    phosphate = SBU("Phosphate", sample.atom_labels[4],
                     sample.atom_coordinates[4])

    # TODO(pboyd): add random choice of link and match metal with linker

    for i in range(len(metal.links)):
        testings = align(phosphate.coordinates, phosphate.links[0],
                     phosphate.bondpoints[0], -1*metal.links[i],
                     metal.bondpoints[i])
        write_xyz("rotation_test%s"%str(i), phosphate.sitelabel, testings)
    write_xyz("original", phosphate.sitelabel, phosphate.coordinates)
    write_xyz("metal", metal.sitelabel, metal.coordinates)

if __name__ == '__main__':
    main()

