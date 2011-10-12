#!/usr/bin/env python

"""
GenStruct -- Generation of Structures for fapping.
"""

__version__ = "$Revision$"

import subprocess
import sys
import numpy
from numpy import array
from elements import WEIGHT, ATOMIC_NUMBER
import sample

class SBU(object):
    """
    """

    def __init__(self, type, sitelabel, coordinates, mass=None, numb=None, connect=None):
        self.type = type
        self.coordinates = array(coordinates)
        self.mass = mass
        self.numb = numb
        self.connect = connect 
        self.nbonds = None 
        self.bonding = None
        self.sitelabel = sitelabel
        self.unsat = []
        self.links = []
        self.bondpoints = []

    def check_unsat(self):
        """
        checking for unsaturated bonds
        """
        unsat = []
        # need to improve robustness of this code
        for i in range(len(self.nbonds)):
            if ((self.nbonds[i] <= 2)and(self.sitelabel[i] != "H")):
                unsat.append(i)
        self.unsat = unsat

    def add_links(self):
        """
        checks for certain unsaturated atoms and adds links accordingly
        """

        bondpoint, bondvect = [], []
        for indx, atom in enumerate(self.sitelabel):
            # test for carboxylate
            if atom == "C":
                carb_test = self.carboxylate_test(indx)
                if carb_test[0]:
                    pt = carb_test[1][0] + carb_test[1][1] \
                         - self.coordinates[indx]
                    vect = carb_test[1][0] + carb_test[1][1] \
                           - 2 * self.coordinates[indx]

                    bondpoint.append(pt)
                    bondvect.append(vect)
            #TODO(pboyd): add other tests for different 
            # coordinating atoms
        self.links = bondvect
        self.bondpoints = bondpoint
                    

    def carboxylate_test(self, atom):
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

def coord_matrix(coord):
    """
    generate a coordination matrix for a list of x,y,z coordinates
    """
    connect = numpy.zeros((len(coord),len(coord)))
    bonding = []
    numbonds = []
    for i in range(len(coord)):
        bonding.append([])
        numbonds.append(0)
        for j in range(len(coord)):
            dist = length(coord[i],coord[j])
            #TODO(pboyd): add a tolerance for the distance, depending on atom-atom interaction
            tol = 1.6
            if (i != j) and (dist <= tol):
                connect[i,j] = dist
                connect[j,i] = dist
                numbonds[i] += 1
                bonding[i].append(j)
            else:
                connect[i,j] = 0.
                connect[j,i] = 0.
    return (connect, bonding, numbonds)

def length(coord1, coord2):
    vector = coord2 - coord1
    return numpy.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2])

def align():
    """
    align two vectors and adjust the coordinates
    """
def write_xyz(label, atoms, coords):
    xyzfile = open('%s.xyz' % label, 'w')
    xyzfile.write('%i\ncoordinate dump\n'% len(coords))
    for i in range(len(coords)):
        xyzfile.write('%s%12.5f%12.5f%12.5f\n' % (atoms[i],
                coords[i][0], coords[i][1], coords[i][2]) )
    
def main():
    """Default if run as an executable"""

    metal = SBU("Zn",65.39,30)
    linkr = SBU("BTC", sample.atom_labels, sample.atom_coordinates)
#   put all of the following lines within SBU

    structcalc = coord_matrix(linkr.coordinates)
    linkr.connect = structcalc[0]
    linkr.bonding = structcalc[1]
    linkr.nbonds = structcalc[2]
    COM = linkr.centre_of_mass()
    linkr.check_unsat()
    linkr.add_links()

    write_xyz(linkr.type, linkr.sitelabel, linkr.coordinates)
    write_xyz("bond_points", ["X"]* len(linkr.bondpoints), \
            linkr.bondpoints)

if __name__ == '__main__':
    main()

