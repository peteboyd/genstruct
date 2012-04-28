#!/usr/bin/env python
import numpy as np
from optparse import OptionParser
from operations import *
from elements import WEIGHT
class File(object):
    def __init__(self, file):
        self.bondingatoms = []
        self.coordinates = []
        self.connect_vector = []
        self.anglevect = []
        self.connectpoints = []
        self.connectangles = []
        self.atomlabel = []
        self.connect = []
        self.bonding = []
        self.nbonds = []
        self.symmetrytype = []
        self.special_bond = []
        self.from_pdb(file)

    def from_pdb(self, filename):
        """Reads atom coordinates of an SBU from a pdb file"""

        filename = filename + ".pdb"
        pdbfile = open(filename)
        pdblines = pdbfile.readlines()
        for line in pdblines:
            if line.lower().startswith('remark'):
                if line[7:15].strip().lower().startswith('index'):
                    self.index = int(line[16:].strip())
                if line[7:15].strip().lower().startswith('name'):
                    self.name = line[16:].strip()
            if line.lower().startswith('name'):
                self.name = line[5:].strip()
            if line.lower().startswith('atom'):
                label = line[12:16].strip()
                self.atomlabel.append(label)
                self.coordinates.append([float(line[30:38]), 
                    float(line[38:46]), float(line[47:54])])
                if (label == "X") or (label == "G"):
                    # error if no symmetry type is included.
                    self.symmetrytype.append(int(line[80:85]))
                    angles = [float(i) for i in line[86:].split()]
                    if len(angles) == 0:
                        self.connectangles.append([0.])
                    else:
                        self.connectangles.append(angles)
        pdbfile.close()
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

        # add linking points to the SBU.
        # TODO(pboyd): add possible linking types
        self.add_connect_vector()
        # determine the number of H's in the SBU
        self.hydrogens = [i for i in range(len(self.atomlabel)) if 
                        self.atomlabel[i] == "H"]
    def add_connect_vector(self):
        """
        checks for certain atoms and adds topological connecting
        points.  X = the connection point, Y = the vector designed
        to orient the bond, Z = the bond vector(Z-X)
        """

        purgeatoms = []
        for indx, atom in enumerate(self.atomlabel):
            # test for "X" atoms (considered a linking point)
            if (atom == "X") or (atom == "G"):
                # there should only be a "Y" and a "Z" bonding to an "X"
                # TODO(pboyd): add error checking if Y and Z not on X
                # also, if no X atom
                if (atom == "G"):
                    # flag this bond as a special case
                    # this is the index for the bond.
                    self.special_bond.append(len(self.connectpoints))
                for i in self.bonding[indx]:
                    if self.atomlabel[i] == "Y":
                        # angle vector required to align SBUs
                        vect = vect_sub(self.coordinates[i],
                                        self.coordinates[indx])
                        self.anglevect.append(vect)
                        purgeatoms.append(i)
                    elif self.atomlabel[i] == "Z":
                        vect = vect_sub(self.coordinates[i],
                                        self.coordinates[indx])
                        self.connect_vector.append(vect)
                        purgeatoms.append(i)

                    else:
                        self.bondingatoms.append(i)
                
                pt = self.coordinates[indx]
                self.connectpoints.append(pt)
                purgeatoms.append(indx)
        # purge any coordinates in self.coordinates belonging to "X"'s
        # and "Y"'s as they are now part of self.connect_vector 
        # and self.connectpoints and self.anglevect
        self.purge(purgeatoms)

        if len(self.anglevect) == 0:
            error("anglevect is wrong")
            sys.exit(0)

    def purge(self, atoms):
        """Remove entries for X and Y atoms in the atomic coordinates"""
        atoms.sort()
        for patm in reversed(atoms):
            # correct for index of bonding atoms when purging 
            # connectivity data
            for idx, atom in enumerate(self.bondingatoms):
                if atom > patm:
                    self.bondingatoms[idx] = atom -  1
            self.coordinates.pop(patm)
            self.atomlabel.pop(patm)
            self.nbonds.pop(patm)
            self.bonding.pop(patm)

        struct = coord_matrix(self.atomlabel, self.coordinates)
        self.connect = struct[0]
        self.bonding = struct[1]
        self.nbonds = struct[2]

    def COM_shift(self):
        """
        Shifts all the coordinates such that the centre of mass (COM)
        is at the origin.
        This will make it easier for geometry shifts and rotations
        """
        for icoord, xyz in enumerate(self.coordinates):
            self.coordinates[icoord] = vect_sub(xyz, self.COM)
    def centre_of_mass(self):
        """
        calculates the centre of mass:
        sum(mass*coordinate) / sum(masses)
        
        or top / bottom

        In this calculation the bottom value is also stored
        as the self.mass
        """
        top = zeros1
        bottom = 0.
        for i in range(len(self.atomlabel)):
            massvect = scalar_mult(WEIGHT[self.atomlabel[i]], self.coordinates[i])
            top = vect_add(top, massvect)
            bottom += WEIGHT[self.atomlabel[i]]
        self.mass = bottom
        return scalar_div(bottom, top)


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
            dist = vect_length(coord[i],coord[j])
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
    # G is the metal - metal bond for linked chains
    elif("G" in (atom1,atom2) and 
        (("Y" not in (atom1,atom2))and("Z" not in (atom1,atom2)))):
        return 1.35
    elif("X" in (atom1, atom2) and 
        (("Y" not in (atom1,atom2))and("Z" not in (atom1,atom2)))):
        return 0.97
    elif("G" in (atom1,atom2) or "X" in (atom1, atom2))and \
        (("Y" in (atom1,atom2))or("Z" in (atom1,atom2))):
        return 1.1
    elif("Y" in (atom1,atom2))and("G" not in (atom1,atom2) 
            or "X" not in (atom1,atom2)):
        return 0.
    elif("Z" in (atom1,atom2))and("G" not in (atom1,atom2) 
            or "X" not in (atom1,atom2)):
        return 0.

def main():
    usage="usage: %prog 'vector string' 'vector string'"
    parser = OptionParser(usage=usage)
    (junk, file) = parser.parse_args()
    fname = file[0].strip(".pdb")
    pdbfile = File(fname)

    line = "[coordinates]\n"
    for i in range(len(pdbfile.coordinates)):
        line = line + "%5s"%(pdbfile.atomlabel[i])
        line = line + "%12.3f%8.3f%8.3f\n"%(tuple(pdbfile.coordinates[i]))
    line = line + "[connectivity]\n"
    for i in range(len(pdbfile.connectpoints)):
        line = line+ "%5i"%pdbfile.symmetrytype[i]
        line = line + "%12.3f%8.3f%8.3f%12.3f%8.3f%8.3f%12.3f%8.3f%8.3f\n"%(
                tuple(pdbfile.connectpoints[i] + pdbfile.connect_vector[i] + pdbfile.anglevect[i]))

    print line
if __name__ == "__main__":
    main()


