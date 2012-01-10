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
import bookkeeping
from random import random, uniform, randrange, choice
from datetime import date

# Constants

RAD2DEG = 180./np.pi

class Generate(object):
    """
    This will contain all the wrappers for iterative and non iterative
    structure generations.
    """
   
    def __init__(self, choice):
        
        self.bondrules = []
        self.anglrules = []
        self.generation()

    def generation(self):
        """
        Provides rules for generating a MOF from a subset of SBUs,
        these rules will be altered (randomly??? or iteratively???)
        a few times regardless of build success. 
        """
        dataset = [SBU("Copper Paddlewheel", xsdfile="Cu_paddlewheel"),
                   SBU("BTC", xsdfile="btc_stripped", ismetal=False),
                   SBU("BDC", xsdfile="bdc_stripped", ismetal=False),
                   SBU("linear_Fe", xsdfile="linear_Fe")]

        done = False
        self.set_rules(dataset)
        count = 0
        while not done:
            # call the Structure class with a set of SBU's and rules to
            # attempt to build the MOF
            struct = Structure(dataset, self.bondrules, self.anglrules)
            struct.build()
            if not struct.goodstruct:
                self.change_rules(dataset)
            else:
                done = True
            count += 1
            if count == 5:
                done = True

    def change_rules(self, dataset):
        """
        Big random change of the bonding and angle rules for building 
        a MOF.  
        """
        # TODO(pboyd): use different keys to determine the type of change
        # ie. bigrandom, smallrandom, bigiterative, smalliterative etc..

        # reset the bond rules
        self.bondrules = []
        self.anglrules = []
        organic = [i for i in range(len(dataset)) if not dataset[i].ismetal]
        org_count = len(organic)
        metal = [i for i in range(len(dataset)) if dataset[i].ismetal]
        met_count = len(metal)
        for number, unit in enumerate(dataset):
            if unit.ismetal:
                numbonds = len(unit.connectpoints)
                bondstr = "" 
                for i in range(numbonds):
                    # Random bond choice
                    bondstr = bondstr + str(organic[randrange(org_count)])
                # anglestr typically ignored for the metals, this will
                # be applied to the organic linkers
                anglestr = ""
                # angle choice random as is the initial.
                for angles in unit.connectangles:
                    numangles = len(angles)
                    anglestr = anglestr + str(randrange(numangles))
            else:
                numbonds = len(unit.connectpoints)
                numangle = len(unit.connectangles)
                bondstr = "" 
                for i in range(numbonds):
                    # Random bond choice
                    bondstr = bondstr + str(metal[randrange(met_count)])
                anglestr = ""
                # angle choice random as is the initial.
                for angles in unit.connectangles:
                    numangles = len(angles)
                    anglestr = anglestr + str(randrange(numangles))

            self.bondrules.append(bondstr)
            self.anglrules.append(anglestr)

        return

    def set_rules(self, dataset):
        """
        initial generation of building rules for new MOF
        """
        # currently make bonding only between metal and organic SBUs
        organic = [i for i in range(len(dataset)) if not dataset[i].ismetal]
        org_count = len(organic)
        metal = [i for i in range(len(dataset)) if dataset[i].ismetal]
        met_count = len(metal)
        for number, unit in enumerate(dataset):
            if unit.ismetal:
                numbonds = len(unit.connectpoints)
                # alternate between SBUs for initial config
                # FIXME(pboyd): note, problems if the number of SBUs
                # is greater than the number of bonds on the metal
                bondstr = "%i" * numbonds %(tuple(
                            [organic[(i)%org_count] for i in range(numbonds)]))
                # anglestr typically ignored for the metals, this will
                # be applied to the organic linkers
                anglestr = ""
                for angles in unit.connectangles:
                    numangles = len(angles)
                    # Random? choosing of initial angles
                    anglestr = anglestr + str(randrange(numangles))
            else:
                numbonds = len(unit.connectpoints)
                numangle = len(unit.connectangles)
                # alternate between SBUs for initial config
                bondstr = "%i" * numbonds %(tuple(
                            [metal[(i)%met_count] for i in range(numbonds)]))
                anglestr = ""
                for angles in unit.connectangles:
                    numangles = len(angles)
                    # Random? choosing of initial angles
                    anglestr = anglestr + str(randrange(numangles))

            self.bondrules.append(bondstr)
            self.anglrules.append(anglestr)

        return



class SBU(object):
    """
    """

    def __init__(self, type=None, atomlabel=None, coordinates=None,
                 pdbfile=None, xsdfile=None, ismetal=None):
        #self.log = bookkeeping.Log()

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

    def plane_test(self, vector):
        """
        Check to see if the vector is planar with the other
        cell vectors.
        """
        return planar_test(array([np.zeros(3)]+[x for x in self.cell]+[vector]))

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
 
    def __init__(self, sbu_array, bondrules, anglrules):
        self.pbc = Cell() 
        # the list mof keeps track of all the SBUs
        self.mof = []
        # initialize logging for this MOF
        self.log = bookkeeping.Log()
        # iterators for the sbu and it's connection point 
        self.sbuind, self.bondind = 0, 0
        # list to keep track of bonding between SBU's
        self.connectivity = []
        # list of SBU's
        self.sbu_array = sbu_array
        for i, j in enumerate(self.sbu_array):
            j.index = i
        self.bondrules = bondrules
        self.anglrules = anglrules

        self.goodstruct = False
        self.terminate = False

        self.sbuind = 0
        self.bondind = 0

        print self.bondrules, self.anglrules

    def build(self):
        """Exhaustively builds MOF with metals and linkers"""
    
        self.mof.append(copy.deepcopy(self.sbu_array[0]))
        self.connectivity.append([None]*len(self.mof[0].connectpoints))
        self.xyz_debug()
        dump = coordinate_dump(self.mof)
        write_xyz("history", dump[1], dump[0], self.pbc.cell)
   
        # progressive building of MOF, will complete when all
        # bonds are saturated
        while not self.terminate:
            if self.pbc.index == 2:
                self.pbc.complete_box()
            # check full saturation condition
            if self.saturated():
                self.mof_reorient()
                self.goodstruct = True
                self.terminate = True
            # try growing the mof, if fails then go back and change
            # bond rules
            self.add_new_sbu(self.sbuind, self.bondind)
            newsbu = len(self.mof) - 1
            self.sbu_check(newsbu)
            self.increment()

        if self.goodstruct:
            # Apply pbc's to all coordinates
            self.final_coords()

        # dump framework into an xyz file
        dump = coordinate_dump(self.mof)
        write_xyz("FINAL", dump[1], dump[0], self.pbc.cell, self.pbc.origins)
        write_pdb("FINAL", dump[1], dump[0], self.pbc.acell)

    def add_new_sbu(self, exind, exbond):
        """adds a new sbu to the growing MOF"""
        if self.connectivity[exind][exbond] is None:

            # choose an SBU to attach according to the rules of
            # the SBU exind
            sbutype = self.choose_SBU(exind, exbond)
            # following are lists of possible bonds and angles.
            sbubond = self.choose_bond(exind, sbutype)
            sbuangle = self.choose_angle(exind, exbond, sbutype, sbubond)
            # attempt to attach the selected sbu to the designated bond,
            # if this is not achieved, return with negative to the 
            # outter process.
            self.mof.append(copy.deepcopy(self.sbu_array[sbutype]))
            newind = len(self.mof) - 1
            self.connectivity.append([None]*len(self.mof[newind].connectpoints))
            iter = 0
            done = False
            while not done:
                # try to attach the bond and angles listed
                newbond = sbubond[iter]
                angle = sbuangle[iter] 
                self.log.info(
                    "added %s, sbu %i, bond %i to SBU %i, %s bond %i"
                    %(self.mof[newind].type, newind, newbond, exind,
                      self.mof[exind].type, exbond))
                self.connectivity[exind][exbond] = newind
                self.connectivity[newind][newbond] = exind

                self.xyz_debug()
                dump = coordinate_dump(self.mof)
                write_xyz("history", dump[1], dump[0], self.pbc.cell, self.pbc.origins)

                self.sbu_align(exind, exbond, newind, newbond)

                self.xyz_debug()
                dump = coordinate_dump(self.mof)
                write_xyz("history", dump[1], dump[0], self.pbc.cell, self.pbc.origins)
                # rotate by Y vector
                self.bond_align(exind, exbond, newind, newbond, angle) 

                self.xyz_debug()
                dump = coordinate_dump(self.mof)
                write_xyz("history", dump[1], dump[0], self.pbc.cell, self.pbc.origins)
                if self.overlap(newind):
                    self.log.info("Overlap found with SBU %i"%newind)
                    self.connectivity[exind][exbond] = None 
                    self.connectivity[newind][newbond] = None 
                    iter+=1
                else:
                    done = True
                if iter == len(sbubond):
                    # no structure made, raise flag
                    self.log.info("No structure made according to "+\
                            "bonding rules.")
                    self.goodstruct = False
                    self.terminate = True
                    done = True

    def choose_SBU(self, exind, exbond):
        """
        chooses the SBU and bond to add according to the rules
        designated for SBU exind
        """
        sbu1 = self.mof[exind].index 
        sbu2 = int(self.bondrules[sbu1][exbond])
        return sbu2

    def choose_bond(self, exind, sbutype):
        """
        returns a list of possible bonds to bond to the sbu exind
        """
        sbu1 = self.mof[exind].index
        bondrule = self.bondrules[sbutype]
        return [i for i in range(len(bondrule)) if int(bondrule[i]) == sbu1]

    def choose_angle(self, exind, exbond, sbutype, sbubond):
        """
        returns a list of angles associated with the bonds possible
        for binding of SBU exind with SBU sbutype
        """
        # angles will typically be listed with the organic SBU
        # hence metal SBUs will need to find the accompanying 
        # organic bonding angle

        extype = self.mof[exind].index
        bondrule = self.bondrules[sbutype]
        if self.mof[exind].ismetal:
            # adding organic linker, take angles from it's
            # own anglerules
            anglrule = self.anglrules[sbutype]
            array = [self.sbu_array[sbutype].connectangles[i][
                    int(anglrule[i])] for i in sbubond]
        else:
            # adding metallic cluster, take angle from
            # the organic unit
            anglrule = str(self.anglrules[extype][exbond])
            array = [self.mof[exind].connectangles[exbond][
                     int(anglrule)] for i in range(len(sbubond))]
        return array 
    
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

    def increment(self):
        """
        Increments the sbu and/or the bond index based on the 
        values given
        """
        if self.bondind == len(self.mof[self.sbuind].connectpoints) - 1:
            self.sbuind += 1
            self.bondind = 0
        else:
            self.bondind += 1

        return
    def saturated(self):
        """
        returns False if at least one bond is unsaturated
        """
        for imof, mof in enumerate(self.connectivity):
            if [no_bond for no_bond in mof if no_bond is None]:
                return False
        return True

    def sbu_check(self, sbu_ind1):
        """
        Checks an SBU for 
        1. bonding locally,
        2. bonding by existing pbc
        """
        # exclude bonds in scan which are already bonded in the sbu
        bonds = [i for i in range(len(self.connectivity[sbu_ind1])) 
                 if self.connectivity[sbu_ind1][i] is None]
        # exclude the sbu in question when scanning the other sbus
        # (it generally will not bond to itself!)
        mofscan = [i for i in range(len(self.mof)) if i is not sbu_ind1]

        # loop over bonds in sbu1, inner loop over all other
        # sbus and their bonds
        for sbu_ind2 in mofscan:
            for nbond2 in range(len(self.connectivity[sbu_ind2])):
                for nbond1 in bonds:
                    # TODO(pboyd): investigate if an inverted matrix with a 
                    # zero vector will perform the same as the projection 
                    # method below.
                    connecting_point1 = self.apply_pbc(
                            self.mof[sbu_ind1].connectpoints[nbond1])
                    connecting_point2 = self.apply_pbc(
                            self.mof[sbu_ind2].connectpoints[nbond2])
                    if self.is_valid_bond(sbu_ind1, nbond1, sbu_ind2, 
                                          nbond2, connecting_point1,
                                          connecting_point2):
                        print "valid bond reached"
                        # local bond, no periodic boundaries considered.
                        if self.local_attach_sbus(
                                self.mof[sbu_ind1].connectpoints[nbond1],
                                self.mof[sbu_ind2].connectpoints[nbond2]): 
                            print "local bond found between sbu %i bond %i, and sbu %i bond %i"%(
                                sbu_ind1, nbond1, sbu_ind2, nbond2)
                            self.add_bond(sbu_ind1, nbond1, sbu_ind2,
                                          nbond2, True)
                            return 
                        # remote bonding via application of periodic 
                        # boundaries.
                        if self.remote_attach_sbus(
                                self.mof[sbu_ind1].connect_vector[nbond1],
                                connecting_point1, 
                                self.mof[sbu_ind2].connect_vector[nbond2],
                                connecting_point2):
                            print "remote bond found between sbu %i bond %i, and sbu %i bond %i"%(
                                sbu_ind1, nbond1, sbu_ind2, nbond2)
                            self.add_bond(sbu_ind1, nbond1, sbu_ind2, 
                                          nbond2, False)
                            return
                        # remote bond and addition of new boundary
                        # condition.
                        if self.pbc.index < 2 and self.new_pbc(
                                self.mof[sbu_ind1].connect_vector[nbond1],
                                self.mof[sbu_ind1].connectpoints[nbond1], 
                                self.mof[sbu_ind2].connect_vector[nbond2],
                                self.mof[sbu_ind2].connectpoints[nbond2]):
                            self.add_bond(sbu_ind1, nbond1, sbu_ind2, 
                                          nbond2, False)
                            bond_vector = self.mof[sbu_ind1].connectpoints[nbond1] -\
                                          self.mof[sbu_ind2].connectpoints[nbond2]
                            print "new cell found between sbu %i bond %i, and sbu %i bond %i"%(
                                sbu_ind1, nbond1, sbu_ind2, nbond2), bond_vector
                            self.pbc.add_vector(bond_vector)
                            self.pbc.origins[self.pbc.index] = \
                                    self.mof[sbu_ind2].connectpoints[nbond2][:]
                            # When a new pbc is added, there should be
                            # a scan for other possible bonds with 
                            # existing sbus
                            return

                        print "didn't do anything!", connecting_point2 - connecting_point1

    def add_bond(self, sbu1, bond1, sbu2, bond2, islocal):
        """ Joins two sbus together"""
        if islocal:
            self.connectivity[sbu1][bond1] = sbu2
            self.connectivity[sbu2][bond2] = sbu1
        else:
            self.connectivity[sbu1][bond1] = -1*sbu2
            self.connectivity[sbu2][bond2] = -1*sbu1

    def is_valid_bond(self, sbu1, bond1, sbu2, bond2, 
                      connecting_point1, connecting_point2):
        """
        Checks for anti-parallel bonding vectors between two SBUs,
        makes sure they satisfy bonding requirements i.e.
        they are not bonded already and they are linking to the 
        appropriate SBU
        """
        if (self.connectivity[sbu2][bond2] is not None):
            return False

        if self.mof[sbu1].ismetal == self.mof[sbu2].ismetal:
            return False

        # FIXME(pboyd): this returns a valid bond if two connecting
        # points are close to one another, this overlooks any anti
        # parallel vector considerations or aligned anglevects...

        if  points_close(connecting_point1, 
                         connecting_point2, 1.5):
            return True

        return anti_parallel_test(self.mof[sbu1].connect_vector[bond1],
                                  self.mof[sbu2].connect_vector[bond2])
        
    def local_attach_sbus(self, connecting_point1, connecting_point2):
        """
        Checks to see if the bond is local and apply boundary 
        corrections if pbc vectors are present
        """
        # TODO(pboyd): additional checks, parallel connect_vectors,
        # parallel anglevects
        test = points_close(connecting_point1, connecting_point2, 1.5)
        return test 

    def remote_attach_sbus(self, connect_vector1, connecting_point1, 
                           connect_vector2, connecting_point2):
        """
        Join two SBUs if they have anti-parallel bonding vectors,
        are pointing away from each other, and existing pbc vectors
        can join them.
        """
        # tolerance for matching vector lengths.
        tol = 0.5
        # if vectors are pointing away, they can be
        # joined by a periodic boundary
        #if self.pointing_away(connect_vector1, connecting_point1, 
        #                      connect_vector2, connecting_point2):

        # check if a local bond can be applied after applying pbc
        # conditions.  NOTE: this is still considered a "remote"
        # attachment since pbc conditions are used.
        if self.local_attach_sbus(connecting_point1, connecting_point2):
            return True

        bond_vector = connecting_point1 - connecting_point2

        if self.pbc.index >= 0:
            proj_a = project(bond_vector, self.pbc.cell[0])
            test_a = np.allclose(length(self.pbc.cell[0]), 
                                 length(proj_a), atol=tol)
            null_a = np.allclose(proj_a, np.zeros(3), atol=tol)
            # case 1: the bond between two sbus can already be
            # established by a single boundary condition
            for i in range(self.pbc.index+1):
                if self.pbc.line_test(i, bond_vector):
                    if np.allclose(length(self.pbc.cell[i]), 
                                   length(bond_vector), atol=tol):
                        return True

        if self.pbc.index >= 1:
            proj_b = project(bond_vector, self.pbc.cell[1])
            test_b = np.allclose(length(self.pbc.cell[1]), 
                                 length(proj_b), atol=tol)
            null_b = np.allclose(proj_b, np.zeros(3), atol=tol)
            # case 2: the bond between two sbus is at the corner of
            # two periodic boundaries.
            if self.pbc.index == 1:
                if test_a and test_b:
                    return True

        if self.pbc.index == 2:
            proj_c = project(bond_vector, self.pbc.cell[2])
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

    def new_pbc(self, connect_vector1, connecting_point1,
                connect_vector2, connecting_point2):
        """
        Determine if the sbus can form a new boundary vector 
        """
        bond_vector = connecting_point1 - connecting_point2 
        if self.pbc.valid_vector(bond_vector):
            if self.pbc.index < 2:
                return True

        return False

    def pbc_scan(self):
        """
        Scans the growing MOF for new attachments with pbc vectors.
        If the result is over-lap then the vector is bad.
        """
        # TODO(pboyd): Do not permanently shift vectors by the periodic
        # boundaries, but when scannning new SBUs, implement the pbc 
        # correction to the scan.
        boundaries = [self.pbc.cell[x] for x in range(self.pbc.index+1)]
        boundaries = array(boundaries)
        tempcoords, tempconnect = [], []
        for i, j in enumerate(boundaries):
            for k, m in enumerate(self.mof):
                for n, o in enumerate(m.coordinates):
                    tempcoords.append(vectorshift(o,j))
                for p, q in enumerate(m.connectpoints):
                    tempconnect.append(vectorshift(q,j))

            #for p in range(len(tempcoords)/2):
            #    for q in range(len(tempcoords)/2, len(tempcoords)):
            #        if points_close(tempcoords[p],tempcoords[q]):
            #            self.log.info("Overlap found with pbc's")
                        #self.pbc_break(i)
            #            return

    def apply_pbc(self, vector):
        """
        applies periodic boundary conditions to a vector
        """
        A = self.pbc.cell[0]
        B = self.pbc.cell[1]
        C = self.pbc.cell[2]
        if np.allclose(vector, np.zeros(3)):
            return vector
        if self.pbc.index == 0:
            vector = vectorshift(vector, A)
            #vector = vector
        elif self.pbc.index == 1:
            ab = dot(A, B)
            proj_a = project(vector, self.pbc.cell[0])
            proj_b = project(vector, self.pbc.cell[1])
            a = length(proj_a) / length(self.pbc.cell[0])
            b = length(proj_b) / length(self.pbc.cell[1])
            if anti_parallel_test(proj_a, self.pbc.cell[0]):
                a = -1*a
            if anti_parallel_test(proj_b, self.pbc.cell[1]):
                b = -1*b
            #vector = vector
            vector = vector - (math.floor(a)*A) - (math.floor(b)*B)
        elif self.pbc.index == 2:
            bca = dot(B, cross(C,A))
            a = dot(C, cross(vector, B)) / bca
            b = dot(C, cross(vector, A)) / (-1.*bca)
            c = dot(B, cross(vector, A)) / bca
            a = a - math.floor(a)
            b = b - math.floor(b)
            c = c - math.floor(c)
            vector = (a)*A + (b)*B + (c)*C
            #vector = vector
        return vector

    def pbc_break(self, ind):
        """
        Removes a periodic boundary, removes 
        bonds from SBUs, resets MOF index
        """
        # FIXME(pboyd): Not done yet.
        for i in self.connectivity:
            self.connectivity[i[1]] = None
        self.bondind = 0
        self.sbuind = 0

    def pointing_away(self, connect_vector1, connecting_point1,
                      connect_vector2, connecting_point2):
        """
        Checks if two vectors are pointing away from each other
        """
        tail1 = connecting_point1 -\
                connect_vector1 
        head1 = connecting_point1

        tail2 = connecting_point2 -\
                connect_vector2 
        head2 = connecting_point2 

        if (length(tail1, tail2) < length(head1, head2)):
            return True
        
        return False
    
    def overlap(self, struct):
        """
        checks if one of the SBUs in structure (istruct) is
        overlapping with the other SBUs
        """
        # distance tolerance in angstroms
        disttol = .5
        exclrange = range(len(self.mof))
        exclrange.remove(struct)
        overlapping = False
        for sbu in exclrange:
            for coords in self.mof[sbu].coordinates:
                coords = self.apply_pbc(coords)
                for xyz in self.mof[struct].coordinates:
                    xyz = self.apply_pbc(xyz)
                    if (length(xyz, coords) <= disttol):
                        print "overlap found with SBU %i and SBU %i"%(
                                struct, sbu)
                        overlapping = True

        return overlapping

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

    def bond_align(self, sbu1, bond1, sbu2, bond2, rotangle):
        """
        Align two sbus by their orientation vectors
        """

        sbu1 = self.mof[sbu1]
        sbu2 = self.mof[sbu2]

        axis = normalize(sbu2.connect_vector[bond2])
        angle = calc_angle(sbu1.anglevect[bond1],
                           sbu2.anglevect[bond2])

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

        # Change sbu1 and sbu2 from indices to reference to the
        # class
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
        if anti_parallel_test(proj_vect, root_vector):
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
    return math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + 
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
    dot12 = vect1[0]*vect2[0] + vect1[1]*vect2[1] + vect1[2]*vect2[2]
    dist11 = math.sqrt(vect1[0]*vect1[0] + vect1[1]*vect1[1] +\
             vect1[2]*vect1[2])
    dist22 = math.sqrt(vect2[0]*vect2[0] + vect2[1]*vect2[1] +\
             vect2[2]*vect2[2])
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

    #TODO(pboyd): the above matrix is the transpose of the rotation
    # matrix.  This has caused me a lot of grief.  Fix it
    # permanently.
    matrix = np.transpose(matrix)
    return matrix

def planar_test(coordinates, tol=None):
    """
    Tests whether four or more points are within a plane.
    Returns False if only one of the points is co-planar
    """

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
            test = np.inner(newvect, np.cross(planevector1, planevector2))
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

def anti_parallel_test(vector1, vector2, tol=None):
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
        tol = 1e-2 
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
    v0 = vector1[0] * vector2[0]
    v1 = vector1[1] * vector2[1]
    v2 = vector1[2] * vector2[2]
    return (v0 + v1 + v2) 

def cross(vector1, vector2):
    """
    returns the cross product of two three dimensional vectors
    """
    v0 = (vector1[1] * vector2[2]) - (vector1[2] * vector2[1])
    v1 = -1.*((vector1[0] * vector2[2]) - (vector1[2] * vector2[0]))
    v2 = (vector1[0] * vector2[1]) - (vector1[1] * vector2[0])
    return array([v0,v1,v2])

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

def select_dataset():
    """
    selects a subset of SBUs from the database and passes it on for 
    structure generation
    """
    # TODO(pboyd): set up this part
    return

def main():
    """Default if run as an executable"""

    choice = "iterative"
    genstruct = Generate(choice)

#    struct = Structure()
#    struct.build()
if __name__ == '__main__':
    main()
