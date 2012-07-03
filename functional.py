#!/usr/bin/env python
import sys
import io
from ConfigParser import ConfigParser
from scipy.spatial import distance
from operations import *

class Functional_groups(object):
    """
    Functional groups contain connectivity info, file reading etc.
    """
    def __init__(self):
        """
        Read in library of functional groups, build arrays.
        """
        #lib_file = "functional_groups.lib"
        lib_file = "group.lib"
        self.groups = ConfigParser()
        self.groups.read(lib_file)
        self.name = {}
        self.atoms = {}
        self.coordinates = {}
        self.connect_vector = {}
        self.connect_align = {}
        self.bond_length = {}
        self.connect_points = {}
        self.table = {}
        self.populate_arrays()

    def populate_arrays(self):
        """
        Populate the self.atoms, self.coordinates, self.connect_vectors
        arrays.
        """
        fnl_groups = self.groups.sections()
        for group in fnl_groups:
            lines = io.BytesIO(self.groups.get(group, "atoms").strip("\n")).readlines()
            tmpatoms = []
            tmpcoord = []
            for line in lines:
                line = line.strip("\n")
                tmpatoms.append(line.split()[0])
                tmpcoord.append([float(i) for i in line.split()[1:4]])
            idx = self.groups.getint(group,"index")
            self.atoms[idx] = tmpatoms
            self.coordinates[idx] = tmpcoord
            self.connect_vector[idx] = [float(i) for i in 
                    self.groups.get(group, "orientation").split()]
            self.bond_length[idx] = self.groups.getfloat(group, "carbon_bond")
            self.connect_align[idx] = [float(i) for i in 
                    self.groups.get(group, "normal").split()]
            self.connect_points[idx] = self.groups.getint(group, "connection_point")
            self.name[idx] = self.groups.get(group, "name")
            self.coord_matrix(idx)
    def coord_matrix(self, idx):
        """
        Generate a coordination matrix for the SBU's coordinates.
        """
        numatms = len(self.atoms[idx])
        # populate empty connectivity table
        self.table[idx] = [[0.] * numatms for i in xrange(numatms)]
        distmatrx = distance.cdist(self.coordinates[idx], 
                                   self.coordinates[idx])
        for i in range(numatms):
            for j in range(i+1, numatms):
                tol = bond_tolerance(self.atoms[idx][i], 
                                     self.atoms[idx][j])
                if (i != j) and (distmatrx[i,j] <= tol):
                    self.table[idx][i][j] = distmatrx[i,j]
                    self.table[idx][j][i] = distmatrx[j,i]
        return 

        
def main():
    fnlgrp = Functional_groups()

if __name__ == "__main__":
    main()

