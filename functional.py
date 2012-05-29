#!/usr/bin/env python
import sys
import io
from ConfigParser import ConfigParser

class Functional_groups(object):
    """
    Functional groups contain connectivity info, file reading etc.
    """
    def __init__(self):
        """
        Read in library of functional groups, build arrays.
        """
        lib_file = "functional_groups.lib"
        self.groups = ConfigParser()
        self.groups.read(lib_file)
        self.name = {}
        self.atoms = {}
        self.coordinates = {}
        self.connect_vector = {}
        self.connect_align = {}
        self.bond_length = {}
        self.connect_points = {}
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
def main():
    fnlgrp = Functional_groups()

if __name__ == "__main__":
    main()

