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
import os
import io
import itertools
import ConfigParser
from config import Options
from operations import *
from numpy import array
from elements import WEIGHT
from bookkeeping import * 
from random import random, uniform, randrange, choice
from datetime import date
from logging import warning, debug, error, info, critical
from scipy.spatial import distance
# Constants

RAD2DEG = 180./np.pi
DEG2RAD = np.pi/180.

class Generate(object):
    """
    This will contain all the wrappers for iterative and non iterative
    structure generations.
    """
   
    def __init__(self, database):
       
        # store the database locally
        self.database = database.database
        # count the number of structures generated
        self.nstructs = 0
        # the max number of mofs allowable for a particular set of
        # linkers.
        self.max_mofs = 6
        # temporary array to store up to self.max_mofs for a given
        # combination of SBUs
        self.complete_mofs = []
        # moflib is the list containing MOF structures for extensive
        # branching
        self.moflib = []
        # string history stores how MOFs were built
        self.stringhist = {} 
        # max number of sbus in a structure
        self.nsbumax = 200
        # list to store strings
        self.strings = []
        self.outdir = "output/"
        self.bondtypes = {}
        self.gen_dir()
        # list of bonding with bondtypes for successful structures
        self.bondhist = []
        # bonding rules for the bond rule generation method
        self.bondrules = []
        # counter for number of structures generated
        self.mofcount = 0

    def gen_dir(self):
        """creates directory for putting new MOF files"""
        try:
            os.mkdir(self.outdir)
        except:
            pass

    def database_generation(self):
        """Generates MOFs from a database of SBUs"""
        stopwatch = Time()
        stopwatch.timestamp()
        # generate combinations of SBUs
        build = self.set_combinations(3)
        # filter out SBU combinations of metals > 1
        build = self.filter(build)
        info("Genstruct will try %i combinations."%(len(build)))
        for idx, subset in enumerate(build):
            #check for possible MOF topologies
            buildlist = self.determine_topology(subset)
            for name in buildlist:
                self.complete_mofs = []
                basestructure = Structure(
                        [copy.deepcopy(self.database[i])
                        for i in subset],name)
                self.unique_bondtypes(basestructure)
                self.exhaustive_generation(basestructure)
                # add functional groups to built MOFs
                #if len(self.complete_mofs) > 0:
                #    info("Applying functional groups...")
                #    self.apply_functional_groups()

        stopwatch.timestamp()
        info("Genstruct finished. Timing reports %f seconds."%(stopwatch.timer))

    def apply_functional_groups(self):
        """
        Applies functional groups to built MOFs.
        """
        # determine which hydrogens to replace, what their
        # connecting atoms are (Important for overlap considerations)
        # Tally which SBU's in the MOF have hydrogens to replace
        sbus_withH = [sbu for mof in self.complete_mofs for 
                sbu in mof.sbu_array if len(sbu.hydrogens) > 0]
        dic = {}
        # isolate unique SBUs
        for sbu in sbus_withH:
            dic[sbu.index] = sbu

        sbus_withH = dic.values()
        # generate a set of hydrogen replacements for each sbu in SBUs
        # to a limit of 1000.  After this many, just quit.
        hydrogenlist = []
        for i in sbus_withH:
            #itertools chain, combinations
            s = i.hydrogens
            hydrogenlist = \
            sorted(itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1)))
        # combine 
        itertools.product(hydrogens)
        print hydrogenlist

    def determine_topology(self, indices):
        """
        Returns a list of topological names which appear across all
        SBU indices.
        """
        setlist = [set(self.database[i].connect_vector) for i in indices]
        union = set.intersection(*setlist)
        union = [st for st in union if st != "metallic"]
        return union 

    def set_combinations(self, n):
        """ Generates all combinations of a database with length n """
        indices = [i for i in range(len(self.database))]
        return list(itertools.combinations_with_replacement(indices, n))
    
    def filter(self, list):
        """ 
        Filters all combinations where there is more than one
        metal in the set.
        """
        return [j for j in list if self.metcount(j) == 1]

    def metcount(self, entry):
        """ returns the number of metals in the entry """
        mcount = 0
        for element in entry:
            if self.database[element].metal:
                mcount += 1 
        return mcount

    def sort(self, array):
        """Sort an SBU array based on indices."""
        newarray = []
        indices = [i.index for i in array]
        indices.sort()
        for i in indices:
            # MAKE SURE your SBU's have unique indices, the following
            # line will remove multiple copies of the same index.
            newarray.append([j for j in array if j.index == i][0])
        return newarray

    def apply_strings(self, dataset, indices):
        """
        (re)Builds a MOF based on a set of strings found in self.stringhist
        """
        struct = Structure(dataset)
        self.random_insert(0, dataset, struct)
        # check for periodic bonds with itself.
        struct.sbu_check(0)
        for string in self.stringhist[tuple(indices[:-1])]:
            ints = [int(i) for i in string.split("-")]
            debug("Applying %s"%(string))
            debug("Number of SBUs in structure: %i"%
                    (len(struct.mof)))
            sbu2 = len(struct.mof)
            struct.apply_string(string)
            struct.join_sbus(ints[0], ints[1], sbu2, ints[3], True)
            if struct.overlap_allcheck():
                return False
            else:
                struct.sbu_check(sbu2)
                struct.stringhist.append(string)

            struct.complete_box()
            if struct.saturated() and struct.complete_box():
                info("Structure generated!")
                struct.getfractionals()
                struct.mof_reorient()
                final_coords = struct.getfinals()
                dump = struct.coordinate_dump()
                self.nstructs += 1
                filename = "%06i_struct"%(self.nstructs)
                for i in indices:
                    filename = filename + "_%i"%(i)
                write_pdb(self.outdir+filename, dump[0], final_coords, 
                          struct.acell)
                write_xyz(self.outdir+filename, dump[0], final_coords, 
                    struct.cell, struct.origins)
                self.structcounter += 1
                # stop the timer
                return True

    def unique_bondtypes(self, structure):
        """
        Determine all the unique possible bonding which can occur
        with the dataset
        """
        self.bondtypes={}
        dataset = structure.sbu_array
        datalength = len(dataset)
        # outter loops: sbus
        bondtype = 0
        sbu_pairs = list(
                itertools.combinations_with_replacement(
                    range(datalength), 2))
        bonds = []
        for pair in sbu_pairs:
            bonds = self.possible_bonds(pair, structure)
        return bonds

    def possible_bonds(self, sbu_inds, structure):
        """ determine all pairs of bonds between two SBUs """
        sbus = [structure.sbu_array[sbu] for sbu in sbu_inds]
        name = structure.name
        bonding = []

        bonds = [[],[]]; flag_name = [[],[]]; flag_sym = [[],[]]

        standard_bonds = {}
        # add standard bonds associated with "name"
        for ind, sbu in enumerate(sbus):
            nbonds = len(sbu.connect_points[name])
            bonds[ind] += range(nbonds)
            flag_name[ind] += [name]*nbonds
            flag_sym[ind] += [i for i in sbu.connect_sym[name]]
            if name == "pcu":
                if sbu.connect_points.has_key("metallic"):
                    nbonds = len(sbu.connect_points["metallic"])
                    bonds[ind] += range(len(bonds[ind]), 
                                    len(bonds[ind])+nbonds)
                    flag_name[ind] += ["metallic"]*nbonds
                    flag_sym[ind] += [i for i in 
                            sbu.connect_sym["metallic"]]

            if sbu.metal and sbu.connect_points.has_key("intermetal"):
                nbonds = len(sbu.connect_points["intermetal"])
                bonds[ind] += range(nbonds)
                flag_name[ind] += ["intermetal"]*nbonds
                flag_sym[ind] += [i for i in 
                            sbu.connect_sym["intermetal"]]
        
        # generate list of valid bonds
        complete_set = [(i, j) for i in bonds[0] for j in bonds[1]]
        for bond in complete_set:
            # determine if this bond can be formed
            if sbus[0].metal == sbus[1].metal:
                if flag_name[0][bond[0]] == "intermetal" and \
                    flag_name[1][bond[1]] == "intermetal":
                    if flag_sym[0][bond[0]] != flag_sym[1][bond[1]]:
                        # add symmetry type
                        angle = 0
                        key = ("intermetal", sbus[0].index, 
                                flag_sym[0][bond[0]], sbus[1].index, 
                                flag_sym[1][bond[1]], angle)
                        value = ("intermetal", sbus[0].index, bond[0],
                                 sbus[1].index, bond[1])

                        self.bondtypes.setdefault(key, [])
                        if value not in self.bondtypes[key]:
                            self.bondtypes[key].append(value)
                        
            else:
                if flag_name[0][bond[0]] == flag_name[1][bond[1]]:
                    # range over all angles
                    for angle in range(len(structure.connect_angles)):
                        key = (flag_name[0][bond[0]], sbus[0].index, 
                                flag_sym[0][bond[0]], sbus[1].index, 
                                flag_sym[1][bond[1]], angle)
                        # Note: for the organic linker, the "metallic"
                        # bonds replace the regular bonds, they do not
                        # add to the total bond count as they do for
                        # the metals, so the value must be adjusted 
                        # to accomodate this.
                        bval1 = bond[0]
                        bval2 = bond[1]
                        if flag_name[0][bond[0]] == "metallic":
                            # note the reference value here assumes
                            # no difference between the number of bonds
                            # for the "metallic" and normal bonding
                            if not sbus[0].metal:
                                bval1 = bond[0] % len(
                                        sbus[0].connect_points[name])
                            elif not sbus[1].metal:
                                bval2 = bond[1] % len(
                                        sbus[1].connect_points[name])
                        value = (flag_name[0][bond[0]], sbus[0].index, 
                                 bval1, sbus[1].index, bval2, angle)
                        self.bondtypes.setdefault(key, [])
                        if value not in self.bondtypes[key]:
                            self.bondtypes[key].append(value)
                
        return bonding

    def bondtype_exists(self, arrayt):
        """
        check for existance of unique bond type
        """
        if arrayt in [i for i in self.bondtypes.itervalues()]:
            return True
        else:
            return False
            

    def generate_stringset(self, structure):
        """
        Encode the set of strings to build a MOF
        """
        self.strings = []
        count = 0
        string = "0-0-0-0-0"
        maxstring = self.setmax(structure)
        done=False
        while not done:
            if count == 1e4:
                done = True    
            self.strings.append(string)
            string = self.iterate_string(string, maxstring)
            count += 1
        return

    def random_insert(self, isbu, dataset, struct):
        """
        Randomly places a SBU anywhere
        """
        struct.mof.append(copy.deepcopy(dataset[isbu]))
        struct.connectivity.append([None]*len(dataset[isbu].connectpoints))
        # TODO(pboyd): randomly perturb the coordinates

    def exhaustive_generation(self, structure):
        """
        generates MOFs via string iteration. flag indicates what
        connectivity to call during the building process.
        """
        structcount = 0
        self.moflib = []
        self.bondhist = []
        self.generate_stringset(structure)
        self.moflib.append(structure)
        structure.sbu_check(0)
        # Start the string off at "0-0-0-0-0"
        iterstring = 0
        indexset = set([sbu.index for sbu in structure.sbu_array])
        done = False
        while not done:
            #apply the string
            for struct in self.moflib:
                string = self.strings[iterstring]
                if self.valid_string(string, struct):
                    copystruct = copy.deepcopy(struct)
                    sbu = len(copystruct.mof) 
                    copystruct.apply_string(string)
                    if copystruct.overlap_allcheck():
                        pass
                    else:
                        bondtype = self.determine_bondtype(
                                string, copystruct)
                        copystruct.sbu_check(sbu)
                        copystruct.bondhist.append(bondtype)
                        self.bondhist.append(copystruct.bondhist)
                        copystruct.stringhist.append(string)
                        self.moflib.append(copystruct)
                        
                    if copystruct.saturated() and copystruct.complete_box():
                        # check to see if all SBU's are included in 
                        # the structure.
                        base = [i.index for i in copystruct.sbu_array]
                        sets = [i[3] for i in copystruct.bondhist]
                        sets += [i[1] for i in copystruct.bondhist]
                        if len(set(sets)) < len(set(base)):
                            pass
                        else:
                            #///////////////////////////////////////
                            # TODO(pboyd):put into separate function
                            copystruct.final_coords()
                            copystruct.mof_reorient()
                            copystruct.get_cartesians()
                            self.mofcount += 1
                            pdbfile = self.outdir + "%05i_structure"%(
                                    self.mofcount)
                            for sbu in copystruct.sbu_array:
                                pdbfile += "_%i"%(sbu.index)
                            pdbfile += "_%s"%(copystruct.name)
                            copystruct.write_pdb(pdbfile)
                            # export the MOF for functional group
                            # placement
                            self.complete_mofs.append(copy.deepcopy(copystruct))
                            #///////////////////////////////////////
                            info("Structure Generated!")
                            copystruct.origins = zeros3[:]
                            #copystruct.xyz_debug()
                            structcount += 1

                    # Give up if number of structures becomes too big.
                    if len(self.moflib) > 12000:
                        done = True
                        break
                    # constraints and conditions related to specific
                    # topologies
                    if copystruct.name == "pcu" or \
                       copystruct.name == "nbo":
                        if len(indexset) > 2:
                            if structcount == 3:
                                done = True
                                break
                        else:
                            if structcount == 1:
                                done = True
                                break
                    else:
                        if structcount == 1:
                            done = True
                            break

            #iterate the string
            iterstring += 1

    def valid_struct(self, string, ind):
        """
        determine if the structure needs to be deleted.
        """
        ints = [int(i) for i in string.split("-")]
        if ints[0] >= len(self.moflib[ind].mof):
            return False
        return True

    def generate_mofstrings(self, dataset):
        """
        Generate sequences of strings to build MOFs
        """

        stringset = []
        # start by placing one in space by random_insert
        stringset.append("0")
        string = "0-0-0-0-0"
        # iterate through strings with knowledge of bonds, types, etc..
        done = False
        iter = 0
        while not done:
            iter += 1

            # keep track of which SBU# is which type
            for strings in stringset:
                # apply string if it is valid 
                if self.teststring(string, strings):
                    strings.append(string)

            # iterate string
            done = True

        return stringset

    def valid_string(self, string, struct):
        """
        determines whether a string will be valid for the next
        addition of an SBU.
        """
        name = struct.name
        idx = [int(i) for i in string.split("-")]
        if idx[0] >= len(struct.mof):
            return False
        sbu1 = struct.mof[idx[0]]
        sbu2 = struct.sbu_array[idx[2]]
        sbu_list = (sbu1, sbu2)

        # TODO(pboyd): if the topology is 'tbo' then there are max 13?
        # SBU's in the structure.  If the string goes higher than this
        # ignore.
        # Similarly, 'nbo' goes to 8.
        if struct.name == "tbo":
            if len(struct.connectivity) > 13:
                return False
        #elif struct.name == "nbo":
        #    if len(struct.connectivity) > 9:
        #        return False
        inter_metal_test = [sbu.connect_points.has_key("intermetal")
                for sbu in sbu_list]
        pcu_metal_test = [sbu.connect_points.has_key("metallic") and \
                sbu.metal and name == "pcu" for sbu in sbu_list]
        # establish some upper bounds for which the string can't surpass
        max_bonds = [len(i.connect_points[name])-1 for i in sbu_list]
        for i in range(len(sbu_list)):
            if pcu_metal_test[i]:
                max_bonds[i] += len(sbu_list[i].connect_points["metallic"])
            if inter_metal_test[i]:
                max_bonds[i] += len(sbu_list[i].connect_points["intermetal"]) 
        upper_bounds = [len(struct.mof)-1, max_bonds[0], 
                len(struct.sbu_array)-1, max_bonds[1], 
                len(struct.connect_angles)-1]

        # boundary tests
        for ind, bound in enumerate(upper_bounds):
            if idx[ind] > bound:
                return False

        # check if the bonds are already occupied
        if struct.connectivity[idx[0]][idx[1]] is not None:
            return False

        # specific bond tests

        # if the bond value is not in self.bondtype then it is not 
        # allowed.
        bond = self.determine_bondtype(string, struct)
       
        if not bond:
            return False

        # check if the bond with the same symmetry has been made
        structhist = struct.bondhist[:]
        structhist.append(bond)
        if structhist in self.bondhist:
            return False

        return True 

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
        
        return "%i-%i-%i-%i-%i"%(tuple(liststring))
 
    def setmax(self, structure):
        """
        Set maximum string
        """
        data = structure.sbu_array
        name = structure.name

        maxangle, maxbond = 0, 0
        maxsbu = len(data) - 1
        for sbu in data:
            if sbu.connect_angles.has_key(name):
                anglelen = len(sbu.connect_angles[name])-1
            else:
                anglelen = len(sbu.connect_angles["default"])-1
            numbonds = len(sbu.connect_points[name])-1
            if name == "pcu" and sbu.metal and \
                    sbu.connect_points.has_key("metallic"):
                numbonds += len(sbu.connect_points["metallic"])-1
            if sbu.metal and sbu.connect_points.has_key("intermetal"):
                numbonds += len(sbu.connect_points["intermetal"])-1
                
            maxangle = max(maxangle, anglelen)
            maxbond = max(maxbond, numbonds)
        return([self.nsbumax, maxbond, maxsbu, maxbond, maxangle])


    def determine_bondtype(self, string, struct):
        """
        Determine the type of bond formed using the string and the
        structure.
        """
        idx = [int(i) for i in string.split("-")]
        sbutype1 = struct.mof[idx[0]].index
        sbutype2 = struct.sbu_array[idx[2]].index
        bond1 = idx[1]
        bond2 = idx[3]
        angletype = idx[4]
        # determine the type of bonding taking place.  Ie what is the
        # type of bond joining the two SBUs as determined by the SBU
        # already placed.
        bondname = struct.vectorflag[idx[0]][idx[1]]
        # since index, bondtype is unique, if the bond exists it is
        # in one of the following two forms in the dictionary
        type1 = (bondname, sbutype1, bond1, sbutype2, bond2, 
                 angletype)
        type2 = (bondname, sbutype2, bond2, sbutype1, bond1, 
                 angletype)
        bondtype = [k for k, v in self.bondtypes.iteritems() 
                    if type1 in v or type2 in v]
        # if the bondtype doesn't match then most likely it's a metal - metal or 
        # organic - organic bond.
        if len(bondtype) == 0:
            return None
        return bondtype[0]

class SBU(object):
    """
    Secondary Building Unit class to keep all information
    regarding the modular units used to build MOF structures.
    """

    def __init__(self, data):
        self.vectors = []
        self.points = []

        # each connective bond is assigned a type.  If bonds are 
        # different they should have different types.
        self.name = None
        # index is used to keep track of the SBU's order when building
        # new structures.
        self.index = 0
        # internal index to keep track of the bonding rules
        self.type = None
        self.atom_label = {}
        self.coordinates = {}
        # store history of functional group replacements
        self.fnlgrpchoice = []
        self.mass = 0.
        # each bond is assigned a symmetry type.
        self.connect_sym = {}
        # vectors of connection points for SBU
        self.connect_vector = {}
        # terminal points of connecting vectors for SBU
        self.connect_points = {}
        # alignment vectors used to align SBUs once bonded together
        self.connect_align = {}
        # angles for possible SBU - SBU bond rotation
        self.connect_angles = {}
        # index of the atoms which form connecting bonds with
        # other SBUs
        self.bondingatoms = []
        self.COM = {}

        # store a connectivity table
        self.table = []
        # store indices which are Hydrogen atoms
        self.hydrogens = []
        # store indices of atoms bonded to Hydrogen atoms.
        # in the same order.
        self.hydrogens_connect = []

        self.init_arrays(data)

        self.centre_of_mass()
        self.COM_shift()
        self.find_hydrogen_atoms()

    def find_hydrogen_atoms(self):
        """
        Finds hydrogen atoms and the atom connected to it (probably C)
        """
        self.coord_matrix()
        # reference the default coordinates
        name = "default"
        numatms = len(self.atom_label[name])
        self.hydrogens = [ind for ind in range(numatms) if
                self.atom_label[name][ind] == "H"]
       

        # assumes only one atom bonded to hydrogen
        self.hydrogens_connect = [0] * len(self.hydrogens)
        for ind, hydrogen in enumerate(self.hydrogens):
            # determine bonded atom index
            for atom, dist in enumerate(self.table[hydrogen]):
                if dist > 0.:
                    atname = self.atom_label[name][atom]
                    self.hydrogens_connect[ind] = atom

    def coord_matrix(self):
        """
        Generate a coordination matrix for the SBU's coordinates.
        Currently only does 'default' coordinates.
        """
        # assume the 'default' arrangment is transferrable.
        name = "default"
        numatms = len(self.atom_label[name])
        # populate empty connectivity table
        self.table = [[0.] * numatms for i in xrange(numatms)]
        # bonding depreciated?
        bonding = [[]] * numatms

        distmatrx = distance.cdist(self.coordinates[name], 
                                   self.coordinates[name])
        fuckthis = []
        for i in range(numatms):
            for j in range(i+1, numatms):
                tol = bond_tolerance(self.atom_label[name][i], 
                                     self.atom_label[name][j])
                if (i != j) and (distmatrx[i,j] <= tol):
                    fuckthis.append((i, j, distmatrx[i,j]))
                    self.table[i][j] = distmatrx[i,j]
                    self.table[j][i] = distmatrx[j,i]
                    bonding[i].append(j)
                    bonding[j].append(i)
        return 

    def init_arrays(self, data):
        """ Converts data stream to arrays. """
        config = ConfigParser.ConfigParser()
        config.readfp(io.BytesIO(data))
        # TODO(pboyd): test for population of default sections
        # TODO(pboyd): this block should call some default values
        # instead of if/else statements.
        if config.has_section("properties"):
            if config.has_option("properties", "name"):
                self.name = config.get("properties", "name")
            else:
                self.name = None
            if config.has_option("properties", "metal"):
                self.metal = config.getboolean("properties", "metal")
            else:
                # Default False for metal
                self.metal = False

            if config.has_option("properties", "index"):
                self.index = config.getint("properties", "index")

        if config.has_section("coordinates"):
            section = "coordinates"
            for name in config.options(section):
                lines = io.BytesIO(config.get(section, name)).readlines()
                lines = [line.rstrip() for line in lines if line.rstrip()]
                # append atom label dictionary
                [self.atom_label.setdefault(name, []).append(
                    line.rstrip().split()[0]) for line in lines]
                [self.coordinates.setdefault(name, []).append(
                    [float(elem) for elem in line.rstrip().split()[1:4]]) 
                    for line in lines]
            
        if config.has_section("connectivity"):
            section = "connectivity"
            for name in config.options(section):
                lines = io.BytesIO(config.get(section, name)).readlines()
                lines = [line.rstrip() for line in lines if line.rstrip()]
                [self.connect_sym.setdefault(name, []).append(
                    int(line.rstrip().split()[0])) for line in lines]
                [self.connect_points.setdefault(name, []).append(
                    [float(elem) for elem in line.rstrip().split()[1:4]])
                    for line in lines]
                [self.connect_vector.setdefault(name, []).append(
                    [float(elem) for elem in line.rstrip().split()[4:7]])
                    for line in lines]
                [self.connect_align.setdefault(name, []).append(
                    [float(elem) for elem in line.rstrip().split()[7:10]])
                    for line in lines]

        if config.has_section("angles"):
            section = "angles"
            for name in config.options(section):
                lines = io.BytesIO(config.get(section, name)).readlines()
                angles = [float(j) for i in lines for j in i.split()]
                self.connect_angles[name] = angles
        else:
            self.connect_angles.setdefault("default", []).append(0.0)

    def choose_fnl_sites(self, fnlgrp):
        """Return an array of hydrogen sites to be replaced with
        functional groups"""
        if len(self.hydrogens) == 0:
            return []

        done = False
        trial = 0
        while not done:
            string = "%i"%(fnlgrp.index)
            sites = []
            num = randrange(len(self.hydrogens))
            for i in range(num+1):
                h = randrange(len(self.hydrogens))
                while self.hydrogens[h] in sites:
                    h = randrange(len(self.hydrogens))
                string = string + "-%i"%(self.hydrogens[h])
                sites.append(self.hydrogens[h])
            # check to make sure this arrangement of fnlgrps hasn't
            # been used yet.
            if self.fnl_group_test(string):
                self.fnlgrpchoice.append(string)
                done = True
            else:
                trial += 1
            if trial == 40:
                return None
        return sites

    def fnl_group_test(self, fnlstr):
        """Checks to see if a combination of functional group
        replacements has already been done"""
        # convert string to integers
        ints = [int(i) for i in fnlstr.split("-")]
        # iterate through old fnlgroup strings
        for old in self.fnlgrpchoice:
            oldints = [int(i) for i in old.split("-")]
            if oldints[0] == ints[0]:
                comb1 = ints[1:]
                comb2 = oldints[1:]
                comb1.sort()
                comb2.sort()
                if comb1 == comb2:
                    return False
        return True

    def add_functional_group(self, fnlgrp, sites):
        """Adds functional group to the SBU"""
        # return if no hydrogens to add to structure
        if len(self.hydrogens) == 0:
            return
        info("%i functional groups added"%(len(sites)))
        # randomly choose the H's to switch and switch them
        for i in sites:
            self.switch(i, copy.copy(fnlgrp))
        sites.sort()
        # remove the hydrogens from the SBU
        for hydrogen in reversed(sites):
            self.coordinates.pop(hydrogen) 
            self.atomlabel.pop(hydrogen)
            for idx, atom in enumerate(self.bondingatoms):
                if hydrogen < atom:
                    self.bondingatoms[idx] = atom - 1 

    def switch(self, hydrogen, fnlgrp):
        """ switch a hydrogen for a functional group"""
        if len(self.bonding[hydrogen]) != 1:
            error("Hydrogen bonded to more than one atom!")

        atom = self.bonding[hydrogen][0]
        info("adding functional group to atom %i"%(atom))
        diffvect = vect_sub(self.coordinates[hydrogen], 
                            self.coordinates[atom])
        bondvector =  scalar_mult(-1., normalize(diffvect))
        fnlgrp_vect = fnlgrp.connect_vector

        self.points.append(self.coordinates[atom])
        self.vectors.append(bondvector)

        axis = rotation_axis(fnlgrp_vect, bondvector)
        angle = calc_angle(bondvector, fnlgrp_vect)
        self.rotation(fnlgrp, axis, angle)

        # add to the SBU
        shift = vect_sub(self.coordinates[atom], fnlgrp.connectpoint)

        fnlgrp.coordinates = [vect_sub(vect_add(i,shift),
                            fnlgrp.connect_vector) for i in 
                            fnlgrp.coordinates]

        self.coordinates = self.coordinates + fnlgrp.coordinates
        for i in fnlgrp.atomlabel:
            self.atomlabel.append(i)
        self.points.append(vect_add(fnlgrp.connectpoint, shift))
        self.vectors.append(fnlgrp.connect_vector)

    def rotation(self, fnlgrp, axis, angle):
        """
        rotates a functional group to be properly added to the SBU
        """
        R = rotation_matrix(axis, angle)
        # R rotates these points about the origin, hence a correction
        # must be applied if the position of rotation is not the 
        # origin, using rotpoint * (I - R) 
        rotpoint = fnlgrp.connectpoint
        #rotpoint = np.zeros(3)
        C = matrx_mult(rotpoint, (matrx_sub(zeros3, R)))
        fnlgrp.coordinates = [vect_add(matrx_mult(i, R), C) 
                                for i in fnlgrp.coordinates]
        # the following assumes one connectpoint and one connectvector
        fnlgrp.connectpoint = vect_add(matrx_mult(fnlgrp.connectpoint, R), C)
        fnlgrp.connect_vector = matrx_mult(fnlgrp.connect_vector, R)

    def centre_of_mass(self):
        """
        calculates the centre of mass:
        sum(mass*coordinate) / sum(masses)
        
        or top / bottom

        In this calculation the bottom value is also stored
        as the self.mass
        """
        for name in self.atom_label:
            top = zeros1[:]
            bottom = 0.
            for i in range(len(self.atom_label[name])):
                massvect = \
                scalar_mult(WEIGHT[self.atom_label[name][i]], 
                        self.coordinates[name][i])
                top = vect_add(top, massvect)
                bottom += WEIGHT[self.atom_label[name][i]]
            self.mass = bottom
            self.COM[name] = scalar_div(bottom, top)
        return

    def COM_shift(self):
        """
        Shifts all the coordinates such that the centre of mass (COM)
        is at the origin.
        This will make it easier for geometry shifts and rotations
        """
        for name in self.coordinates:
            for icoord, xyz in enumerate(self.coordinates[name]):
                self.coordinates[name][icoord] = vect_sub(xyz, 
                                                self.COM[name])

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
                     "atom_vector%12.5f%12.5f%12.5f\n"
        atomformat = "%s%12.5f%12.5f%12.5f\n"
        
        line = []
        atomcount = 0
        for ibond, bondpt in enumerate(self.points):
            atomcount += 1
            line.append(bondformat%(tuple(["F"]+list(bondpt)+
                                        list(self.vectors[ibond]))))
        for icoord, coord in enumerate(self.coordinates): 
            atomcount += 1
            line.append(atomformat%(tuple([self.atomlabel[icoord]]+list(coord))))

        xyzfile.write("%i\nDebug\n"%atomcount)
        xyzfile.writelines(line)
        xyzfile.close()

class Structure(object):
    """
    Structure contains all the machinery necessary to build a MOF from
    simple organic and metallic building units
    """
 
    def __init__(self, sbu_array, name):
        # list of SBU's
        self.sbu_array = sbu_array
        self.name = name
        # TODO(pboyd): Store coordinates and connect_* values 
        # of SBUs which have been successfully added to the 
        # structure.
        self.coordinates = []
        self.atoms = []
        self.connect_points = []
        self.connect_align = []
        self.connect_vector = []
        self.connect_sym = []
        # list to keep track of bonding between SBU's
        self.connectivity = []
        self.vectorflag = []

        self.connect_angles = []
        self.store_angles()
        # keep track of all the SBUs
        self.mof = []
        # keep track of the most recent string
        self.string = ""
        # maxsbu is the max number of sbus before termination..
        self.maxsbu = 99

        self.bondhist = []
        self.fcoords = []
        # periodic stuff
        self.cell = zeros3[:]
        self.origins = zeros3[:]
        # cell parameters
        self.acell = zeros6[:]
        # cell with unit vectors
        self.ncell = zeros3[:] 
        # inverted cell
        self.icell = zeros3[:]
        # index to keep track of number of vectors added
        self.pbcindex = -1
        # history of how the structure was made
        self.stringhist = []
        self.seed(name, 0)

    def store_angles(self):
        for sbu in self.sbu_array:
            self.connect_angles += list(sbu.connect_angles.itervalues())[0]
        self.connect_angles = list(set(self.connect_angles))

    def seed(self, name, ind):
        """Seed the structure with the identity of the building topology"""
        self.mof.append(self.sbu_array[ind])
        self.import_sbu(ind, name)

    def reset(self):
        sbuarray = self.sbu_array
        self.__init__(sbuarray)

    def apply_rule(self, key, pair):
        """
        applies the bonding rule [pair] to the sbu defined by [key].
        """
        sbu = key[0]
        bond = key[1]
        angle = key[2]

        new_sbu = pair[0]
        new_bond = pair[1]
        new_angle = pair[2]

        # angles are typically placed on the organic unit, so choose
        # the angle from the organic unit.
        # Exceptions to this are when a metal cluster bonds to a 
        # metal cluster.
        if self.mof[sbu].metal:
            angle = new_angle 
        else: 
            angle = angle
        self.add_new_sbu(sbu, bond, new_sbu, new_bond, angle)

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

    def sbu_check(self, sbu, tol=None):
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
            tol = 1 
        # return dictionary of anti-parallel, unsaturated bonds within MOF
        # note, may be a problem if bonds are close but not exactly
        # anti-parallel
        bonds = self.antiparallel_bonds(sbu)
        # generate vectors for all the eligible bonds 
        bondvects = self.gen_bond_vectors(sbu, bonds)
        bondlist = [(i,tuple(j)) for i in bonds for j in bonds[i]]

        if not bonds:
            info("No new bonds formed with SBU %i"%(sbu))
            return

        for bond in bondlist:
            #info("non adj vector: (%6.3f,%6.3f,%6.3f)"%tuple(bondvects[bond]))
            adjusted_vect = self.min_img(bondvects[bond])
            #info("min img vector: (%6.3f, %6.3f, %6.3f)"%tuple(adjusted_vect))
            if length(adjusted_vect) <= tol:
                # bonds are close
                try:
                    del bonds[bond[0]]
                    info("Local bond found between "+
                        "SBU %i, bond %i, "%(sbu, bond[0])+
                        "and SBU %i, bond %i."%(bond[1][0],bond[1][1]))
                    self.join_sbus(sbu, bond[0], bond[1][0], 
                               bond[1][1], True)
                except:
                    info("Tried to bond SBU %i bond %i to SBU %i "
                         %(sbu, bond[0], bond[1][0]) +
                         "bond %i, but it has already been bonded."
                         %(bond[1][1]))
            elif self.periodic_vector(sbu, bond, bondvects[bond]):
                try:
                    del bonds[bond[0]]
                    if self.pbcindex < 2:
                        info("New periodic boundary formed between "+
                        "SBU %i, bond %i, "%(sbu, bond[0])+
                        "and SBU %i, bond %i "%(bond[1][0],bond[1][1])+
                        " with the vector (%f, %f, %f)."%(
                        tuple(bondvects[bond])))
                        self.join_sbus(sbu, bond[0], bond[1][0], 
                              bond[1][1], True)
                        self.add_vector(bondvects[bond])
                        self.origins[self.pbcindex] = \
                            self.connect_points[sbu][bond[0]][:]
                except:
                    info("Tried to bond SBU %i bond %i to SBU %i "
                         %(sbu, bond[0], bond[1][0]) +
                         "bond %i, but it has already been bonded."
                         %(bond[1][1]))

    def valid_connect(self, sbu1, bond1, sbu2, bond2):
        """
        Determine if the two SBU's can be joined together to form
        a bond
        """

        sbu1c = self.mof[sbu1]
        sbu2c = self.mof[sbu2]

        bondtype = set([self.vectorflag[sbu1][bond1], 
                       self.vectorflag[sbu2][bond2]])
        if len(bondtype) <= 1:
            if "metallic" in bondtype:
                if sbu1c.metal != sbu2c.metal:
                    return True
            elif "intermetal" in bondtype:
                if sbu1c.metal and sbu2c.metal:
                    if self.connect_sym[sbu1][bond1] != \
                            self.connect_sym[sbu2][bond2]:
                        return True
            else:
                # general case: make sure the indices are not the same
                # and that it is a metal - organic link
                test = (sbu1c.metal != sbu2c.metal) and (sbu1 != sbu2)
                if test:
                    return True
        return False

    def import_sbu(self, ind, name):
        """
        Import all relevant information from the SBU class to 
        the Structure class arrays: 
        self.coordinates
        self.atoms
        self.connect_points
        self.connect_align
        self.connect_vector
        self.connect_sym
        self.connectivity
        self.vectorflag
        """
        # check for default or metallic names
        size = len(self.connectivity)
        if name == "metallic":
            self.coordinates.append(self.sbu_array[ind].coordinates[name][:])
            self.atoms.append(self.sbu_array[ind].atom_label[name][:])
        else:
            self.coordinates.append(self.sbu_array[ind].coordinates["default"][:])
            self.atoms.append(self.sbu_array[ind].atom_label["default"][:])
        nconnect = len(self.sbu_array[ind].connect_vector[name])
        self.connectivity.append([None]*nconnect)
        self.vectorflag.append([name]*nconnect)
        self.connect_points.append(self.sbu_array[ind].connect_points[name][:])
        self.connect_vector.append(self.sbu_array[ind].connect_vector[name][:])
        self.connect_align.append(self.sbu_array[ind].connect_align[name][:])
        self.connect_sym.append(self.sbu_array[ind].connect_sym[name][:])
        # test for addition of metallic coordination bonds
        if name == "pcu" and self.sbu_array[ind].metal:
            if self.sbu_array[ind].connect_vector.has_key("metallic"):
                nconnect = len(self.sbu_array[ind].connect_vector["metallic"])
                self.connectivity[size] += [None]*nconnect
                self.vectorflag[size] += ["metallic"]*nconnect
                self.connect_points[size] += self.sbu_array[ind].connect_points["metallic"][:]
                self.connect_vector[size] += self.sbu_array[ind].connect_vector["metallic"][:]
                self.connect_align[size] += self.sbu_array[ind].connect_align["metallic"][:]
                self.connect_sym[size] += self.sbu_array[ind].connect_sym["metallic"][:]

        # test for metal - metal bonds
        if self.sbu_array[ind].metal and \
            self.sbu_array[ind].connect_vector.has_key("intermetal"):
            nconnect = len(self.sbu_array[ind].connect_vector["intermetal"])
            self.connectivity[size] += [None]*nconnect
            self.vectorflag[size] += ["intermetal"]*nconnect
            self.connect_points[size] += self.sbu_array[ind].connect_points["intermetal"]
            self.connect_vector[size] += self.sbu_array[ind].connect_vector["intermetal"]
            self.connect_align[size] += self.sbu_array[ind].connect_align["intermetal"]
            self.connect_sym[size] += self.sbu_array[ind].connect_sym["intermetal"]

    def periodic_vector(self, sbu1, bond, bondvector):
        """
        check to see if two unsaturated connecting points can be
        joined by a periodic vector.  There is a parallel test for
        the Y vector (alignment vector) as well, so more unique
        structures can be obtained.
        """
        # establish which SBUs are tyring to be joined
        bond1 = bond[0]
        sbu2 = bond[1][0]
        bond2 = bond[1][1]

        # their anti-parallel X-vectors has already been established.
        # checking now for aligned Y-vectors (alignment).
        vector1 = normalize(self.connect_align[sbu1][bond1])
        vector2 = normalize(self.connect_align[sbu2][bond2])

        # test by dot product, the absolute value should = 1
        dotprod = abs(dot(vector1, vector2))
        dottest = np.allclose(dotprod, 1., atol=0.1)

        # test by cross product, the vector should be (0,0,0)
        xprod = cross(vector1, vector2)
        xtest = np.allclose(xprod, zeros3[:], atol=0.1)
        # FIXME(pboyd): temporary fix for the fact that this constraint
        # won't work for chain-type structures like In3+
        if dottest and xtest:
            # check to see if the vector can be added to the 
            # existing lattice vectors
            if self.valid_vector(bondvector):
                return True
        return False

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

        if self.pbcindex >= 0:
            proj_a = project(bondvector, self.cell[0])
            test_a = np.allclose(length(self.cell[0]), 
                                 length(proj_a), atol=tol)
            null_a = np.allclose(proj_a, np.zeros(3), atol=tol)
            # case 1: the bond between two sbus can already be
            # established by a single boundary condition
            for i in range(self.pbcindex+1):
                if self.line_test(i, bondvector):
                    if np.allclose(length(self.cell[i]), 
                                   length(bondvector), atol=tol):
                        return True

        if self.pbcindex >= 1:
            proj_b = project(bondvector, self.cell[1])
            test_b = np.allclose(length(self.cell[1]), 
                                 length(proj_b), atol=tol)
            null_b = np.allclose(proj_b, np.zeros(3), atol=tol)
            # case 2: the bond between two sbus is at the corner of
            # two periodic boundaries.
            if self.pbcindex == 1:
                if test_a and test_b:
                    return True

        if self.pbcindex == 2:
            proj_c = project(bondvector, self.cell[2])
            test_c = np.allclose(length(self.cell[2]), 
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
                connecting_point1 = self.connect_points[sbu1][bond1[0]]
                connecting_point2 = self.connect_points[sbu2[0]][sbu2[1]]
                vect = vect_sub(connecting_point2, connecting_point1)
                # key is (bond, [sbu, bond])
                vectors[(bond1[0],tuple(sbu2))] = vect
        return vectors

    def antiparallel_bonds(self, sbu1):
        """
        Returns a dictionary of parallel bonds with each unsaturated
        bond of sbu1.  The format is {bond :[sbu#, bond#]}
        """
        # exclude bonds in sbu1 which are already bonded
        # remove sbu1 from scan
        mofscan = [i for i in range(len(self.mof))]
        bonds = [i for i in range(len(self.connectivity[sbu1])) 
                 if self.connectivity[sbu1][i] is None]
        bonding_dic = {}
        for sbu2 in mofscan:
            bondpair = [(i, j) for i in range(len(self.connectivity[sbu1]))
                               for j in range(len(self.connectivity[sbu2]))
                               if self.connectivity[sbu1][i] is None and
                               self.connectivity[sbu2][j] is None]
            for pair in bondpair:
                if self.valid_connect(sbu1, pair[0], sbu2, pair[1]):
                    vect1 = self.connect_vector[sbu1][pair[0]]
                    vect2 = self.connect_vector[sbu2][pair[1]]
                    # note the tolerance here will dictate how many
                    # structures will be found
                    if antiparallel_test(vect1, vect2, tol=0.25):
                        bonding_dic.setdefault(pair[0],[]).append([sbu2,pair[1]])
        return bonding_dic

    def overlap_allcheck(self):
        """
        checks for all pairwise atomistic overlaps in a growing MOF
        """
        # tolerance in angstroms
        tol = 0.3 
        # TODO(pboyd): include minimum image conv. when determining
        # distances.
        for idx1 in range(len(self.connectivity)-1):
            for idx2 in range(idx1+1, len(self.connectivity)):
                distmat = distance.cdist(self.coordinates[idx1], 
                                        self.coordinates[idx2])
                check = [(i,j) for i in range(len(distmat)) for 
                        j in range(len(distmat[i])) if 
                        distmat[i][j] <= tol]
                if len(check) > 0:
                    for ovlp in check:
                        warning(
                        "overlap found between SBU %i (%s), atom %i, "%
                        (idx1, self.mof[idx1].name, ovlp[0]) + "%s and SBU %i (%s),"
                        %(self.atoms[idx1][ovlp[0]], idx2, self.mof[idx2].name) +
                        "atom %i, %s."%(ovlp[1], self.atoms[idx2][ovlp[1]]))
                    return True
        return False

    def bonded(self, idx1, idx2):
        """
        determine if sbu index 1 and sbu index2 is bonded.
        """
        # only need to check to see if one of the sbu's is connected
        # to the other.
        bonds = [isbond for isbond in self.connectivity[idx1] if 
                 idx2 == isbond]
        if bonds:
            return True
        return False

    def bonded_atoms(self, idx1, idx2):
        """
        Determine which atoms form the bond between sbu's idx1
        and idx2.
        """
        sbu1 = self.mof[idx1]
        sbu2 = self.mof[idx2]
        bndindx1 = [idx for idx in range(len(self.connectivity[idx1]))
                if self.connectivity[idx1][idx] == idx2]
        bndindx2 = [idx for idx in range(len(self.connectivity[idx2]))
                if self.connectivity[idx2][idx] == idx1]
        return (sbu1.bondingatoms[bndindx1[0]], 
                sbu2.bondingatoms[bndindx2[0]])

    def overlap(self, sbu1, tol=None):
        """
        checks if one of the SBUs in structure (istruct) is
        overlapping with the other SBUs
        """
        # distance tolerance in angstroms
        if tol is not None:
            tol = tol
        else:
            tol = 1.0

        mofrange = [i for i in range(len(self.mof)) if i != sbu1]
        for sbu2 in mofrange:
            sbu2 = self.mof[sbu2]
            for atom2 in sbu2.coordinates:
                for atom1 in self.mof[sbu1].coordinates:
                    vect = vect_sub(atom1, atom2)
                    dist = self.apply_pbc(vect) 
                    if length(dist) <= tol:
                        return True
        return False

    def coordinate_dump(self):
        """ 
        Dumps all the atom labels, xyz coordinates and what SBU index
        they belong to
        """
        sbuind = []
        coords = [i for j in self.mof for i in j.coordinates]
        atoms = [i for j in self.mof for i in j.atomlabel]
        for sbu in range(len(self.mof)):
            sbuind.append([sbu for i in range(len(self.mof[sbu].coordinates))]) 
        return (atoms, coords, sbuind)

    def add_new_sbu(self, sbu1, bond1, sbutype2, bond2, iangle):
        """adds a new sbutype2 to the growing MOF"""
        # determine what type of bond is being formed
        name = self.vectorflag[sbu1][bond1]
        angle = self.connect_angles[iangle]
        # TODO(pboyd): change this from a copy of the SBU class to an
        # "import" of the relevant data.
        self.import_sbu(sbutype2, name)
        self.mof.append(copy.deepcopy(self.sbu_array[sbutype2]))
        sbu2 = len(self.mof) - 1

        info(
            "Added %s, SBU %i, bond %i to SBU %i, %s bond %i."
            %(self.mof[sbu2].name, sbu2, bond2, sbu1,
              self.mof[sbu1].name, bond1))

        self.join_sbus(sbu1, bond1, sbu2, bond2, True)
        # TODO(pboyd): add option to debug and apply if true.
        #self.xyz_debug()
        #dump = self.coordinate_dump()
        #write_xyz("history", dump[0], dump[1], self.cell, self.origins)

        # align sbu's by Z vector
        self.sbu_align(sbu1, bond1, sbu2, bond2)

        #self.xyz_debug()
        #dump = self.coordinate_dump()
        #write_xyz("history", dump[0], dump[1], self.cell, self.origins)

        # rotate by Y vector
        self.bond_align(sbu1, bond1, sbu2, bond2, angle) 

        #self.xyz_debug()
        #dump = self.coordinate_dump()
        #write_xyz("history", dump[0], dump[1], self.cell, self.origins)
    
    def bond_align(self, sbu1, bond1, sbu2, bond2, rotangle):
        """
        Align two sbus by their orientation vectors
        """

        axis = normalize(self.connect_vector[sbu2][bond2])
        angle = calc_angle(self.connect_align[sbu2][bond2],
                           self.connect_align[sbu1][bond1])

        if np.allclose(angle, 0., atol=0.01):
            return
        # origin of rotation
        rotpoint = self.connect_points[sbu2][bond2] 

        # test to see if the rotation is in the right direction
        R = rotation_matrix(axis, angle)
        test_vect = matrx_mult(self.connect_align[sbu2][bond2], R)
        test_angle = calc_angle(self.connect_align[sbu1][bond1],
                           test_vect)
        # FIXME(pboyd): the tolerance below is quite large,
        if not np.allclose(test_angle, 0., atol=1e-1):
            axis = scalar_mult(-1.,axis)
        # rotation of SBU
        self.sbu_rotation(sbu2, rotpoint, axis, angle + rotangle)

    def sbu_align(self, sbu1, bond1, sbu2, bond2):
        """
        Align two sbus, were sbu1 is held fixed and
        sbu2 is rotated and shifted to bond with sbu1
        """

        # the order of the cross product is important to 
        # signify which vector remains stationary (statvect)
        # and which will rotate.  This axis will determine the
        # proper rotation when applying the rotation matrix

        axis = rotation_axis(self.connect_vector[sbu2][bond2], 
                            scalar_mult(-1., 
                                self.connect_vector[sbu1][bond1]))
        angle = calc_angle(self.connect_vector[sbu2][bond2], 
                           scalar_mult(-1.,
                               self.connect_vector[sbu1][bond1]))
        # origin of rotation
        rotpoint = centre_of_mass(self.atoms[sbu2], self.coordinates[sbu2])

        while (np.allclose(axis, np.zeros(3), atol=1e-4)):
            randaxis = normalize([random(), random(), random()])
            randangle = uniform(0., np.pi/3.)
            self.sbu_rotation(sbu2, rotpoint, randaxis, randangle)
            axis = rotation_axis(self.connect_vector[sbu2][bond2], 
                             scalar_mult(-1.,self.connect_vector[sbu1][bond1]))
            angle = calc_angle(self.connect_vector[sbu2][bond2], 
                           scalar_mult(-1.,self.connect_vector[sbu1][bond1]))
        # rotate sbu2
        self.sbu_rotation(sbu2, rotpoint, axis, angle)

        # shift the coordinates by the shift vectors 
        shiftvector = vect_sub(self.connect_points[sbu1][bond1],
                               self.connect_points[sbu2][bond2])
        # bond points are then shifted by the shiftvector
        self.connect_points[sbu2] = [vect_add(i,shiftvector) for i 
                              in self.connect_points[sbu2]]
        self.coordinates[sbu2] = [vect_add(i,shiftvector) for i in
                            self.coordinates[sbu2]]

    def getfinals(self):
        """
        get final coordinates from fractionals.
        """
        A = self.cell[0]
        B = self.cell[1]
        C = self.cell[2]
        finals = []
        for i in self.fcoords:
            vect = zeros1[:]
            vect = vect_add(vect, scalar_mult(i[0], A))
            vect = vect_add(vect, scalar_mult(i[1], B))
            vect = vect_add(vect, scalar_mult(i[2], C))
            finals.append(vect)

        return finals


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
        C = matrx_mult(rotpoint, matrx_sub(identity3, R))
        self.coordinates[sbu] = [vect_add(matrx_mult(i, R), C) 
                            for i in self.coordinates[sbu]]
        self.connect_points[sbu] = [vect_add(matrx_mult(i, R), C)
                            for i in self.connect_points[sbu]]
        # connect_vector and anglevect are ALWAYS centered at the origin!!
        self.connect_vector[sbu] = [matrx_mult(i, R) for i in 
                                    self.connect_vector[sbu]]
        self.connect_align[sbu] = [matrx_mult(i, R) for i 
                                    in self.connect_align[sbu]]
    def get_cartesians(self):
        """
        convert fractional coordinates back to cartesian coordinates
        and store in self.coordinates.
        """
        atmcnt = 0
        for isbu in range(len(self.coordinates)):
            for icoord in range(len(self.coordinates[isbu])):
                coord = matrx_mult(self.fcoords[atmcnt], self.cell)
                self.coordinates[isbu][icoord] = coord[:]
                atmcnt += 1

    def final_coords(self):
        """
        convert all coordinates to the shifted coordinates derived
        from the periodic boundary conditions
        """
        firstsbu = self.mof[0]
        cellcentre = zeros1[:]
        for k in self.cell:
            cellcentre = vect_add(cellcentre, scalar_div(2., k))
        shiftvect = vect_sub(cellcentre, firstsbu.COM["default"])
        for sbu in range(len(self.coordinates)):
            for ibond, bond in enumerate(self.connect_points[sbu]):
                self.connect_points[sbu][ibond] = vect_add(shiftvect, bond)
            for icoord, coord in enumerate(self.coordinates[sbu]):
                self.coordinates[sbu][icoord] = vect_add(shiftvect, coord)
        self.fcoords = [self.fractional(coord) for sbu in
                        self.coordinates for coord in sbu]
        return

    def fractional(self, vector):
        
        """
        returns fractional coordinates based on the periodic boundary
        conditions
        """
        (a,b,c) = matrx_mult(vector, self.icell)
        a = a - math.floor(a)
        b = b - math.floor(b)
        c = c - math.floor(c)
        return [a,b,c]

    def mof_reorient(self):
        """
        Re-orients the cell vectors and the coordinates such that the
        first cell vector points in the x cartesian axis and the 
        second cell vector points in the xy cartesian plane
        """

        #FIXME(pboyd): change all coordinates to fractional,
        # rotate the pbc cell and the connect_vectors, anglevect
        # then re-apply the fractional coordinates.

        xaxis = [1.,0.,0.]
        yaxis = [0.,1.,0.]
        # first: rotation to the x-axis
        x_rotangle = calc_angle(self.cell[0], xaxis)
        x_rotaxis = rotation_axis(self.cell[0], xaxis)
        RX = rotation_matrix(x_rotaxis, x_rotangle)
        self.cell = [matrx_mult(i, RX) for i in self.cell]

        # second: rotation to the xy - plane
        projx_b = vect_sub(self.cell[1], project(self.cell[1], xaxis))
        xy_rotangle = calc_angle(projx_b, yaxis)
        xy_rotaxis = xaxis
        RXY = rotation_matrix(xy_rotaxis, xy_rotangle)
        # test to see if the rotation is in the right direction
        testvect = matrx_mult(projx_b, RXY)
        testangle = calc_angle(testvect, yaxis)
        if not np.allclose(testangle, 0., atol = 1e-3):
            RXY = rotation_matrix(scalar_mult(-1.,xy_rotaxis), xy_rotangle)

        self.cell = [matrx_mult(i, RXY) for i in self.cell]
        # third change: -z to +z direction
        if self.cell[2][2] < 0.:
            # invert the cell
            self.cell[2][2] = -1. * self.cell[2][2]
            self.fcoords = [vect_sub([1., 1., 1.], i) 
                            for i in self.fcoords]
        return

    def saturated(self):
        """
        returns False if at least one bond is unsaturated
        """
        for imof, mof in enumerate(self.connectivity):
            if [no_bond for no_bond in mof if no_bond is None]:
                return False
        return True

    def disjoin_sbus(self, sbu1, bond1, sbu2, bond2):
        """ deletes a bond between SBU's """
        self.connectivity[sbu1][bond1] = None
        self.connectivity[sbu2][bond2] = None


    def join_sbus(self, sbu1, bond1, sbu2, bond2, islocal):
        """ Joins two sbus together"""
        if islocal:
            self.connectivity[sbu1][bond1] = sbu2
            self.connectivity[sbu2][bond2] = sbu1
        else:
            self.connectivity[sbu1][bond1] = -1*sbu2
            self.connectivity[sbu2][bond2] = -1*sbu1
        
    def min_img(self, vector):
        """
        applies periodic boundary conditions to a vector
        """

        # FIXME(pboyd): the pbcs when the index = 1 appear broken
        A = self.cell[0]
        B = self.cell[1]
        C = self.cell[2]
        if np.allclose(vector, zeros1[:]):
            return vector
        if self.pbcindex == 0:
            # project the bond vector on the single periodic vector
            # then determine if it needs to be shifted.
            projection = project(vector, A)
            fractional = length(projection) / length(A)
            if antiparallel_test(projection, A):
                # shift vector is a multiple of the periodic boundary
                # to reduce the vector to it's minimum image
                shift_vector = scalar_mult(-1.*round(fractional), A)
            else:
                shift_vector = scalar_mult(round(fractional), A)
            # one dimensional shift is just a vector subtraction
            vector = vect_sub(vector, shift_vector)

        elif self.pbcindex == 1:
            proj_a = project(vector, A)
            proj_b = project(vector, B)
            a = length(proj_a) / length(A)
            b = length(proj_b) / length(B)
            a = round(a)
            b = round(b)
            if antiparallel_test(proj_a, self.cell[0]):
                a = -1.*a
            if antiparallel_test(proj_b, self.cell[1]):
                b = -1.*b
            vector = vect_sub(vector, scalar_mult(a,A)) 
            vector = vect_sub(vector, scalar_mult(b,B))

        elif self.pbcindex == 2:
            (a, b, c) = matrx_mult(vector, self.icell)
            a = a - round(a)
            b = b - round(b)
            c = c - round(c)

            #vector2 = matrx_mult([a,b,c], self.cell)
            vector = scalar_mult(a,A)
            vector = vect_add(vector, scalar_mult(b,B))
            vector = vect_add(vector, scalar_mult(c,C))
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

        line = []
        #atomcount = 1
        #line.append(xyzcellfmt2%(tuple(["C"] + [0.,0.,0.] +
        #            list(self.cell[0]) + list(self.cell[1])+ list(self.cell[2]))))
        atomcount = 3
        for icell, cell in enumerate(self.cell):
            line.append(xyzcellfmt1%(tuple(["C"] + list(self.origins[icell])
                        + list(cell))))
        
        for sbu in range(len(self.connect_points)):
            for ibond in range(len(self.connect_points[sbu])):
                atomcount += 1
                line.append(xyzbondfmt%(
                    tuple(["F"]+list(self.connect_points[sbu][ibond])+
                    list(self.connect_vector[sbu][ibond])+
                    list(self.connect_align[sbu][ibond]))))
            for icoord in range(len(self.coordinates[sbu])): 
                atomcount += 1
                line.append(xyzatomfmt%(tuple([self.atoms[sbu][icoord]]
                    +list(self.coordinates[sbu][icoord]))))

        xyzfile.write("%i\nDebug\n"%atomcount)
        xyzfile.writelines(line)
        xyzfile.close()

    def write_pdb(self, label):
        # re-orient the structure such that a -> x, b -> xy
        pdbfile = open('%s.pdb' % label, 'w')
        today = date.today()
        lines = []
        lines.append("%6s   %s\n"% ('REMARK',
            "GenStruct! created by: Aliens." ))
        atoms = [atom for sbu in self.atoms for atom in sbu]
        coordinates = [coord for sbu in self.coordinates for coord
                       in sbu]

        lines.append("%6s   Created: %s \n" % ('REMARK', 
                     today.strftime("%A %d %B %Y")))

        line = ["CRYST1"]
        [line.append(param) for param in self.acell]
        lines.append(pdbcellfmt%tuple(line))
        for atom in range(len(atoms)):
            line = []
            line.append("ATOM")
            # atom serial number 7 - 11
            line.append(atom+1)
            # atom name 13 - 16
            line.append(atoms[atom])
            # alternate location indicator 17
            line.append(" ")
            # residue name 18 - 20
            line.append("MOL")
            # chain ID 22
            line.append("X")
            # Residue sequence number 23 - 26
            line.append(1)
            # Insertion of residue code 27
            line.append(" ")
            # Coordinates 31 - 38, 39 - 46, 47 - 54
            [line.append(coord) for coord in coordinates[atom]]
            # Occupancy 55 - 60
            line.append(1.)
            # Temperature factor 61 - 66
            line.append(0.)
            # Element symbol (right justified) 77 - 78
            line.append(atoms[atom])
            # Charge 79 - 80
            line.append(0)
            lines.append(pdbatmfmt%(tuple(line))) 
        lines.append("TER")
        pdbfile.writelines(lines)
        pdbfile.close()


    def valid_vector(self, vect):
        """
        checks a vector against existing cell vectors
        """
        for i in range(self.pbcindex+1):
            if self.line_test(i, vect, tol=0.2):
                info("parallel boundary vector exists!")
                return False
        if self.pbcindex == 1:
            if self.plane_test(vect, tol=0.2):
                info("vector lies in the same plane!")
                return False
        if self.pbcindex == 2:
            # This may change in the future, if changeable 
            # periodic boundaries are considered.
            return False
        return True

    def add_vector(self, vect):
        self.pbcindex +=1
        self.cell[self.pbcindex] = vect
        self.ncell[self.pbcindex] = normalize(vect)
        if self.pbcindex == 2:
            self.cellparams()
            self.invert_cell()

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
        for i in range(self.pbcindex+1):
            for j in range(i+1, self.pbcindex+1):
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
            tol = 1e-1
        nvect = normalize(vector)
        ncell = normalize(self.cell[icell])
        xprod = cross(nvect,ncell)
        return np.allclose(xprod, np.zeros(3), atol=tol)

    def invert_cell(self):
        """ get the inverted cell for pbc corrections"""
        det = determinant(self.cell)

        a = self.cell[0][0]; b = self.cell[0][1]; c = self.cell[0][2]
        d = self.cell[1][0]; e = self.cell[1][1]; f = self.cell[1][2]
        g = self.cell[2][0]; h = self.cell[2][1]; k = self.cell[2][2]

        A = (e*k - f*h); B = (c*h - b*k); C = (b*f - c*e)
        D = (f*g - d*k); E = (a*k - c*g); F = (c*d - a*f)
        G = (d*h - e*g); H = (g*b - a*h); K = (a*e - b*d)

        self.icell = [[A/det, B/det, C/det],
                      [D/det, E/det, F/det],
                      [G/det, H/det, K/det]]
        return

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
            return True

class Database(object):
    """
    class to read in files and interpret them
    """
    def __init__(self, file, options):
        self.options = options
        self.database = []
        self.readfile(file)

    def readfile(self, filename):
        """populate databases with SBUs"""
        file = open(filename, "r")
        counter = 0
        body = ""
        # TODO(pboyd):Check for inappropriate inputs (no default sections, 
        # no names for subsections, duplicate names etc.)
        for line in file:
            body = body + line
            if "[start]" in line:
                counter += 1
        file.close()
        info("There are %i building units in the file: %s"
                %(counter, filename))
        #TODO(pboyd): include some index call to build a single
        #structure of choice SBUs
        body = io.BytesIO(body)
        for sbu in range(counter):
            sbuchunk = self.parselines(body)
            self.database.append(SBU(sbuchunk))

    def parselines(self, file):
        """
        Takes a chunk of the file between the lines [start] and [end]
        """
        chunk = "" 
        # ensure that the first line is [start]
        if "[start]" in file.readline():
            endread = False
            while not endread:
                line = file.readline()
                if "[end]" in line:
                    endread = True
                else:
                    chunk = chunk + line
        return chunk

def vectorshift(shift_vector, root_vector):
    """Shift a vector to within the bounds of a root vector"""

    if(np.allclose(zeros3[:], root_vector)):
        return shift_vector
    else:
        # break function if the shift_vector is zero
        if np.allclose(shift_vector, zeros3[:]):
            return shift_vector
        # project the vector to be shifted onto the root vector
        proj_vect = project(shift_vector, root_vector)

        # express the projected vector as a fraction of the root vector
        fraction = length(proj_vect) / length(root_vector)
        if antiparallel_test(proj_vect, root_vector):
            fraction = -1. * fraction
        # return the shift_vector which has been shifted such that
        # it's projection lies within the bounds of the root_vector
        return vect_sub(shift_vector, 
                scalar_mult(math.floor(fraction),root_vector))

def determinant(matrix):
    """ calculates the determinant of a 3x3 matrix"""
    a = matrix[0][0]; b = matrix[0][1]; c = matrix[0][2]
    d = matrix[1][0]; e = matrix[1][1]; f = matrix[1][2]
    g = matrix[2][0]; h = matrix[2][1]; i = matrix[2][2]
    det = a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h
    return det

def length(coord1, coord2=None):
    """ 
    Returns the length between two vectors.
    If only one vector is specified, it returns
    the length of that vector from the origin
    """
    if coord2 is not None:
        coord2 = coord2
    else:
        coord2 = zeros1[:]
    dist = distance.cdist([coord1], [coord2], 'euclidean')[0][0]
    return dist 

def rotation_axis(vect1, vect2):
    """
    axis of rotation for two vectors (90 degrees from both)
    """
    vect1 = normalize(vect1)
    vect2 = normalize(vect2)
    return normalize(cross(vect1, vect2))

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

    matrix = [[
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
        np.cos(angle) + uz * uz * (1 - np.cos(angle))]]
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
    planevector1 = normalize(vect_sub(coordinates[1],coordinates[0]))

    index = 1
    while np.allclose(sum(planevector1), 0., atol=tol):
        index += 1
        planevector1 = normalize(vect_sub(coordinates[index], 
                                           coordinates[0]))
    # test points for colinearity with the first vector
    for coord in coordinates[index+1:]:
        if(not sum(coord) == 0.):
            if(not linear_test(array([coordinates[0],coordinates[1],coord]))):
                planevector2 = normalize(vect_sub(coord, coordinates[0]))
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
            newvect = normalize(vect_sub(coord, coordinates[0]))

            # test = 0 if co-planar
            test = dot(newvect, cross(planevector1, planevector2))
            if (not np.allclose(0., test, atol=tol)):
                return False
    
    return True
        

def normalize(vector):
    """changes vector length to unity"""
    if sum(inner(vector,vector)) == 0.:
        return vector
    return scalar_div(math.sqrt(dot(vector, vector)), vector)

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
    vector1 = normalize(vect_sub(coordinates[1], coordinates[0]))

    for point in coordinates[2:]:
        if sum(point) != 0.:
            vector2 = normalize(vect_sub(point, coordinates[0]))
            crossprod = cross(vector2,vector1)
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
    return scalar_mult(length(vector1) * np.cos(angle), unit2)

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
        tol = 0.2 
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

def centre_of_mass(atoms, coordinates):
    """
    calculates the centre of mass:
    sum(mass*coordinate) / sum(masses)
    
    or top / bottom

    In this calculation the bottom value is also stored
    as the self.mass
    """
    top = zeros1[:]
    bottom = 0.
    for i in range(len(atoms)):
        massvect = scalar_mult(WEIGHT[atoms[i]], coordinates[i])
        top = vect_add(top, massvect)
        bottom += WEIGHT[atoms[i]]
    return scalar_div(bottom, top)

def main():
    """Default if run as an executable"""

    # initialize logging 
    Log()
    open('history.xyz', 'w')
    open('debug.xyz', 'w')
    options = Options()
    sbutest = Database("wilmerdatabase", options)
    genstruct = Generate(sbutest)
    genstruct.database_generation()

if __name__ == '__main__':
    main()
