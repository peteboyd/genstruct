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
import itertools
from operations import *
from numpy import array
from elements import WEIGHT, ATOMIC_NUMBER
import sample
from bookkeeping import Log, Time
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
   
    def __init__(self):
        
        # count the number of structures generated
        self.nstructs = 0
        # count the number of structures generated with a particular
        # set of linkers
        self.structcounter = 0
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
    def gen_dir(self):
        """creates directory for putting new MOF files"""
        try:
            os.mkdir(self.outdir)
        except:
            pass

    def bondrule_generation(self):
        """ Generates MOFs based on a set of bonding rules for each SBU"""

        database = Database()
        database.build_databases()
        stopwatch = Time()
        stopwatch.timestamp()

        metals = self.sort([SBU(pdbfile=i) for i in database.metals])
        orgncs = self.sort([SBU(pdbfile=i, ismetal=False) 
                            for i in database.organic])
        fnlgrp = [None] + self.sort([Functional_group(file=i) 
                                    for i in database.fnlgrps])

        for i in metals:
            for j in orgncs:
                for m in fnlgrp:
                    # set up the bond rules for these SBUs
                    self.set_bond_rules([i, j])
                    done = False
                    idx = 0
                    while not done:
                        dataset = [copy.deepcopy(sbu) for sbu in [i, j]]
                        for idx, sbu in enumerate(dataset):
                            sbu.type = idx
                        rule = self.bondrules[idx]
                        while not self.apply_bondrule(dataset, rule):
                            info("Incrementing bond rule")
                            idx+=1
                            rule = self.bondrules[idx]
                        done = True

        stopwatch.timestamp()
        info("Finished. Time lapse : %f seconds"%stopwatch.timer)

    def set_bond_rules(self, dataset):
        """
        Populates a dictionary of bond rules for reference when
        building a MOF
        """
        # generate SBU, bond, angle data
        bonddata = [(i, j, k) for i in range(len(dataset)) for j in 
                    range(len(dataset[i].connectpoints)) for k in 
                    range(len(dataset[i].connectangles[j]))]
        # generate all pairwise bonds including self - self bonds
        bondpairs = list(itertools.combinations(bonddata, 2))
        bondpairs += list(itertools.izip(bonddata,bonddata))

        # remove forbidden bonds
        bondpairs = [pairs for pairs in bondpairs if 
                     self.bonding_allowed(pairs, dataset)]
        alternates = []
        for pair in bondpairs:
            alternates.append((pair[1], pair[0]))
        bondpairs = bondpairs + alternates
        dic = {}
        for i in bondpairs:
            dic.setdefault(i[0],[]).append(i[1])
        max = len(dic.keys())
        # set up a comprehensive list of bondpairs
        reduced_rules = [i for i in itertools.combinations(bondpairs, max)]

        reduced_rules = self.remove_duplicates(reduced_rules)
        # remove rules which don't contain specifications for all
        # bonds in the SBUs
        #bonds = [(i, j) for i in range(len(dataset)) for j in 
        #         range(len(dataset[i].connectpoints))]
        #reduced_rules = [i for i in reduced_rules if 
        #                    self.complete_rules(bonds, i)]
        # set up a dictionary of bondpairs (redundancy included for 
        # simplification later on.
        rules = [{} for i in range(len(reduced_rules))]
        for idx, rule in enumerate(reduced_rules):
            for pair in rule:
                rules[idx].setdefault(pair[0],[]).append(pair[1])
        self.bondrules = rules
        # NOTE: rules should contain inconsistencies

    def complete_rules(self, bonds, rule):
        """
        Check if rule contains all the bonds in the list "bonds"
        """
        crapshoot = [(i[0], i[1]) for i in itertools.chain.from_iterable(rule)]
        for i in bonds:
            if i not in crapshoot:
                return False
        return True

    def remove_duplicates(self, bondcombo):
        """
        Remove entries in the bonding list where the same bond is
        referenced.
        """
        return [i for i in bondcombo if not self.dup_flag(i)]

    def dup_flag(self, seq):
        """ Return True if duplicates of first entry in list """
        seen = set()
        seen_add = seen.add
        for x in seq:
            # awkward use of if statement, but the double negative
            # seems to be the only way this works.
            if x[0] not in seen and not seen_add(x[0]):
                pass
            else:
                return True
        return False

    def bonding_allowed(self, pair, dataset):
        """ Determine if bonding is allowed for the bond pair """
        sbu1 = pair[0][0]
        bond1 = pair[0][1]
        angle1 = pair[0][2]

        sbu2 = pair[1][0]
        bond2 = pair[1][1]
        angle2 = pair[1][2]

        # special bonding
        if dataset[sbu1].symmetrytype[bond1] < 0:
            if dataset[sbu2].symmetrytype[bond1] > 0:
                return False
            else:
                # check if bonds are compatible
                return True
        if dataset[sbu1].ismetal == dataset[sbu2].ismetal:
            return False

        return True

    def database_generation(self):
        """Generates MOFs from a database of SBUs"""
        database = Database()
        database.build_databases()
        stopwatch = Time()
        stopwatch.timestamp()
        # generate SBUs
        metals = [SBU(pdbfile=i) for i in database.metals]
        orgncs = [SBU(pdbfile=i, ismetal=False) for i in database.organic]
        fnlgrp = self.sort([Functional_group(file=i) for i in database.fnlgrps])

        orgncs = self.sort(orgncs)
        metals = self.sort(metals)
        #fnlgrp = [None] + fnlgrp
        # Tom asked for a no functionalization routine.  
        fnlgrp = [None]

        for i in metals:
            for j in orgncs:
                for k in fnlgrp:
                    if k == None:
                        indices = [i.index, j.index, 0]
                    else:
                        indices = [i.index, j.index, k.index]
                    done = False
                    self.unique_bondtypes([i, j])
                    while not done:
                        met = copy.deepcopy(i)
                        org = copy.deepcopy(j)
                        fnl = copy.deepcopy(k)
                        if fnl is not None:
                            # chose a random number of H's to switch
                            # with fnl group
                            sites = j.choose_fnl_sites(fnl)
                            if sites is None:
                                done = True
                            else:
                                # attach functional groups
                                org.add_functional_group(fnl, sites)
                        dataset = [met, org]
                        # build MOF
                        if self.stringhist.get((i.index,j.index)) is not None:
                            if self.stringhist.get((i.index, j.index)) == "Bad":
                                done = True
                            else:
                                self.apply_strings(dataset, indices)
                        else:
                            self.exhaustive_generation(dataset, indices)
                        # terminate if no functional groups were added
                        if fnl is None:
                            done = True
                        # terminate if up to 4 structures were
                        # generated with functional groups. 
                        if self.structcounter == 4:
                            self.structcounter = 0
                            done = True

        stopwatch.timestamp()
        info("Genstruct finished. Timing reports %f seconds."%(stopwatch.timer))

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
                #struct.final_coords()
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

    def unique_bondtypes(self, dataset):
        """
        Determine all the unique possible bonding which can occur
        with the dataset
        """
        self.bondtypes={}
        datalength = len(dataset)
        # outter loops: sbus
        bondtype = 0
        for idx1 in range(datalength):
            for idx2 in range(idx1, datalength):
                sbut1 = dataset[idx1]
                sbut2 = dataset[idx2]
                # middle loops: their bond types
                for bdx1, bondt1 in enumerate(sbut1.symmetrytype):
                    for bdx2, bondt2 in enumerate(sbut2.symmetrytype):
                        if self.bondtype_valid(idx1, idx2, bdx1, bdx2, dataset):
                            # inner loops: their angles
                            for iang1, angt1 in enumerate(sbut1.connectangles[bdx1]):
                                for iang2, angt2 in enumerate(sbut2.connectangles[bdx2]):
                                    # check if the bond type already exists
                                    arrayt = [[idx1, bondt1, iang1],
                                        [idx2, bondt2, iang2]]
                                    if not self.bondtype_exists(arrayt):
                                        bondtype += 1
                                        self.bondtypes[bondtype] = [
                                            [idx1, bondt1, iang1],
                                            [idx2, bondt2, iang2]]

    def bondtype_valid(self, idx1, idx2, bnd1, bnd2, dataset):
        """Checks for bond type validity."""
        sbu1 = dataset[idx1]
        sbu2 = dataset[idx2]

        if (bnd1 in sbu1.special_bond) and (bnd2 in sbu2.special_bond):
            if sbu1.symmetrytype[bnd1] != sbu2.symmetrytype[bnd2]:
                return True
        elif (bnd1 not in sbu1.special_bond) and (bnd2 not in 
                sbu2.special_bond):
            if not (sbu1.ismetal == sbu2.ismetal):
                return True
        return False


    def bondtype_exists(self, arrayt):
        """
        check for existance of unique bond type
        """
        if arrayt in [i for i in self.bondtypes.itervalues()]:
            return True
        else:
            return False
            

    def generate_stringset(self, dataset):
        """
        Encode the set of strings to build a MOF
        """
        stringlist = []
        self.strings = []
        count = 0
        string = "0-0-0-0-0"
        maxstring = self.setmax(dataset)
        done=False
        while not done:
            if count == 1e5:
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

    def apply_bondrule(self, dataset, rule):
        """
        build a MOF based on pre-defined rules. Return True if structure
        was built, otherwise False.
        """
        struct = Structure(dataset)
        self.random_insert(0, dataset, struct)

        # iterate through SBUs
        for sid, sbu in enumerate(struct.mof):
            # iterate through bonds
            for bid, bond in enumerate(sbu.connectpoints):
                cart = choice(sbu.coordinates)
                if struct.connectivity[sid][bid] is None:
                    # determine bonding rule
                    tid = struct.mof[sid].type
                    key = [(tid, bid, aid) for aid in 
                        range(len(sbu.connectangles[bid])) if 
                        rule.has_key((tid,bid,aid))]
                    # random choice of bonding partner if > 1 in the 
                    # bondrule dictionary
                    pair = choice(rule[key[0]])
                    # TODO(pboyd): attach sbus, check for overlaps
                    # check for bonding 
                    # (with rules applied: [strict, leanient])
                    newsbu = len(struct.mof)
                    struct.apply_rule((sid,bid,aid), pair)
                    if struct.overlap_allcheck():
                        # end and try a different rule
                        return False 
                    else:
                        struct.join_sbus(sid, bid, newsbu, pair[1], True)
                    struct.sbu_check(newsbu)
                if struct.saturated() and struct.complete_box():
                    struct.final_coords()
                    struct.mof_reorient()
                    final_coords = struct.getfinals()
                    dump = struct.coordinate_dump()
                    filename = "%i_struct"%(self.nstructs)
                    for i in dataset:
                        filename = filename + "_%i"%(i.index)
                    info("Structure generated!")
                    return True
                    #write_pdb(self.outdir+filename, dump[0], final_coords, 
                    #          newstruct.acell)
                    #newstruct.xyz_debug(self.outdir+filename)
                elif len(struct.mof) == 30:
                    info("Number of SBU's reached 50 without forming"
                            + " a MOF. Exiting..")
                    return False

    def exhaustive_generation(self, dataset, indices):
        """
        generates MOFs via string iteration.
        The code will first generate a list of strings, intelligently
        built to generate structures based on the dataset.  Then each
        set of strings will be attempted until 5 structures are built
        or it runs out of strings.
        """
        # call the Structure class with a set of SBU's and rules to
        # attempt to build the MOF
        self.moflib = []
        self.bondhist = []
        self.moflib.append(Structure(dataset))
        # Start the string off at "0-0-0-0-0"
        string = "0-0-0-0-0"
        self.random_insert(0, dataset, self.moflib[0])
        done = False
        # count the number of strings iterated 
        stringcount = 0
        maxstring = self.setmax(dataset)
        while not done:
            newlist = []
            # list of MOF structures to purge from self.moflib
            purge = []
            for ind, struct in enumerate(self.moflib):
                # check to see if the string will apply to the
                # structure.
                if not self.valid_struct(string, ind):
                    purge.append(ind)
                    stringtype = False
                else:
                    newstruct = copy.deepcopy(struct)
                    stringtype = self.valid_string(string, newstruct,
                                dataset)
                if stringtype == "True":
                    dir = [int(i) for i in string.split("-")]
                    newstruct.sbu_check(dir[0])
                    info("Applying %s"%(string))
                    debug("Number of SBUs in structure: %i"%(
                          len(newstruct.mof)))
                    debug("Size of MOF array: %i"%(len(self.moflib)+
                                                    len(newlist)))
                    sbu2 = len(newstruct.mof)
                    newstruct.apply_string(string)
                    # Temp join SBUs so that the overlap algorithm
                    # can recognize the sbus as being bound together.
                    newstruct.join_sbus(dir[0], dir[1], sbu2,
                                            dir[3], True)
                    #if newstruct.bad_addition(string, newstruct):
                    if newstruct.overlap_allcheck():
                        newstruct.disjoin_sbus(dir[0], dir[1], sbu2,
                                                dir[3])
                    else:
                        newstruct.sbu_check(sbu2)
                        # store the bonding type in history lists
                        bondtype = self.determine_bondtype(string,
                                struct, dataset)
                        newstruct.bondhist.append(bondtype)
                        newstruct.stringhist.append(string)
                        self.bondhist.append(newstruct.bondhist)
                        newlist.append(newstruct)
               
                if newstruct.saturated() and newstruct.complete_box():
                    newstruct.final_coords()
                    newstruct.mof_reorient()
                    final_coords = newstruct.getfinals()
                    dump = newstruct.coordinate_dump()
                    self.nstructs += 1
                    filename = "%06i_struct"%(self.nstructs)
                    for i in indices:
                        filename = filename + "_%i"%(i)
                    self.structcounter += 1
                    info("Structure generated!")
                    info("%s"%(newstruct.stringhist))
                    write_pdb(self.outdir+filename, dump[0], final_coords, 
                              newstruct.acell)
                    #newstruct.xyz_debug(self.outdir+filename)
                    self.stringhist[tuple(indices[:-1])] = newstruct.stringhist
                    # just terminate
                    done = True
                    break
                # Terminate if the number of possible mof structures 
                # gets too big.
                if len(self.moflib) >= 400:
                    self.stringhist[tuple(indices[:-1])] = "Bad"
                    done = True
                    break
                # Terminate if the size of the MOF gets too big.
                if len(newstruct.mof) > 20:
                    self.stringhist[tuple(indices[:-1])] = "Bad"
                    done = True
                    break
            for i in reversed(purge):
                self.moflib.pop(i)
            if len(self.moflib) == 0:
                self.stringhist[tuple(indices[:-1])] = "Bad"
                done = True

            stringcount+=1
            string = self.iterate_string(string, maxstring)
            [self.moflib.append(i) for i in newlist]

    def valid_struct(self, string, ind):
        """
        determine if the structure needs to be deleted.
        """
        ints = [int(i) for i in string.split("-")]
        if ints[0] >= len(self.moflib[ind].mof):
            return False
        return True


    def branched_generation(self, dataset, indices):

        # generate a complete list of strings to sample
        self.generate_stringset(dataset)
        # correction to bondhist such that the routines work as they
        # do for the exhaustive_generation routine
        self.bondhist = [[]]
        structcount = 0
        # initialize timing
        stopwatch = Time()
        # only one structure generated (a "leaf")
        done = False
        struct = Structure(dataset)
        self.random_insert(0, dataset, struct)
        struct.xyz_debug()
        dump = struct.coordinate_dump()
        write_xyz("history", dump[0], dump[1], struct.cell, struct.origins)
        step = 0
        # start the timer
        while not done:
            string = self.strings[step]
            stringtype = self.valid_string(string, struct, dataset)
            ints = [int(i) for i in string.split("-")]
            if stringtype == "True":
                debug("Applying %s"%(string))
                debug("Number of SBUs in structure: %i"%
                    (len(struct.mof)))
                sbu2 = len(struct.mof)
                struct.apply_string(string)
                if struct.bad_addition(string, struct):
                    warning("SBU overlap occured")
                    # remove the SBU from the list, try an incremented
                    # string.
                    struct.mof.pop()
                    struct.connectivity.pop()

                else:
                    bondtype = self.determine_bondtype(string, struct, dataset)
                    struct.bondhist.append(bondtype)
                    self.bondhist[0].append(bondtype)
                    struct.join_sbus(ints[0], ints[1], sbu2, ints[3], True)
                    struct.sbu_check(sbu2)
                    struct.stringhist.append(string)

                struct.complete_box()
                if struct.saturated() and struct.complete_box():
                    #struct.final_coords()
                    struct.getfractionals()
                    struct.mof_reorient()
                    final_coords = struct.getfinals()
                    dump = struct.coordinate_dump()
                    self.nstructs += 1
                    filename = "%i_struct"%(self.nstructs)
                    for i in indices:
                        filename = filename + "_%i"%(i)
                    write_pdb(self.outdir+filename, dump[0], final_coords, 
                              struct.acell)
                    write_xyz(self.outdir+filename, dump[0], final_coords, 
                        struct.cell, struct.origins)
                    # stop the timer
                    self.stringhist[tuple(indices[:-1])] = struct.stringhist
                    info("Structure generated!")
                    done = True

                elif struct.saturated() and not struct.complete_box():
                    info("No structure generated.  Not 3-periodic")
                    done = True
            step += 1
            if stringtype == "Backup":
                done = True
                # re-initialize the structure
                #struct = Structure(dataset)
                # re-apply the existing strings
                #self.readjust(struct, it, dataset)

            if step == 64000:
                done = True

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

    def teststring(self, string, strings):
        """
        Check to see if the string can be applied to the set of
        strings
        """


    def readjust(self, struct, it, dataset):
        """
        Restart the structure generation up to the point where
        the string was successful
        """
        # Reset struct to a blank slate
        struct.mof.append(copy.deepcopy(dataset[0]))
        struct.connectivity.append([None]*len(dataset[0].connectpoints))
        for oldstring in self.stringhist[:it]:
            oldstring = self.strings[oldstring]
            newsbu = len(struct.mof)
            oldints = [int(i) for i in oldstring.split("-")]
            # Rebuild with the good strings
            struct.apply_string(oldstring)
            struct.join_sbus(oldints[0], oldints[1], 
                             newsbu, oldints[3], True)
            struct.sbu_check(newsbu)


    def valid_string(self, string, struct, database):
        """
        determines whether a string will be valid for the next
        addition of an SBU.
        """
        ints = [int(i) for i in string.split("-")]
        # all possibilities are exhausted without joining two SBUs
        # then go back to the previous SBU and change it.
        # TODO(pboyd): this needs to be hashed out.

        if ints[0] >= len(struct.mof):
            return "Backup"
      
        # test to see if the string is at the last bond in the
        # structure.
        if ints[1] > len(struct.mof[ints[0]].connectpoints)-1:
            bondcount = [i for i in struct.connectivity[ints[0]] if
                         i is not None]
            # this means that all bonds have been tried but none have
            # succeeded in placing a new SBU.
            if len(bondcount) <= 1:
                return "Backup"
           
            # this means that the string has gone too far.
            if ints[1] >= len(struct.mof[ints[0]].connectpoints):
                return "False"

        # TODO(pboyd): if no bond is made for one connection point 
        # (ie, the string eventually goes beyond the bond without 
        # attaching a SBU) then terminate the generation
        if (ints[2] == len(database)-1) and \
            (ints[3] == len(database[ints[2]].connect_vector) - 1):
            if struct.connectivity[ints[0]][ints[1]] is None:
                #pass
                # technically illegal because it skips a series of
                # MOF sampling but speeds up the algorithm
                return "Backup"

        if ints[2] >= len(database):
            return "False"

        if ints[3] >= len(database[ints[2]].connectpoints):
            return "False"
        
        if struct.connectivity[ints[0]][ints[1]] is not None:
            return "False"
        # return false if not a special bond
        if(ints[1] in struct.mof[ints[0]].special_bond) and \
             (ints[3] in database[ints[2]].special_bond):
            if struct.mof[ints[0]].ismetal == database[ints[2]].ismetal:
                # return false if the same bond
                 if (struct.mof[ints[0]].symmetrytype[ints[1]] 
                     == database[ints[2]].symmetrytype[ints[3]]):
                     return "False"
            else:
                return "False"
        elif(ints[1] in struct.mof[ints[0]].special_bond) != \
                (ints[3] in database[ints[2]].special_bond):
            return "False"

        elif(ints[1] not in struct.mof[ints[0]].special_bond) and \
                (ints[3] not in database[ints[2]].special_bond):
            if struct.mof[ints[0]].ismetal == database[ints[2]].ismetal:
                return "False"

        if struct.mof[ints[0]].ismetal:
            if ints[4] >= len(database[ints[2]].connectangles[ints[3]]):
                return "False"
        else:
            if ints[4] >= len(struct.mof[ints[0]].connectangles[ints[1]]):
                return "False"

        bondtype = self.determine_bondtype(string, struct, database)
        # check to see if this type of bond has already been made
        # in a structure.
        structhist = struct.bondhist[:]
        structhist.append(bondtype)
        if structhist in self.bondhist:
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
        
        return "%i-%i-%i-%i-%i"%(tuple(liststring))
 
    def setmax(self, dataset):
        """
        Set maximum string
        """
        maxangle, maxbond = 0, 0
        maxsbu = len(dataset) - 1
        for sbu in dataset:
            anglelen = max([(len(i)-1) for i in sbu.connectangles])
            numbonds = len(sbu.connectpoints)-1
            maxangle = max(maxangle, anglelen)
            maxbond = max(maxbond, numbonds)
            # TESTING PURPOSES
            #maxbond = 1
            #maxangle =1
        return([self.nsbumax, maxbond, maxsbu, maxbond, maxangle])


    def determine_bondtype(self, string, struct, dataset):
        """
        Determine the type of bond formed using the string and the
        structure.
        """
        ints = [int(i) for i in string.split("-")]
        sbutype1 = struct.mof[ints[0]].index
        sbutype2 = ints[2]
        bondtype1 = struct.mof[ints[0]].symmetrytype[ints[1]]
        bondtype2 = dataset[ints[2]].symmetrytype[ints[3]]
        angletype = ints[4]
        if struct.mof[ints[0]].ismetal:
            bonding = [[sbutype1, bondtype1, 0],[sbutype2, bondtype2, angletype]]
        else:
            bonding = [[sbutype1, bondtype1, angletype],[sbutype2, bondtype2, 0]]

        bonding.sort()
        bondtype = [k for k, v in self.bondtypes.iteritems() if v == bonding]
        # if the bondtype doesn't match then most likely it's a metal - metal or 
        # organic - organic bond.
        if len(bondtype) == 0:
            return None
        return bondtype[0]

class SBU(object):
    """
    """

    def __init__(self, pdbfile=None, ismetal=None):
        self.vectors = []
        self.points = []
        # each connective bond is assigned a type.  If bonds are 
        # different they should have different types.
        self.symmetrytype = []
        # angles for possible rotations
        self.connectangles = []
        self.name = None
        # index is used to keep track of the SBU's order when building
        # new structures.  This should be found in the pdb file.
        self.index = 0
        # internal index to keep track of the bonding rules
        self.type = None
        self.atomlabel = []
        self.coordinates = []
        if pdbfile is not None:
            self.from_pdb(pdbfile)
        # store history of functional group replacements
        self.fnlgrpchoice = []
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
        # index of the atoms which form connecting bonds with
        # other SBUs
        self.bondingatoms = []
        # angle vectors used to align SBUs once bonded together
        self.anglevect = []
        # store special bonding indices here
        self.special_bond = []

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
        # determine the number of H's in the SBU
        self.hydrogens = [i for i in range(len(self.atomlabel)) if 
                        self.atomlabel[i] == "H"]

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

    def from_xsd(self, filename):
        """
        DEPRECIATED
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
        top = zeros1
        bottom = 0.
        for i in range(len(self.atomlabel)):
            massvect = scalar_mult(WEIGHT[self.atomlabel[i]], self.coordinates[i])
            top = vect_add(top, massvect)
            bottom += WEIGHT[self.atomlabel[i]]
        self.mass = bottom
        return scalar_div(bottom, top)

    def COM_shift(self):
        """
        Shifts all the coordinates such that the centre of mass (COM)
        is at the origin.
        This will make it easier for geometry shifts and rotations
        """
        for icoord, xyz in enumerate(self.coordinates):
            self.coordinates[icoord] = vect_sub(xyz, self.COM)

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

class Functional_group(object):
    """
    Class to contain all information about a functional group
    """
    
    def __init__(self, file=None):
        self.index = None
        self.name = None
        self.atomlabel = []
        self.coordinates = []
        self.connectpoint = np.zeros(3)
        self.connect_vector = np.zeros(3) 

        if (file is not None):
            self.from_pdb(file)
            structure = coord_matrix(self.atomlabel, self.coordinates)
            self.evaluate_bonding(structure[1])

    def evaluate_bonding(self, bonding):
        """Determines the connectpoint, connect_vector"""
        purgeatom = []
        for it, bond in enumerate(bonding):
            if self.atomlabel[it] == "J":
                vect = vect_sub(self.coordinates[it],
                                self.coordinates[bond[0]])
                self.connect_vector = vect
                self.connectpoint = self.coordinates[bond[0]]
                purgeatom.append(it)

        for j in purgeatom:
            self.atomlabel.pop(j)
            self.coordinates.pop(j)

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
            if line.lower().startswith('atom'):
                self.atomlabel.append(line[12:16].strip())
                self.coordinates.append([float(line[30:38]), \
                    float(line[38:46]), float(line[47:54])])
        pdbfile.close()


class Structure(object):
    """
    Structure contains all the machinery necessary to build a MOF from
    simple organic and metallic building units

    """
 
    def __init__(self, sbu_array):
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

        self.bondhist = []
        self.fcoords = []

        # periodic stuff
        self.cell = np.zeros((3,3))
        self.origins = np.zeros((3,3)) 
        # cell parameters
        self.acell = np.zeros(6)
        # cell with unit vectors
        self.ncell = np.zeros((3,3))
        # inverted cell
        self.icell = zeros3
        # index to keep track of number of vectors added
        self.pbcindex = -1

        # history of how the structure was made
        self.stringhist = []
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
        if self.mof[sbu].ismetal:
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

    def bad_addition(self, string, struct):
        """apply checks to see if the new addition is valid"""
        strlist = string.split("-")

        # assumes sbu2 is the newest addition to the MOF
        sbu1 = int(strlist[0])
        bond1 = int(strlist[1])
        sbu2 = len(self.mof) - 1
        bond2 = int(strlist[3])
        # checks overlap over all atoms
        if self.overlap_allcheck():
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
            tol = 1 
        # TODO(pboyd):check for alignment (angle_vect + angles)??

        sbu1 = self.mof[sbu_ind1]
        # local bonding tolerance is low
        for sbuind, sbu in enumerate(self.mof):
            for bondind, bond in enumerate(sbu.connectpoints):
                for bondind1, bond1 in enumerate(
                        sbu1.connectpoints):
                    test = self.connectivity[sbu_ind1][bondind1] is None \
                            and self.connectivity[sbuind][bondind] is None
                    if sbuind != sbu_ind1 and test:
                        point1 = self.apply_pbc(bond1)
                        point2 = self.apply_pbc(bond)
                        if length(vect_sub(point1,point2)) <= tol:
                            info("Local bond found between " +
                                "SBU %i, bond %i, "%(sbu_ind1, bondind1)+
                                "and SBU %i, bond %i."%(sbuind, bondind)+
                                " Found with the application of pbc's.")
                            self.join_sbus(sbu_ind1, bondind1, sbuind, 
                                            bondind, True)
        # return dictionary of anti-parallel, unsaturated bonds within MOF
        # note, may be a problem if bonds are close but not exactly
        # anti-parallel
        bonds = self.antiparallel_bonds(sbu_ind1)
        # generate vectors for all the eligible bonds 
        bondvects = self.gen_bond_vectors(sbu_ind1, bonds)

        bondlist = [(i,tuple(j)) for i in bonds for j in bonds[i]]
        # !!!problems if connection site forms multiple possible bonds!!!
        # check to see if any bonds are formed
        if not bonds:
            info("No new bonds formed with SBU %i"%(sbu_ind1))
            return

        for bond in bondlist:
            if length(bondvects[bond]) <= tol:
                # bonds are close
                info("Local bond found between "+
                        "SBU %i, bond %i, "%(sbu_ind1, bond[0])+
                        "and SBU %i, bond %i."%(bond[1][0],bond[1][1]))
                self.join_sbus(sbu_ind1, bond[0], bond[1][0], 
                               bond[1][1], True)
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
            # test to see if the two sbu's can be joined
            if self.periodic_vector(sbu_ind1, bond, bondvects[bond]):
                if self.connectivity[sbu_ind1][bond[0]] is None and \
                        self.connectivity[bond[1][0]][bond[1][1]] is \
                        None and self.pbcindex < 2:
                    info("New periodic boundary formed between "+
                    "SBU %i, bond %i, "%(sbu_ind1, bond[0])+
                    "and SBU %i, bond %i "%(bond[1][0],bond[1][1])+
                    " with the vector (%f, %f, %f)."%(tuple(bondvects[bond])))
                    self.join_sbus(sbu_ind1, bond[0], bond[1][0], 
                              bond[1][1], True)
                    self.add_vector(bondvects[bond])
                    self.origins[self.pbcindex] = \
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
            if self.connectivity[sbu_ind1][bond[0]] is None and \
               self.connectivity[bond[1][0]][bond[1][1]] is None:
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
                if self.remote_attach_sbus(bondvects[bond]):
                #if self.remote_attach_sbus(pbc_shifted):
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

    def periodic_vector(self, sbu1, bond, bondvector):
        """
        check to see if two unsaturated connecting points can be
        joined by a periodic vector.  There is a parallel test for
        the Y vector (alignment vector) as well, so more unique
        structures can be obtained.
        """
        # establish which SBUs are tyring to be joined
        sbu1 = self.mof[sbu1]
        bond1 = bond[0]
        sbu2 = self.mof[bond[1][0]]
        bond2 = bond[1][1]

        # their anti-parallel X-vectors has already been established.
        # checking now for aligned Y-vectors (alignment).
        vector1 = normalize(sbu1.anglevect[bond1])
        vector2 = normalize(sbu2.anglevect[bond2])

        # test by dot product, the absolute value should = 1
        dotprod = abs(dot(vector1, vector2))
        dottest = np.allclose(dotprod, 1., atol=0.1)

        # test by cross product, the vector should be (0,0,0)
        xprod = cross(vector1, vector2)
        xtest = np.allclose(xprod, np.zeros(3), atol=0.1)
        # FIXME(pboyd): temporary fix for the fact that this constraint
        # won't work for chain-type structures like In3+
        padd_test = "paddlewheel" in self.sbu_array[0].name
        if padd_test:
            if dottest and xtest:
            # check to see if the vector can be added to the 
            # existing lattice vectors
                if self.valid_vector(bondvector):
                    return True
        else:
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
                connecting_point1 = self.mof[sbu1].connectpoints[bond1[0]]
                connecting_point2 = self.mof[sbu2[0]].connectpoints[sbu2[1]]
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
        bonds = [i for i in range(len(self.connectivity[sbu1])) 
                 if self.connectivity[sbu1][i] is None]
        # remove sbu1 from scan
        mofscan = [i for i in range(len(self.mof))]
        bonding_dic = {}
        
        for sbu2 in mofscan:
            bonds2 = [i for i in range(len(self.connectivity[sbu2]))
                      if self.connectivity[sbu2][i] is None]
            for bond2 in bonds2:
                for bond1 in bonds:
                    if (self.validbond(sbu1, sbu2, bond1, bond2)):
                        vect1 = self.mof[sbu1].connect_vector[bond1]
                        vect2 = self.mof[sbu2].connect_vector[bond2]
                        # note the tolerance here will dictate how many
                        # structures will be found
                        if antiparallel_test(vect1, vect2, tol=0.25):
                            bonding_dic.setdefault(bond1,[]).append([sbu2,bond2])
        return bonding_dic

    def validbond(self, idx1, idx2, bnd1, bnd2):
        """Checks for bond type validity."""
        sbu1 = self.mof[idx1]
        sbu2 = self.mof[idx2]

        if (bnd1 in sbu1.special_bond) and (bnd2 in sbu2.special_bond):
            if sbu1.symmetrytype[bnd1] != sbu2.symmetrytype[bnd2]:
                return True
        elif (bnd1 not in sbu1.special_bond) and (bnd2 not in 
                sbu2.special_bond):
            if not (sbu1.ismetal == sbu2.ismetal):
                return True
        return False

    def overlap_allcheck(self):
        """
        checks for all pairwise atomistic overlaps in a growing MOF
        """
        # tolerance in angstroms
        tol = 1.5
        # TODO(pboyd): include minimum image conv. when determining
        # distances.
        for idx1, sbu in enumerate(self.mof):
            for idx2 in range(idx1+1, len(self.mof)):
                bonded_atoms = (None, None)
                if self.bonded(idx1, idx2):
                    # return atoms which form bond
                    bonded_atoms = self.bonded_atoms(idx1, idx2)
                sbu2 = self.mof[idx2]
                distmat = distance.cdist(sbu.coordinates, sbu2.coordinates)
                #distmat = self.min_imgconv(distmat)
                check = [(i,j) for i in range(len(distmat)) for 
                        j in range(len(distmat[i])) if 
                        distmat[i][j] <= tol and (bonded_atoms != (i,j))]
                if len(check) > 0:
                    for ovlp in check:
                        warning(
                        "overlap found between SBU %i (%s), atom %i, "%
                        (idx1, sbu.name, ovlp[0]) + "%s and SBU %i (%s),"
                        %(sbu.atomlabel[ovlp[0]], idx2, sbu2.name) +
                        "atom %i, %s."%(ovlp[1], sbu2.atomlabel[ovlp[1]]))
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

    def min_imgconv(self, distmat):
        """
        Takes a distance matrix and reduces the distances to the 
        minimum image convention.
        """
        # Step 1 convert distances to fractional coordinates.  If 
        # there are no periodic boundaries forget this entire routine.

        fdist = fractional(distmat)

        # Step 2 apply minimum image convention to fractional coordinates

        #distmat = minimg.

        # Step 3 convert back to cartesian coordinates

        distmat = distmat * cellvectors

        return distmat

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
            %(self.mof[sbu2].name, sbu2, bond2, sbu1,
              self.mof[sbu1].name, bond1))
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

        sbu1 = self.mof[sbu1]
        sbu2 = self.mof[sbu2]

        axis = normalize(sbu2.connect_vector[bond2])
        angle = calc_angle(sbu2.anglevect[bond2],
                           sbu1.anglevect[bond1])

        if np.allclose(angle, 0., atol=0.01):
            return
        # origin of rotation
        rotpoint = sbu2.connectpoints[bond2] 

        # test to see if the rotation is in the right direction
        R = rotation_matrix(axis, angle)
        test_vect = matrx_mult(sbu2.anglevect[bond2], R)
        test_angle = calc_angle(sbu1.anglevect[bond1],
                           test_vect)
        # FIXME(pboyd): the tolerance below is quite large,
        if not np.allclose(test_angle, 0., atol=1e-1):
            axis = scalar_mult(-1.,axis)
        # rotation of SBU
        self.sbu_rotation(sbu2, rotpoint, axis, angle + rotangle)
        angle2 = calc_angle(sbu1.anglevect[bond1],
                           sbu2.anglevect[bond2])
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
                            scalar_mult(-1., sbu1.connect_vector[bond1]))
        angle = calc_angle(sbu2.connect_vector[bond2], 
                           scalar_mult(-1.,sbu1.connect_vector[bond1]))
        # origin of rotation
        rotpoint = sbu2.centre_of_mass()

        while (np.allclose(axis, np.zeros(3), atol=1e-4)):
            randaxis = normalize([random(), random(), random()])
            randangle = uniform(0., np.pi/3.)
            self.sbu_rotation(sbu2, rotpoint, randaxis, randangle)
            axis = rotation_axis(sbu2.connect_vector[bond2], 
                             scalar_mult(-1.,sbu1.connect_vector[bond1]))
            angle = calc_angle(sbu2.connect_vector[bond2], 
                           scalar_mult(-1.,sbu1.connect_vector[bond1]))

        # rotate sbu2
        self.sbu_rotation(sbu2, rotpoint, axis, angle)

        # shift the coordinates by the shift vectors 
        shiftvector = vect_sub(sbu1.connectpoints[bond1],
                               sbu2.connectpoints[bond2])
        # bond points are then shifted by the shiftvector
        sbu2.connectpoints = [vect_add(i,shiftvector) for i 
                              in sbu2.connectpoints]
        sbu2.coordinates = [vect_add(i,shiftvector) for i in
                            sbu2.coordinates]

    def getfinals(self):
        """
        get final coordinates from fractionals.
        """
        A = self.cell[0]
        B = self.cell[1]
        C = self.cell[2]
        finals = []
        for i in self.fcoords:
            vect = zeros1
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
        sbu.coordinates = [vect_add(matrx_mult(i, R), C) 
                            for i in sbu.coordinates]
        sbu.connectpoints = [vect_add(matrx_mult(i, R), C)
                            for i in sbu.connectpoints]
        # connect_vector and anglevect are ALWAYS centered at the origin!!
        sbu.connect_vector = [matrx_mult(i, R) for i in sbu.connect_vector]
        sbu.anglevect = [matrx_mult(i, R) for i in sbu.anglevect]
    
    def final_coords(self):
        """
        convert all coordinates to the shifted coordinates derived
        from the periodic boundary conditions
        """
        firstsbu = self.mof[0]
        cellcentre = zeros1
        for k in self.cell:
            cellcentre = vect_add(cellcentre, scalar_div(2., k))
        shiftvect = vect_sub(cellcentre, firstsbu.COM)
        for mof in self.mof:
            for ibond, bond in enumerate(mof.connectpoints):
                mof.connectpoints[ibond] = vect_add(shiftvect, bond)
            for icoord, coord in enumerate(mof.coordinates):
                mof.coordinates[icoord] = vect_add(shiftvect, coord)
        
        self.getfractionals()
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

    def getfractionals(self):
        """
        stores the fractional coordinates of the mof
        """
        fractionals = [self.fractional(coord) for sbu in self.mof 
                        for coord in sbu.coordinates ]
        self.fcoords = fractionals
                
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
        

    def apply_pbc(self, vector):
        """
        applies periodic boundary conditions to a vector
        """

        # FIXME(pboyd): the pbcs when the index = 1 appear broken
        A = self.cell[0]
        B = self.cell[1]
        C = self.cell[2]
        if np.allclose(vector, zeros1):
            return vector
        if self.pbcindex == 0:
            vector = vectorshift(vector, A)
        elif self.pbcindex == 1:
            proj_a = project(vector, self.cell[0])
            proj_b = project(vector, self.cell[1])
            a = length(proj_a) / length(self.cell[0])
            b = length(proj_b) / length(self.cell[1])
            if antiparallel_test(proj_a, self.cell[0]):
                a = -1.*a
            if antiparallel_test(proj_b, self.cell[1]):
                b = -1.*b
            # round a and b to 4 decimal places.  Precision can cause
            # problems if a and b are really close to 1, but under
            # so math.floor doesn't subtract by 1.
            a = round(a, 4)
            b = round(b, 4)
            vector = vect_sub(vector, scalar_mult(math.floor(a),A)) 
            vector = vect_sub(vector, scalar_mult(math.floor(b),B))

        elif self.pbcindex == 2:
            (a, b, c) = matrx_mult(vector, self.icell)
            #bca = dot(B, cross(C, A))
            #a = dot(C, cross(vector, B)) / bca
            #b = dot(C, cross(vector, A)) / (-1.*bca)
            #c = dot(B, cross(vector, A)) / bca

            #a = a - math.floor(a)
            #b = b - math.floor(b)
            #c = c - math.floor(c)
            
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
                    list(self.cell[0]) + list(self.cell[1])+ list(self.cell[2]))))
        #atomcount = 3
        #for icell, cell in enumerate(self.cell):
        #    line.append(cellformat%(tuple(["C"] + list(self.origins[icell])
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

    def valid_vector(self, vect):
        """
        checks a vector against existing cell vectors
        """
        for i in range(self.pbcindex+1):
            if self.line_test(i, vect, tol=0.2):
                print "parallel boundary vector exists!"
                return False
        if self.pbcindex == 1:
            if self.plane_test(vect, tol=0.2):
                print "vector lies in the same plane!"
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
    def __init__(self):
        #TODO(pboyd): add flags for directories for the metal and
        # organic linkers.
        self.organic = []
        self.metals = []
        self.fnlgrps = []
        # accepted filetypes. Add to this later?
        # !!! make sure they are 3-letter file extensions!!!
        self.accepted_file = [".pdb", ".xsd"]
        # defaults for the organic and metal directories.
        self.org_dir = "./organic"
        self.met_dir = "./metals"
        self.fnl_dir = "./functional_groups"
    def build_databases(self):
        """populate databases with SBUs"""
        # FIXME(pboyd): Calling a class (SBU) within a class (Database) 
        # from another class (Generate) doesn't work. 
        # organic
        database = self.scan(self.org_dir)
        for i in database:
            file = self.org_dir + '/' + i
            if i[-4:] == ".pdb":
                #self.organic.append(SBU("label", pdbfile=file[:-4], 
                #                        ismetal=False))
                self.organic.append(file[:-4])
            elif i[-4:] == ".xsd":
                #self.organic.append(SBU("label", xsdfile=file[:-4],
                #                        ismetal=False))
                self.organic.append(file[:-4])
        # metals
        database = self.scan(self.met_dir)
        for i in database:
            file = self.met_dir + '/' + i
            if i[-4:] == ".pdb":
                #self.metals.append(SBU("label", pdbfile=file[:-4]))
                self.metals.append(file[:-4])
            elif i[-4:] == ".xsd":
                #self.metals.append(SBU("label", xsdfile=file[:-4])) 
                self.metals.append(file[:-4])
        # functional groups
        database = self.scan(self.fnl_dir)
        for i in database:
            file = self.fnl_dir + '/' + i
            if i[-4:] == ".pdb":
                self.fnlgrps.append(file[:-4])
            elif i[-4:] == ".xsd":
                self.fnlgrps.append(file[:-4])

    def scan(self, dir):
        """Scans the directories for SBUs, returns list of SBU files."""
        dir = dir.strip(".")
        cwd = os.getcwd()
        sbufiles = subprocess.Popen(["ls", cwd + dir], shell=False, 
                stdout=subprocess.PIPE).communicate()[0].split("\n")
        # Note: this only works for file extensions with three letters.
        # any other types (1, 2, 4 etc ...) will not be included properly.
        sbufiles = [i for i in sbufiles if i and i[-4:] in self.accepted_file]
        return sbufiles

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
        return vect_sub(shift_vector, 
                scalar_mult(math.floor(fraction),root_vector))

def determinant(matrix):
    """ calculates the determinant of a 3x3 matrix"""
    a = matrix[0][0]; b = matrix[0][1]; c = matrix[0][2]
    d = matrix[1][0]; e = matrix[1][1]; f = matrix[1][2]
    g = matrix[2][0]; h = matrix[2][1]; i = matrix[2][2]
    det = a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h
    return det

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
    # G is the metal - metal bond for linked chains
    elif("G" in (atom1,atom2) and 
        (("Y" not in (atom1,atom2))and("Z" not in (atom1,atom2)))):
        return 1.1
    elif("X" in (atom1, atom2) and 
        (("Y" not in (atom1,atom2))and("Z" not in (atom1,atom2)))):
        return 0.9
    elif("G" in (atom1,atom2) or "X" in (atom1, atom2))and \
        (("Y" in (atom1,atom2))or("Z" in (atom1,atom2))):
        return 1.1
    elif("Y" in (atom1,atom2))and("G" not in (atom1,atom2) 
            or "X" not in (atom1,atom2)):
        return 0.
    elif("Z" in (atom1,atom2))and("G" not in (atom1,atom2) 
            or "X" not in (atom1,atom2)):
        return 0.
    elif("Br" in (atom1,atom2))and("J" in (atom1,atom2)):
        return 2.0
    elif("I" in (atom1,atom2))and("J" in (atom1,atom2)):
        return 2.5
    elif("Cl" in (atom1,atom2))and("J" in (atom1,atom2)):
        return 1.8 
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
        coord2 = zeros1
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

def write_pdb(label, atoms, coords, acell):
    pdbfile = open('%s.pdb' % label, 'w')
    today = date.today()
    lines = []
    atomformat1 = "%-6s%5i %-4s %3s %1s%4i%s   "
    atomformat2 = "%8.3f%8.3f%8.3f%6.2f%6.2f         %-3s%2s\n"
    lines.append("%6s   GenStruct! created by: Anon.\n" % (
                  'REMARK'))
    lines.append("%6s   Created: %s \n%-6s" % ('REMARK', 
                  today.strftime("%A %d %B %Y"), 'CRYST1'))
    lines.append("%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P1\n"
                  % (tuple(acell)))
    for atom in range(len(atoms)):
        lines.append(atomformat1%('ATOM',atom+1,atoms[atom],"MOL","X",1," ") +
                     atomformat2%(tuple(coords[atom]+[1.,0.,atoms[atom],0])))
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

def main():
    """Default if run as an executable"""

    # initialize logging 
    Log()
    open('history.xyz', 'w')
    open('debug.xyz', 'w')
    genstruct = Generate()
    genstruct.database_generation()
    #genstruct.bondrule_generation()
if __name__ == '__main__':
    main()
