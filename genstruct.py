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
import functional
from config import Options
from operations import *
from numpy import array
from elements import WEIGHT
from bookkeeping import * 
from random import random, uniform, randrange, choice
from datetime import date
from logging import warning, debug, error, info, critical
from scipy.spatial import distance
# modules related to symmetry finding
from atoms import Atoms

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
        self.outdir = "output." + database.dbfile + "/"
        self.bondtypes = {}
        self.gen_dir()
        # list of bonding with bondtypes for successful structures
        self.bondhist = []
        # bonding rules for the bond rule generation method
        self.bondrules = []
        # counter for number of structures generated
        self.mofcount = 0
        # functionalization
        self.functional = functional.Functional_groups()

        # CSV file generation
        self.csv = CSV()
        self.csv.set_name(database.dbfile)
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
                if len(self.complete_mofs) > 0:
                    info("Applying functional groups...")
                    #self.apply_functional_groups()

        stopwatch.timestamp()
        info("Genstruct finished. Timing reports %f seconds."%(stopwatch.timer))

        self.csv.write_file()

    def random_hydrogen_selection(self, sbus_withH):
        """
        Returns a tuple with a random selection of H atoms in
        each SBU with hydrogens in it.
        """

        hydrogen_list = []
        for sbu in sbus_withH:

            numhydrogens = len(sbu.hydrogens)
            # limit the number of hydrogen replacements to 10
            # to prevent HUGE file names.
            if numhydrogens > 10:
                numhydrogens = 10
            randnumhyd = 0
            while randnumhyd == 0:
                randnumhyd = randrange(numhydrogens)
            hdic = {}
            while len(hdic.keys()) < randnumhyd:
                hydrogen = choice(sbu.hydrogens)
                hdic[hydrogen] = "?"

            hydrogens = hdic.keys()
            hydrogens.sort()
            hydrogen_list.append(tuple(hydrogens))

        return tuple(hydrogen_list)

    def apply_functional_groups(self):
        """
        Applies functional groups to built MOFs.
        """
        # determine which hydrogens to replace, what their
        # connecting atoms are (Important for overlap considerations)
        sbus_withH = [sbu for mof in self.complete_mofs for 
                sbu in mof.sbu_array if len(sbu.hydrogens) > 0]

        dic = {}
        # isolate unique SBUs
        for sbu in sbus_withH:
            if sbu.metal and sbu.index == 8 or sbu.index == 2:
                pass
            else:
                dic[sbu.internal_index] = sbu

        sbus_withH = dic.values()

        # return if no hydrogens available for functional swapping
        if len(sbus_withH) == 0:
            return
        # generate a set of hydrogen replacements for each sbu in SBUs
        # to a limit of 1000.  After this many, just quit.
        #hydrogenlist = []
        #for sbu in sbus_withH:
        #    #itertools chain, combinations
        #    s = sbu.hydrogens
        #    templist = \
        #    sorted(itertools.chain.from_iterable(
        #        itertools.combinations(s, r) for r in range(len(s)+1)))
        #    templist = [i for i in templist if i]
        #    hydrogenlist.append(templist)
        
        # itertools.product ensures the first entry is from the first 
        # list and the second from the second list.  This way we can 
        # track which SBU's hydrogens we are manipulating.
        #combine = [i for i in itertools.product(*hydrogenlist)]
        
        groups = self.functional.atoms.keys()
        groups.sort()
        # Iterate through the functional groups in self.functional
        for group in groups:
            # Tally which SBU's in the MOF have hydrogens to replace
            mofnum = len(self.complete_mofs)
            ignore = []; choicecount=0
            done = False
            while not done:
                choicecount += 1
                # randomly choose Hydrogens from each SBU
                
                select_H = self.random_hydrogen_selection(sbus_withH)
                trial = 0
                while select_H in ignore:
                    trial += 1
                    select_H = self.random_hydrogen_selection(sbus_withH)
                    if trial == 3000:
                        done = True
                        break
                ignore.append(select_H)
                #Hind = randrange(len(combine))
                #while Hind in ignore:
                #    Hind = randrange(len(combine))
                #    if len(ignore) == len(combine):
                #        done = True
                #        break
                #select_H = combine[Hind]
                #ignore.append(Hind)
                replace_dic = {}
                for idx, sbu in enumerate(sbus_withH):
                    replace_dic[sbu.internal_index] = select_H[idx]
                # select a structure from self.complete_mofs
                rand = randrange(len(self.complete_mofs))
                newstr = copy.deepcopy(self.complete_mofs[rand])
                ovlp = False
                for idx, sbu in enumerate(newstr.mof):
                    sbutype = newstr.storetype[idx]
                    if sbu.internal_index in replace_dic.keys():
                        # Go through each hydrogen listed in 
                        # replace_dic and swap for the functional group
                        for hydrogen in replace_dic[sbu.internal_index]:
                            connect_atom = sbu.hydrogens_connect[
                                sbu.hydrogens.index(hydrogen)]
                            # do a functional group swap
                            fnlcoords = self.functional_swap(newstr,
                                    idx, hydrogen, connect_atom, group)
                            # check for overlaps
                            if self.overlap_fnlcheck(idx, hydrogen,
                                    connect_atom, group,
                                    self.functional.atoms[group],
                                    fnlcoords, newstr):
                                ovlp = True
                                break
                                #TODO(pboyd): include a rotation scheme
                            else:
                                # append the fractional coordinates to
                                # the mof to include in overlap calcs
                                fcoords = [newstr.fractional(coord) for
                                        coord in fnlcoords]
                                newstr.fcoords[idx] = \
                                        newstr.fcoords[idx] + fcoords[:]
                                # append to MOF
                                newstr.coordinates[idx] = \
                                    newstr.coordinates[idx] + \
                                    fnlcoords[:]
                                newstr.atoms[idx] = newstr.atoms[idx] + \
                                        self.functional.atoms[group][:]

                                newstr.atoms_fftype[idx] = \
                                        newstr.atoms_fftype[idx] + \
                                        self.functional.atoms_fftype[group][:]
                                # add to sbu atoms as well
                                sbu.atom_label[sbutype] += \
                                        self.functional.atoms[group][:]
                                sbusize = len(sbu.table[sbutype])
                                # add the entries for the functional
                                # group. appended to the end.
                                fnlconnectatm = self.functional.\
                                        connect_points[group]
                                bond_length = self.functional.\
                                        bond_length[group]
                                for iatm, fnlcnt in enumerate(
                                    self.functional.table[group]):
                                    # append to the array
                                    if iatm == fnlconnectatm:
                                        # add bonding info where the
                                        # connecting carbon is
                                        bondarray = [0.] * (connect_atom)
                                        # single bond connects the
                                        # functional group...
                                        bondarray += [(bond_length, "S")]
                                        bondarray += [0.] * (sbusize - 
                                                         connect_atom - 1)
                                        bondarray += fnlcnt

                                    else:
                                        bondarray = [0.] * sbusize
                                        bondarray += fnlcnt

                                    sbu.table[sbutype].append(
                                                bondarray)
                                # add info to the existing atoms in the
                                # SBU
                            
                                for jatm in range(sbusize):
                                    if jatm == connect_atom:
                                        bondarray = [0.] * (
                                                fnlconnectatm)
                                        # single bond connects the 
                                        # functional group....
                                        bondarray += [(bond_length, "S")]
                                        bondarray += [0.] * (
                                               len(self.functional.
                                                    table[group]) -
                                                fnlconnectatm - 1)
                                    else:
                                        bondarray = [0.] * len(
                                                self.functional.
                                                table[group])
                                    sbu.table[sbutype][jatm] += bondarray
                    if ovlp:
                        break

                if not ovlp:
                    mofnum += 1
                    # remove all the hydrogens
                    for idx, sbu in enumerate(newstr.mof):
                        sbutype = newstr.storetype[idx]
                        if sbu.internal_index in replace_dic.keys():
                            sortedH = list(replace_dic[sbu.internal_index])
                            sortedH.sort()
                            [newstr.atoms_fftype[idx].pop(hydrogen)
                             for hydrogen in reversed(sortedH)]
                            [newstr.coordinates[idx].pop(hydrogen)
                             for hydrogen in reversed(sortedH)]
                            [newstr.atoms[idx].pop(hydrogen) for 
                             hydrogen in reversed(sortedH)]
                            [newstr.fcoords[idx].pop(hydrogen) for 
                             hydrogen in reversed(sortedH)]
                            #switch out the H with functional group
                            #in the connectivity matrix of the SBU
                            #Remove the entry for hydrogen.
                            [sbu.atom_label[newstr.storetype[idx]].pop(
                                hydrogen) for hydrogen in reversed(
                                    sortedH)]
                            for hydrogen in reversed(sortedH):

                                [i.pop(j) for i in sbu.table[sbutype]
                                     for j in range(len(i)) if 
                                     j == (hydrogen)]
                                [sbu.table[sbutype].pop(i) for i in
                                     range(len(sbu.table[sbutype]))
                                     if i == hydrogen]
                     
                        #print "\n\nNEW SBU"
                        #for id, i in enumerate(sbu.table[sbutype]):
                        #    print newstr.atoms[idx][id], i
                    # re-initialize self.sbu_connections{}
                    newstr.sbu_connections = {}
                    newstr.create_connect_table()
                    self.finalize(newstr, group, replace_dic)
                    # write the file
                if mofnum > 5:
                    done = True
                elif choicecount >= 2000:
                    done = True

    def overlap_fnlcheck(self, sbu, hydrogen, connect_atom, 
                           group, atoms, coordinates, struct):
        """
        Overlap determined by transposing MOF coordinates to the 
        minimum image of the centre of mass of the functional group
        then performing a distance check with cdist.
        """
        # scaling factor for vdw radii overlap 
        sf = 0.6

        fnlconnect = self.functional.connect_points[group]
        fnlcom = centre_of_mass(atoms, coordinates)

        # the fractional coordinates are stored
        mof_coords = struct.min_img_shift(fnlcom)
        #mof_atoms = [struct.atoms[isbu][iatom] for isbu in 
        #                range(len(struct.atoms)) for iatom in 
        #                range(len(struct.atoms[isbu]))]
        #coordsa = [xyz for bu in range(len(mof_coords)) for xyz in 
        #           mof_coords[bu]]
        #coordsa = coordsa + coordinates
        #atm = ["Fe" for i in range(len(atoms))]
        #atm = mof_atoms + atm
        #write_xyz("history", atm, coordsa, zeros3)
        overlap = []
        for idx, coords in enumerate(mof_coords):

            distmat = distance.cdist(coordinates, coords)
            deletion = []
            check = [(i, j, distmat[i][j], 
                Radii[atoms[i]]+Radii[struct.atoms[idx][j]])
                for i in range(len(distmat)) for
                    j in range(len(distmat[i])) 
                if distmat[i][j] <= (Radii[atoms[i]] + 
                    Radii[struct.atoms[idx][j]]) * sf]
            if idx == sbu:
                for ich, ch in enumerate(check):
                    if ch[0] == fnlconnect and ch[1] == connect_atom:
                        deletion.append(ich)
                    if ch[1] == hydrogen:
                        deletion.append(ich)

            deletion.sort()
            for i in reversed(deletion):
                check.pop(i)

            #for thing in check:
            #    info("%i, %i"%(idx, sbu))
            #    info("%i, %i, %i"%(fnlconnect, connect_atom, hydrogen))
            #    info("%i, %i, %f"%(thing[0], thing[1], distmat[thing[0]][thing[1]]))
            overlap = overlap + check

        # apply the conditions of the specific hydrogen and connect atom
        # for the sbu it's attached to. (HARD because fractional coords 
        # are not listed according to SBU)
        
        if len(overlap) > 0:
            return True

        return False

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

    def functional_swap(self, struct, idx, hydrogen, connect_atom, group):
        """ switch a hydrogen for a functional group"""

        sbu = struct.mof[idx]
        hyd_coord = struct.coordinates[idx][hydrogen]
        connect_coord = struct.coordinates[idx][connect_atom]
        hyd_vector = vect_sub(hyd_coord, connect_coord)
        hyd_vector =  normalize(hyd_vector)
        fnl_vector = self.functional.connect_vector[group]

        axis = rotation_axis(fnl_vector, hyd_vector)
        angle = calc_angle(hyd_vector, fnl_vector)
        coordinates = self.rotation(axis, angle, connect_coord, group)
        atm = self.functional.connect_points[group]
        shift = vect_sub(connect_coord, 
                self.functional.coordinates[group][atm])
        bond = scalar_mult(self.functional.bond_length[group], 
                           normalize(hyd_vector))
        shift = vect_add(shift, bond)

        return [vect_add(i, shift) for i in coordinates]
    
    def rotation(self, axis, angle, rotpoint, group):
        """
        rotates a functional group to be properly added to the SBU
        """
        atm = self.functional.connect_points[group]
        rotpoint = self.functional.coordinates[group][atm]
        R = rotation_matrix(axis, angle)
        C = matrx_mult(rotpoint, (matrx_sub(zeros3, R)))
        coordinates = [vect_add(matrx_mult(i, R), C) for i in 
                         self.functional.coordinates[group]]
        return coordinates

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
        #TODO(pboyd): determination of bondtypes should include a 
        #metal and organic with the same index
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

            # if both SBUs are metals, then consider only if there
            # is a bond name called 'intermetal'
            if sbus[0].metal == sbus[1].metal:
                if flag_name[0][bond[0]] == "intermetal" and \
                    flag_name[1][bond[1]] == "intermetal":
                    if flag_sym[0][bond[0]] != flag_sym[1][bond[1]]:
                        # add symmetry type
                        angle = 0
                        # NOTE: changed sbus.index to sbus.internal_index
                        key = ("intermetal", sbus[0].internal_index, 
                                flag_sym[0][bond[0]], 
                                sbus[1].internal_index, 
                                flag_sym[1][bond[1]], angle)
                        value = ("intermetal", sbus[0].internal_index, 
                                 bond[0], 
                                 sbus[1].internal_index, bond[1])

                        self.bondtypes.setdefault(key, [])
                        if value not in self.bondtypes[key]:
                            self.bondtypes[key].append(value)
                        
            # if both SBUs are metals, then consider only if there
            # is a bond name called 'intermetal'
            else:
                if flag_name[0][bond[0]] == flag_name[1][bond[1]]:
                    # range over all angles
                    for angle in range(len(structure.connect_angles)):
                        key = (flag_name[0][bond[0]], 
                                sbus[0].internal_index, 
                                flag_sym[0][bond[0]], 
                                sbus[1].internal_index, 
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
                        value = (flag_name[0][bond[0]], 
                                 sbus[0].internal_index, 
                                 bval1, sbus[1].internal_index, bval2,
                                 angle)
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
        indexset = set([sbu.internal_index for sbu in structure.sbu_array])
        done = False
        while not done:
            #apply the string
            for struct in self.moflib:
                try:
                    string = self.strings[iterstring]
                except:
                    done = True
                if self.valid_string(string, struct):
                    copystruct = copy.deepcopy(struct)
                    sbu = len(copystruct.mof) 
                    copystruct.apply_string(string)
                    #if copystruct.overlap_allcheck():
                    if copystruct.overlap_sbucheck(sbu):
                        copystruct.storetype.pop()
                    else:
                        bondtype = self.determine_bondtype(
                                string, copystruct)
                        copystruct.sbu_check(sbu)
                        copystruct.bondhist.append(bondtype)
                        self.bondhist.append(copystruct.bondhist)
                        copystruct.stringhist.append(string)
                        self.moflib.append(copystruct)
                        # fix for Snurr's big NU100 building units
                        # TODO(pboyd): make this less time consuming
                        for indsbu in range(len(copystruct.mof)):
                            copystruct.sbu_check(indsbu)
                    if copystruct.saturated() and copystruct.complete_box():
                        # check to see if all SBU's are included in 
                        # the structure.
                        base = [i.internal_index for i in 
                                 copystruct.sbu_array]
                        sets = [i[3] for i in copystruct.bondhist]
                        sets += [i[1] for i in copystruct.bondhist]
                        if len(set(sets)) < len(set(base)):
                            pass
                        else:
                            # re-initialize self.sbu_connections{}
                            copystruct.sbu_connections = {}
                            copystruct.create_connect_table()
                            copystruct.final_coords()
                            copystruct.mof_reorient()
                            copystruct.get_cartesians()
                            self.finalize(copystruct, 0)
                            # export the MOF for functional group
                            # placement
                            self.complete_mofs.append(copy.deepcopy(copystruct))
                            info("Structure Generated!")
                            copystruct.origins = zeros3[:]
                            structcount += 1

                    # Give up if number of structures becomes too big.
                    if len(self.moflib) > 12000:
                        done = True
                        break
                    # constraints and conditions related to specific
                    # topologies
                    # just make one for now...
                    if copystruct.name == "garbage":
                    #if copystruct.name == "pcu" or \
                    #   copystruct.name == "nbo":
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

    def finalize(self, struct, idx, hsubs=None):
        """ write the MOF file etc..."""

        # formatting class for calling the symmetry routines
        # TODO(pboyd): this is a bit of a hack job - should clean this.
        self.mofcount += 1
        basename = "str"
        org_ind = []
        for sbu in struct.sbu_array:
            if sbu.metal:
                type = "m"
                met_ind = sbu.index
            else:
                type = "o"
                org_ind.append(sbu.index)

            basename += "_%1s%i"%(type, sbu.index)
            if hsubs is not None:
                if sbu.internal_index in hsubs.keys():
                    sortedH = list(hsubs[sbu.internal_index])
                    sortedH.sort()
                    sortedH = tuple(sortedH)
                    basename += "-sub"
                    for i in sortedH:
                        basename += "-%i"%i
        basename += "_%1s%i"%("f", idx)
        basename += "_%s"%(struct.name)
        filename = self.outdir + basename 
        cif = CIF(struct, struct.master_table)
        cif.write_cif(filename)
        
        sym_name = cif.symmetry.get_space_group_name()
        sym_number = cif.symmetry.get_space_group_number()

        self.csv.add_data(basename,
                metal_index = met_ind,
                organic_index1 = org_ind[0],
                organic_index2 = org_ind[1],
                h_m_symmetry_name=sym_name,
                symmetry_number=sym_number)
        #struct.write_pdb(filename)
        return

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
        elif struct.name == "pto":
            if len(struct.connectivity) > 14:
                return False
        elif struct.name == "nbo":
            if len(struct.connectivity) > 9:
                return False
        elif struct.name == "mtn":
            if len(struct.connectivity) > 9:
                return False
        elif struct.name == "acs":
            if len(struct.connectivity) > 7:
                return False
        elif struct.name == "bcu":
            if len(struct.connectivity) > 6:
                return False
        elif struct.name == "flu":
            if len(struct.connectivity) > 7:
                return False
        elif struct.name == "the":
            if len(struct.connectivity) > 12:
                return False
        elif struct.name == "ubt":
            if len(struct.connectivity) > 7:
                return False
        elif struct.name == "rht":
            if len(struct.connectivity) > 8:
                return False
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
        sbutype1 = struct.mof[idx[0]].internal_index
        sbutype2 = struct.sbu_array[idx[2]].internal_index
        bond1 = idx[1]
        bond2 = idx[3]
        angletype = idx[4]
        # determine the type of bonding taking place.  Ie what is the
        # type of bond joining the two SBUs as determined by the SBU
        # already placed.
        bondname = struct.vectorflag[idx[0]][idx[1]]
        # since internal_index, bondtype is unique, if the bond 
        # exists it is in one of the following two forms in the 
        # dictionary
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
        # index is a global definition which should be defined in the
        # input file.  This should be a unique identifier for a
        # particular SBU so that it can be referenced in a database.
        self.index = 0
        # internal_index is for the purposes of the code, this will 
        # classify each SBU in the dataset individually so that 
        # particular bonding rules can be pre-established between SBUs.
        # this is in case two SBUs have the same index (described above).
        self.internal_index = 0
        # internal index to keep track of the bonding rules
        self.type = None
        # atom_label contains the element symbols
        self.atom_label = {}
        # atoms_fftype contains the forcefield atom types described in
        # the input data
        self.atoms_fftype = {}
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
        self.bondingatoms = {}
        self.COM = {}

        # bondspec is the stored bonding tables under [table] in the
        # input file
        self.bondspecflag = False
        self.bondspec = {}
        # store a connectivity table
        self.table = {} 
        # store indices which are Hydrogen atoms
        self.hydrogens = []
        # store indices of atoms bonded to Hydrogen atoms.
        # in the same order.
        self.hydrogens_connect = []

        self.init_arrays(data)

        self.centre_of_mass()
        self.COM_shift()
        self.find_hydrogen_atoms()
        # identify the atoms in the SBU which are involved in connecting
        # the MOF (for disregarding overlap considerations)
        self.determine_bonding_atoms()

    def determine_bonding_atoms(self):
        """
        Returns the array of atom indices which are involved in
        connecting SBUs.
        """
        # TODO(pboyd): currently ignores bonding between metal SBUs
        # "intermetal" type structures. 
        for name in self.table.keys():
            for idx, bonding in enumerate(self.table[name]):
                if name == "metallic":
                    if self.atom_label[name][idx] == "N":
                        # if the atom is within a certain distance 
                        # of a connect_point, then it is a bonding
                        # atom
                        key = self.connect_points.keys()[0]
                        if [(idx, j) for j in self.connect_points[key]
                                if length(j, self.coordinates[name][idx])
                                    < 1.3]:
                            self.bondingatoms.setdefault(name,[]).\
                                    append(idx)

                        #nbonds = len([i for i in bonding if i != 0.])
                        #if nbonds < 3:
                else:
                    if self.atom_label[name][idx] == "C":
                        # just select one of the bonding parameters at 
                        # random.  They should all start nearly in the 
                        # same spot.
                        param = self.connect_points.keys()
                        try:
                            param.pop(param.index("metallic"))
                        except:
                            pass
                        param = choice(param)
                        if [(idx, j) for j in self.connect_points[param]
                                if length(j, self.coordinates[name][idx])
                                    < 1.1]:
                            self.bondingatoms.setdefault(name,[]).\
                                    append(idx)
                        #nbonds = len([i for i in bonding if i != 0.])
                        #if nbonds < 3:
                        #    self.bondingatoms.setdefault(name,[]).\
                        #            append(idx)
        return

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
            for atom, dist in enumerate(self.table[name][hydrogen]):
                if dist > 0.:
                    atname = self.atom_label[name][atom]
                    self.hydrogens_connect[ind] = atom

    def coord_matrix(self):
        """
        Generate a coordination matrix for the SBU's coordinates.
        """
        for name in self.atom_label.keys():
            numatms = len(self.atom_label[name])
            # populate empty connectivity table
            self.table[name] = [[0.] * numatms for i in xrange(numatms)]
            distmatrx = distance.cdist(self.coordinates[name], 
                                   self.coordinates[name])
            if not self.bondspecflag:
                for i in range(numatms):
                    for j in range(i+1, numatms):
                        tol = bond_tolerance(self.atom_label[name][i], 
                                     self.atom_label[name][j])
                        if (i != j) and (distmatrx[i,j] <= tol):
                            # default to single bonds "S" if not known.
                            self.table[name][i][j] = (distmatrx[i,j], "S")
                            self.table[name][j][i] = (distmatrx[j,i], "S")
            else:
                for i in self.bondspec[name].keys():
                    for j in self.bondspec[name][i]:
                        self.table[name][i][j[0]] = (distmatrx[i,j[0]],
                                                        j[1])
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
                [self.atoms_fftype.setdefault(name, []).append(
                    line.rstrip().split()[1]) for line in lines]
                [self.coordinates.setdefault(name, []).append(
                    [float(elem) for elem in line.rstrip().split()[2:5]]) 
                    for line in lines]
        
        if config.has_section("table"):
            # no need to scan the SBU for distances after this
            section = "table"
            self.bondspecflag = True
            for name in config.options(section):
                lines = io.BytesIO(config.get(section, name)).readlines()
                lines = [line.rstrip() for line in lines if line.rstrip()]
                table = {}
                for line in lines:
                    line = line.split()
                    table.setdefault(int(line[0]), []).append(
                            (int(line[1]), line[2]))
                    table.setdefault(int(line[1]), []).append(
                            (int(line[0]), line[2]))
                self.bondspec[name] = table

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

            self.COM[name] = zeros1[:]

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
        self.fconnect = []
        self.connect_align = []
        self.connect_vector = []
        self.connect_sym = []
        # list to keep track of bonding between SBU's
        self.connectivity = []
        self.vectorflag = []

        self.connect_angles = []
        self.COM = []
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
        
        # store atom force field types
        self.atoms_fftype = []
        # store what type was bonded (ie. metallic, default etc)
        self.storetype = []
        # keeps track of the atoms in each SBU which form bonds with 
        # other SBUs 
        self.sbu_connections = {}
        # master connection table
        self.master_table = {}
        ind=0
        if (name == 'the') and ((self.sbu_array[ind].index == 6)
                or(self.sbu_array[ind].index == 12)) \
                and (sbu_array[ind].metal):
            ind = choice([i for i in range(len(sbu_array))
                              if not sbu_array[i].metal])
        self.seed(name, ind)

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
            tol = 0.6
        # return dictionary of anti-parallel, unsaturated bonds within MOF
        # note, may be a problem if bonds are close but not exactly
        # anti-parallel
        bonds = self.antiparallel_bonds(sbu)
        # generate vectors for all the eligible bonds 
        bondvects = self.gen_bond_vectors(sbu, bonds)
        bondlist = [(i,tuple(j)) for i in bonds for j in bonds[i]]

        if not bonds:
            #info("No new bonds formed with SBU %i"%(sbu))
            return

        for bond in bondlist:
            #info("non adj vector: (%6.3f,%6.3f,%6.3f)"%tuple(bondvects[bond]))
            adjusted_vect = self.min_img(bondvects[bond])
            #info("min img vector: (%6.3f, %6.3f, %6.3f)"%tuple(adjusted_vect))
            #info("length %f"%length(adjusted_vect))
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
        self.COM
        """
        # check for default or metallic names
        size = len(self.connectivity)
        if name == "metallic" and not self.sbu_array[ind].metal:
            self.coordinates.append(self.sbu_array[ind].coordinates[name][:])
            self.atoms.append(self.sbu_array[ind].atom_label[name][:])
            self.COM.append(self.sbu_array[ind].COM[name][:])
            self.storetype.append("metallic")
            self.atoms_fftype.append(self.sbu_array[ind].atoms_fftype[name][:])
        else:
            self.coordinates.append(self.sbu_array[ind].coordinates["default"][:])
            self.atoms.append(self.sbu_array[ind].atom_label["default"][:])
            self.COM.append(self.sbu_array[ind].COM["default"][:])
            self.storetype.append("default")
            self.atoms_fftype.append(self.sbu_array[ind].atoms_fftype["default"][:])

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
        tol = 0.1
        # establish which SBUs are tyring to be joined
        bond1 = bond[0]
        sbu2 = bond[1][0]
        bond2 = bond[1][1]

        # their anti-parallel X-vectors has already been established.
        # checking now for aligned Y-vectors (alignment).
        vector1 = normalize(self.connect_align[sbu1][bond1])
        vector2 = normalize(self.connect_align[sbu2][bond2])
        nbondvect = normalize(bondvector)

        # test by dot product, the absolute value should = 1
        # NOTE: for metals of index 8 and 9, this should NOT
        # be an absolute value
        dotprod = abs(dot(vector1, vector2))
        dottest = np.allclose(dotprod, 1., atol=tol)

        # test by cross product, the vector should be (0,0,0)
        #xprod = cross(vector1, vector2)
        #xtest = np.allclose(xprod, zeros3[:], atol=tol)

        joined_atoms_vect = normalize(self.connect_vector[sbu1][bond1])
        perpcross = cross(joined_atoms_vect, nbondvect)
        perpdot = dot(joined_atoms_vect, nbondvect)
        perptest = np.allclose(length(perpcross), 1., atol=tol) and \
                np.allclose(perpdot, 0., atol=tol)
        #if dottest and xtest and not perptest:
        if dottest and not perptest:
        #if not perptest:
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

    def overlap_sbucheck(self, sbu):
        """
        Checks for pairwise atom overlaps between a particular SBU and
        the rest of the MOF
        """
        #TODO(pboyd): eventually we will need a bonding matrix of 
        # the entire MOF structure instead of these cheap fixes.
        sf = 0.4

        coords = self.coordinates[sbu]
        atoms = self.atoms[sbu]
        # SBU.coordinates is not shifted, only self.coordinates!!!

        #com = self.COM[sbu]

        for atom in coords:
            remaining_coords, remaining_atoms = [], []
            # shift all coordinates to within the min
            # img of the atom coordinates
            if self.pbcindex == 2:
                shifted_coords = self.min_img_shift(atom)
                for index in range(len(self.coordinates)):
                    name = self.storetype[index]
                    try:
                        # get bonding atoms from SBU.  These
                        # atoms will be ignored in the overlap test
                        bonding = self.mof[index].bondingatoms[name]
                    except:
                        # Ba2+ case is special.
                        bonding = []
                    remaining_coords.append([shifted_coords[index][coord]
                            for coord in range(len(self.coordinates[index]))
                            if coord not in bonding and index != sbu])
                    remaining_atoms.append([self.atoms[index][atom]
                        for atom in range(len(self.atoms[index]))
                        if atom not in bonding and index != sbu])
            else:
                for index in range(len(self.coordinates)):
                    name = self.storetype[index]
                    try:
                        # get bonding atoms from SBU.  These
                        # atoms will be ignored in the overlap test
                        bonding = self.mof[index].bondingatoms[name]
                    except:
                        # Ba2+ case is special
                        bonding = []

                    remaining_coords.append([self.coordinates[index][coord]
                        for coord in range(len(self.coordinates[index]))
                        if coord not in bonding and index != sbu])
                    remaining_atoms.append([self.atoms[index][atom]
                        for atom in range(len(self.atoms[index]))
                        if atom not in bonding and index != sbu])

            # make 1 dimensional list
            remaining_coords = [remaining_coords[i][j] for i in 
                range(len(remaining_coords)) for j in
                range(len(remaining_coords[i])) if remaining_coords[i][j]]
            remaining_atoms = [remaining_atoms[i][j] for i in 
                range(len(remaining_atoms)) for j in
                range(len(remaining_atoms[i])) if remaining_atoms[i][j]]
            distmat = distance.cdist(coords, remaining_coords)

            # list check will be populated if a distance falls below
            # the radii of two atoms multiplied by a scaling factor, sf
            check = [(i, j) for i in range(len(distmat)) for j in
                    range(len(distmat[i])) if distmat[i][j] <= 
                    (Radii[atoms[i]]+Radii[remaining_atoms[j]])*sf] 
            if len(check) > 0:
                info("Overlap found")
                return True
        return False

    def overlap_allcheck(self):
        """
        checks for all pairwise atomistic overlaps in a growing MOF
        """
        # tolerance in angstroms
        tol = 1.0 
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

        self.join_sbus(sbu1, bond1, sbu2, bond2, True)
        #self.xyz_debug()
        #dump = self.coordinate_dump()
        #write_xyz("history", dump[0], dump[1], self.cell, self.origins)
        if self.pbcindex == 2:
            self.fcoords.append([self.fractional(coord) for 
                coord in self.coordinates[sbu2]])

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

        while (np.allclose(axis, zeros1, atol=1e-4)):
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
        self.COM[sbu2] = vect_add(self.COM[sbu2],shiftvector)

    def getfinals(self):
        """
        get final coordinates from fractionals.
        """
        A = self.cell[0]
        B = self.cell[1]
        C = self.cell[2]
        finals = []
        for sbu in self.fcoords:
            carts = [matrx_mult(fcoord, self.cell) for fcoord
                    in self.fcoords[sbu]]
            finals = finals + carts
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
        self.COM[sbu] = vect_add(matrx_mult(self.COM[sbu], R), C)
        self.connect_align[sbu] = [matrx_mult(i, R) for i 
                                    in self.connect_align[sbu]]
    def get_cartesians(self):
        """
        convert fractional coordinates back to cartesian coordinates
        and store in self.coordinates.
        """
        for isbu in range(len(self.coordinates)):
            for icoord in range(len(self.coordinates[isbu])):
                coord = matrx_mult(self.fcoords[isbu][icoord], self.cell)
                self.coordinates[isbu][icoord] = coord[:]
            for iconnect in range(len(self.connect_points[isbu])):
                coord = matrx_mult(self.fconnect[isbu][iconnect],
                                   self.cell)
                self.connect_points[isbu][iconnect] = coord[:]

    def final_coords(self):
        """
        convert all coordinates to the shifted coordinates derived
        from the periodic boundary conditions
        """
        firstsbu = self.mof[0]
        cellcentre = zeros1[:]
        # re-initialize the fcoords
        # store connect_point fractionals
        self.fcoords = []
        self.fconnect = []
        for k in self.cell:
            cellcentre = vect_add(cellcentre, scalar_div(2., k))
        shiftvect = vect_sub(cellcentre, self.COM[0])
        for sbu in range(len(self.coordinates)):
            for ibond, bond in enumerate(self.connect_points[sbu]):
                self.connect_points[sbu][ibond] = vect_add(shiftvect, bond)
            for icoord, coord in enumerate(self.coordinates[sbu]):
                self.coordinates[sbu][icoord] = vect_add(shiftvect, coord)
            sbufracts = [self.fractional(coord) for coord in 
                            self.coordinates[sbu]]
            connectfracts = [self.fractional(coord) for coord in 
                             self.connect_points[sbu]]
            self.fcoords.append(sbufracts)
            self.fconnect.append(connectfracts)
        return


    def get_fractionals(self):
        """ Store fractional coordinates. """
        for sbu in range(len(self.coordinates)):
            sbufracts = [self.fractional(coord) for coord in 
                            self.coordinates[sbu]]
            self.fcoords.append(sbufracts)

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
            for sbu in range(len(self.fcoords)):
                self.fcoords[sbu] = [vect_sub([1.,1.,1.], i) for i
                        in self.fcoords[sbu]]
                self.fconnect[sbu] = [vect_sub([1.,1.,1.],i) for i
                        in self.fconnect[sbu]]
        # re-invert cell (is it necessary to invert before this?)
        self.invert_cell()
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

    def create_connect_table(self):
        """
        Determine atoms bonded between SBUs
        """
        # re-initialize table?
        self.master_table = {}
        for outter in range(len(self.connectivity)):
            for bond, inner in enumerate(self.connectivity[outter]):
                self.update_sbu_connections(outter, bond, inner)
    
        atomcount = 0
        # write out the bonding
        for sbu in range(len(self.connectivity)):
            type = self.storetype[sbu]
            for idx1, connect in enumerate(self.mof[sbu].table[type]):
                if (sbu, idx1) in self.sbu_connections.keys():
                    # sum the atoms up that come before the
                    # particular SBU
                    bond = self.sbu_connections[(sbu, idx1)]
                    bondedsbu = bond[0]
                    bondedatm = bond[1]
                    bondedlen = bond[2]
                    # Always assume a Single bond between SBU's
                    # TODO(pboyd): not robust and probably needs
                    # to be changed.
                    bondedtyp = "S"
                    atmsum = 0
                    if bondedsbu > 0:
                        atmsum = sum([len(i) for i in 
                                  self.coordinates[:bondedsbu]])
                    atmtup = [atomcount + idx1, atmsum + bondedatm]
                    atmtup.sort()
                    atmtup = tuple(atmtup)
                    self.master_table[atmtup] = (bondedlen, bondedtyp)

                for idx2, dist in enumerate(connect):
                    if dist != 0.:
                        atmtup = [atomcount + idx1, atomcount + idx2]
                        atmtup.sort()
                        atmtup = tuple(atmtup)
                        self.master_table[atmtup] = dist
            atomcount += len(self.mof[sbu].table[type])

    def update_sbu_connections(self, sbu1, bond1, sbu2):
        """
        Update the connectivity matrix of the MOF
        """
       
        # First check to see which atoms between SBUs are bonded
        # shift by pbc
        if self.pbcindex == 2:
            coords = self.min_img_shift(self.connect_points[sbu1][bond1])
        else:
            coords = self.coordinates[:]

        coord1 = coords[sbu1]
        coord2 = coords[sbu2]
        distmat = distance.cdist(coord1, coord2)

        bond = []
        # There is some leeway for bonding between SBUs so 
        # a scaling factor is introduced to ensure a bond 
        sf = 1.0
        bonddebug = {}
        ncount = 0
        while len(bond) == 0 and ncount < 7:
            for idx1 in range(len(distmat)):
                for idx2 in range(len(distmat[idx1])):
                    # Determine atom types
                    atom1 = self.atoms[sbu1][idx1]
                    atom2 = self.atoms[sbu2][idx2]
                    if distmat[idx1][idx2] < sf * bond_tolerance(
                        atom1, atom2):
                        if sbu1 == sbu2:
                            if idx1 != idx2:
                                # check to make sure the bonding isn't
                                # already recorded in the table
                                table = self.mof[sbu1].table[
                                        self.storetype[sbu1]]
                                if table[idx1][idx2] == 0. and \
                                set((atom1, atom2)) != set("C"):
                                    bonddebug.setdefault(idx1,
                                    []).append(idx2)
                                    bond.append((idx1, idx2,
                                        distmat[idx1][idx2]))
                        else:
                            bonddebug.setdefault(idx1, []).append(idx2)
                            bond.append((idx1, idx2, distmat[idx1][idx2]))
            # increase the scaling factor by 5% 
            sf *= 1.05
            # scale 10 times before giving up
            ncount += 1

        # specify these bonds in the inter-sbu connectivity table
        for i in bond:
            self.sbu_connections[(sbu1, i[0])] = (sbu2, i[1], i[2])
            self.sbu_connections[(sbu2, i[1])] = (sbu1, i[0], i[2])

        if len(bond) == 0:
            print "YIPES"

    def min_img_shift(self, vector, excl=None):
        """
        Shifts all coordinates of a MOF to within the minimum
        image of a vector.
        """
        fvect = matrx_mult(vector, self.icell)
        #fvect = [i - math.floor(i) for i in fvect]
        mof_coords = []
        for sbu in range(len(self.fcoords)):
            if sbu != excl:
                temp_fracts = []
                for fcoord in self.fcoords[sbu]:
                    diff = vect_sub(fvect, fcoord)
                    f = [fcoord[i] + round(diff[i]) 
                         for i in range(3)]
                    temp_fracts.append(f)
                mof_coords.append([matrx_mult(scaled, self.cell) for
                    scaled in temp_fracts])
        # convert temp_fracts to cartesians
        return mof_coords

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
            # store the third entry of the vector 
            # (this should not be altered)
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
            proj1 = project(vect, self.cell[0])
            proj2 = project(vect, self.cell[1])
            planproj = vect_add(proj1, proj2)
            angle = calc_angle(planproj, vect)
            if angle <= 0.4:
                info("Bad angle between vector and plane!")
                return False
        elif self.pbcindex == 2:
            # This may change in the future, if changeable 
            # periodic boundaries are considered.
            return False
        return True

    def add_vector(self, vect):
        self.pbcindex +=1
        self.cell[self.pbcindex] = vect
        self.ncell[self.pbcindex] = normalize(vect)
        if self.pbcindex == 2:
            # store all existing SBUs as fractional coordinates
            # Note: this may slow up the algorithm
            self.cellparams()
            self.invert_cell()
            self.get_fractionals()

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
        self.dbfile = file

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
        self.assign_internal_indices()

    def assign_internal_indices(self):
        """
        Gives each SBU a unique identifier (order index in list) which
        can be referenced outside of the database list
        """
        for ind, sbu in enumerate(self.database):
            sbu.internal_index = ind
        return

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

def rotation_axis(vect1, vect2):
    """
    axis of rotation for two vectors (90 degrees from both)
    """
    vect1 = normalize(vect1)
    vect2 = normalize(vect2)
    return normalize(cross(vect1, vect2))

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

    if len(sys.argv[1:]) > 0:
        dbfile = sys.argv[1]
    else:
        dbfile = "petedatabase"
    sbutest = Database(dbfile, options)
    genstruct = Generate(sbutest)
    genstruct.database_generation()

if __name__ == '__main__':
    main()
