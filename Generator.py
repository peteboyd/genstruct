#!/usr/bin/env python
import itertools
from random import choice
import sys

class Generate(object):
    """Takes as input a sequence of sbus, and returns
    build orders to make structures.
    
    """
    
    def __init__(self, options, sbu_list):
        self.options = options
        self.sbus = SBU_list(sbu_list)
        
    def generate_sbu_combinations(self, N=None):
        if N is None:
            N = self.options.metal_sbu_per_structure + self.options.organic_sbu_per_structure
        for i in itertools.combinations_with_replacement(self.sbus.list, N):
            if self._valid_sbu_combination(i):
                yield tuple(i)
    
    def combinations_from_options(self):
        """Just return the tuples in turn."""
        combs = []
        for combo in self.options.sbu_combinations:
            # first one has to be a metal.
            met = [self.sbus.get(combo[0], _METAL=True)]
            combs.append(tuple(met + [self.sbus.get(i) for i in combo[1:]]))
        return combs
    
    def _valid_sbu_combination(self, sbu_set):
        """Currently only checks if there is just one Metal and the rest organic."""
        return len([i for i in sbu_set if i.is_metal]) == 1
        
    def generate_build_directives(self, sbu, sbus):
        """Requires maximum length of sbu insertions."""

        # insert metal first
        if sbu is None:
            # chose a metal (at random, if more than one)
            sbu = choice([x for x in sbus if x.is_metal])
        
        # expand the current SBU's bonds and establish possible SBU bondings
        # generate exaustive list of sbu combinations.
        for k in self._yield_bonding_sbus(sbu, set(sbus), p=[0 for i in range(self.options.structure_sbu_length)]):
            # pair k to the appropriate sbu,
            # determine if the sbus can be bonded,
            # yield the build directive.
            sbu_bonding_pairs = self._gen_sbu_bonding_pairs([sbu] + self.flatten(k))
            if self._check_operation(sbu_bonding_pairs):
                for j in self._yield_valid_bonds(sbu_bonding_pairs):
                    yield [sbu] + zip(sbu_bonding_pairs,j)
                    # break loop if only the first set of bonding combinations
                    # are requested
                    if self.options.gen_single_bonding_set:
                        break
        
    def _check_operation(self, sbu_bonding_pairs):
        """Return True if the SBUs can be bonded together."""
        #print [(sbu.name, n_sbu.name) for sbu, n_sbu in sbu_bonding_pairs]
        for ind, sbu, n_sbu in sbu_bonding_pairs:
            if not self._valid_sbu_pair(sbu, n_sbu):
                return False
        return True
    
    def _valid_sbu_pair(self, sbu1, sbu2):
        """Determine if the two SBUs can be bonded.  Currently set to
        flag those which are of the same type (organic|metal)
        """
        if sbu1.is_metal != sbu2.is_metal:
            return True
        return False

    def _valid_sbu_pairs(self, sbu):
        """Hacky bit for the recursive SBU bond generator
        """
        def is_valid_tuple(sbu_tuple):
            return all([sbu.is_metal != sbut.is_metal for sbut in sbu_tuple])
        return is_valid_tuple

    def _gen_sbu_bonding_pairs(self, sbu_list):
        """Generate SBU-SBU pairs from the list of SBUs to be inserted.
        These pairs are generated on a first-come-first-serve basis
        depending on how many connect_points belong to the SBU.
        
        """
        
        pair = []
        sbu_cp = 0
        for ind, sbu in enumerate(sbu_list):
            for n_sbu in itertools.islice(sbu_list, sbu_cp + 1, sbu_cp + len(sbu.connect_points) + 1):
                pair.append((ind, sbu, n_sbu))
            sbu_cp += len(sbu.connect_points)                    
        return pair
    
    def flatten(self, s):
        """Returns a flattened list"""
        if s == []:
            return s
        if isinstance(s[0], list):
            return self.flatten(s[0]) + self.flatten(s[1:])
        return s[:1] + self.flatten(s[1:])
    
    def _yield_bonding_sbus(self, sbu, sbus, index=0, p=[]):
        """Return a tuple of SBUs generated exhaustively.
        """
        if index == self.options.structure_sbu_length:
            yield p
        else:
            index += 1
            func = self._valid_sbu_pairs(sbu)
            for iterator in itertools.ifilter(func,
                                            itertools.product(sbus, repeat=(len(sbu.connect_points)))):
                p[index-1] = list(iterator)
                q = self.flatten(p)[index-1]
                for s in self._yield_bonding_sbus(q, sbus, index, p):
                    yield s
        
    def _yield_valid_bonds(self, sbu_pairs):
        """Return the valid connect points between SBUs in the sbu_order"""

        for x in itertools.product(*[itertools.product(*[itertools.product([i], sbu2.connect_points)
                                                        for i in sbu1.connect_points])
                                                        for ind, sbu1, sbu2 in sorted(list(set(sbu_pairs)))]):
            if all([y.special == z.special and not y == z for sbu in x for y, z in sbu]):
                    yield [directive for sbu in x for directive in sbu]
        
    def build_directives_from_options(self, directive):
        """Returns the build directives from the options (if it exists)"""
        return None
        
class SBU_list(object):
    
    def __init__(self, sbu_list):
        self.list = sbu_list
    
    def get(self, identifier, _METAL=False):
        """Produces the SBU with the identifier provided, this filters between
        organic and metal SBUs"""
        for sbu in self.list:
            if sbu.identifier == identifier and sbu.is_metal == _METAL:
                return sbu
        raise Exception("Could not find the SBU with the identifier %s"%(identifier))
    
    def getall(self, identifier):
        """Produces the SBU with the target identifier regardless of being
        metal or organic."""
        for sbu in self.list:
            if sbu.identifier == identifier:
                return sbu
        raise Exception("Could not find the SBU with the identifier %s"%(identifier))

