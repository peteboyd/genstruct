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
        """Currently only checks if there is the correct number of metal 
        SBUs in the combination."""
        return len([i for i in sbu_set if i.is_metal]) == \
                self.options.metal_sbu_per_structure

        
    def generate_build_directives(self, sbu, sbus):
        """Requires maximum length of sbu insertions."""

        # insert metal first
        if sbu is None:
            # chose a metal (at random, if more than one)
            sbu = choice([x for x in sbus if x.is_metal])
        
        # expand the current SBU's bonds and establish possible SBU bondings
        # generate exaustive list of sbu combinations.
        for k in self._yield_bonding_sbus(sbu, set(sbus), 
                p=[0 for i in range(self.options.structure_sbu_length)]):
            # pair k to the appropriate sbu,
            # determine if the sbus can be bonded,
            # yield the build directive.
            print k
            sbu_bonding_pairs = self._gen_sbu_bonding_pairs([sbu] + self.flatten(k))
            if self._check_operation(sbu_bonding_pairs):
                for j in self._yield_valid_bonds(sbu_bonding_pairs):
                    yield [sbu] + zip(sbu_bonding_pairs,j)
        
    def _check_operation(self, sbu_bonding_pairs):
        """Return True if the SBUs can be bonded together."""
        #print [(sbu.name, n_sbu.name) for sbu, n_sbu in sbu_bonding_pairs]
        for ind, sbu, n_sbu in sbu_bonding_pairs:
            if not self._valid_sbu_pair(sbu, n_sbu):
                return False
        return True
    
    def _valid_sbu_pair(self, sbu1, sbu2):
        """Determine if the two SBUs can be bonded.  Currently set to
        flag true if the two sbus contain matching bond flags, otherwise
        if they are a (metal|organic) pair
        """
        spec_cp1s = [cp.special for cp in sbu1.connect_points if cp.special]
        const_cp1s = [cp.constraint for cp in sbu1.connect_points if cp.constraint]
        spec_cp2s = [cp.special for cp in sbu2.connect_points if cp.special]
        const_cp2s = [cp.constraint for cp in sbu2.connect_points if cp.special]

        none_spec_cp1 = [cp.special is None for cp in sbu1.connect_points]
        none_const_cp1 = [cp.constraint is None for cp in sbu1.connect_points]
        none_spec_cp2 = [cp.special is None for cp in sbu2.connect_points]
        none_const_cp2 = [cp.constraint is None for cp in sbu2.connect_points]
        # Return true if they have matching bond_flags
        if any([i == j for i in spec_cp1s for j in const_cp2s]) or \
                any([i == j for i in spec_cp2s for j in const_cp1s]):
            return True

        if any(none_spec_cp1) and any(none_const_cp1) and any(none_spec_cp2) \
                and any(none_const_cp2):    
            if sbu1.is_metal != sbu2.is_metal:
                return True
        return False

    def _valid_sbu_pairs(self, sbu):
        """Hacky bit for the recursive SBU bond generator
        """
        def is_valid_tuple(sbu_tuple):
            # check if one of the sbus contains a special bond equal to 
            # sbus' special bond.
            return all([self._valid_sbu_pair(sbu, sbut) for sbut in sbu_tuple])
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

    def roundrobin(self, *iterables):
        pending = len(iterables)
        nexts = itertools.cycle(iter(it).next for it in iterables)
        while pending:
            try:
                for next in nexts:
                    yield next()
            except StopIteration:
                pending -= 1
                nexts = itertools.cycle(itertools.islice(nexts, pending))

    def _gen_bonding_sbus(self, sbu, sbus, index=0):
        """Returns an iterator which runs over tuples of bonds
        with other sbus."""
        # an iterator that iterates sbu's first, then sbus' connect_points
        ncps = len(sbu.connect_points)
        sbu_repr = list(itertools.product([sbu], sbu.connect_points))

        bond_iter = self.roundrobin(*[itertools.product([s], s.connect_points)
                                    for s in sbus])
        # don't like how this iterates, but will do for now.
        all_bonds = itertools.tee(bond_iter, ncps)
        for bond_set in itertools.product(*all_bonds):
            full_bond_set = list(itertools.izip(sbu_repr, bond_set))
            if all([self._valid_bond_pair(i) for i in full_bond_set]):
                yield [((index, cp1.identifier),(sbu2, cp2)) for 
                        (sbu1, cp1),(sbu2,cp2) in full_bond_set]

    def _valid_bond_pair(self, set):
        """Determine if the two SBUs can be bonded.  Currently set to
        flag true if the two sbus contain matching bond flags, otherwise
        if they are a (metal|organic) pair
        """
        (sbu1, cp1), (sbu2, cp2) = set
        if all([i is None for i in [cp1.special, cp2.special, cp1.constraint, cp2.constraint]]):
            return sbu1.is_metal != sbu2.is_metal

        return (cp1.special == cp2.constraint) and (cp2.special == cp1.constraint)

    def _yield_bonding_sbus(self, sbu, sbus, index=0, p=[]):
        """Return a tuple of SBUs generated exhaustively.
        """
        if index == self.options.structure_sbu_length:
            yield p
        else:
            index += 1
            func = self._valid_sbu_pairs(sbu)
            #for iterator in itertools.ifilter(func,
            #                                itertools.product(sbus, repeat=(len(sbu.connect_points)))):       

            #TODO(pboyd): Probably ignore bonding with the metal-metal cases, since they will likely always form a periodic boundary right at the beginning of the Build.
            for iterator in self._gen_bonding_sbus(sbu, sbus, index):

                p[index-1] = list(iterator)
                q = self.flatten(p)[index-1][1][0]
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
        self._truncate()

    def _truncate(self):
        trunc = []
        for sbu1, sbu2 in itertools.combinations(self.list, 2):
            if sbu1.parent == sbu2.name:
                trunc.append(self.list.index(sbu1)) 
                sbu2.children.append(sbu1)
            elif sbu2.parent == sbu1.name:
                trunc.append(self.list.index(sbu2))
                sbu1.children.append(sbu2)
        for k in reversed(sorted(trunc)):
            del self.list[k] 

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

