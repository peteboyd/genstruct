#!/usr/bin/env python
import numpy as np
from Structure import Structure, Cell
from SecondaryBuildingUnit import SBU
import sys
import os
import copy
import math
from scipy.spatial import distance
import itertools
from logging import info, debug, warning, error, critical
from elements import Radii

class Build(object):
    """Builds MOFs given directives."""
    def __init__(self, options):
        self.options = options
        # store SBU copies
        self.sbus = []
        self.structure = Structure()
        self.periodic_vectors = self.structure.cell
        self.periodic_origins = np.zeros((3,3))
        self.periodic_index = 0
        
    def build_from_directives(self, directives):
        self.sbus = []
        index_type = []
        if self.options.debug_writing:
            self.init_debug()

        for count, operation in enumerate(directives):

            if isinstance(operation, SBU):
                # starting seed
                index_type.append(0)
                self.sbus.append(copy.deepcopy(operation))
                if self.options.debug_writing:
                    self.debug_xyz(operation)            

            elif operation[0][0] not in index_type:
                # this means the SBU wasn't successfully inserted, so all building
                # commands to this SBU should be ignored
                pass
            elif self.sbus[index_type.index(operation[0][0])].\
                        get_cp(operation[1][0].identifier).connected:
                debug("SBU %i already bonded at connect point %i"%(
                            index_type.index(operation[0][0]),operation[1][0].identifier))
            else:
                # sbu1 and connect_point1 must be related to a copied SBU
                # in self.sbus
                debug("Length of structure = %i"%len(self.sbus))
                sbu_ind1 = index_type.index(operation[0][0])
                sbu1 = self.sbus[sbu_ind1]
                connect_point1 = sbu1.get_cp(operation[1][0].identifier)
                sbu2 = copy.deepcopy(operation[0][2])
                connect_point2 = sbu2.get_cp(operation[1][1].identifier)
                debug("Trying SBU %s on (SBU %i, %s) "%(
                    sbu2.name, sbu_ind1, sbu1.name) + "using the bonds (%i, %i)"%(
                    connect_point2.identifier, connect_point1.identifier))
                if self.options.debug_writing:
                    self.debug_xyz(sbu2)
                # perform transformations to sbu2               
                self.translation(sbu1, connect_point1,
                                 sbu2, connect_point2)
                
                if self.options.debug_writing:
                    self.debug_xyz(sbu2)
                    
                self.rotation_z(sbu1, -connect_point1,
                                sbu2, connect_point2)
                
                if self.options.debug_writing:
                    self.debug_xyz(sbu2)                
                self.rotation_y(sbu1, connect_point1,
                                sbu2, connect_point2)

                if self.options.debug_writing:
                    self.debug_xyz(sbu2)
                
                # overlap check
                if self.overlap(sbu2):
                    debug("overlap found")
                    return
                else:
                    # check for periodic boundaries
                    index_type.append(count)
                    self.sbus.append(sbu2)
                    connect_point1.connected = True
                    connect_point1.sbu_bond = (len(self.sbus)-1, connect_point2.identifier)
                    connect_point2.connected = True
                    connect_point2.sbu_bond = (sbu_ind1, connect_point1.identifier)
                    self.bonding_check()
                    if self._completed_structure():
                        # test for periodic overlaps.
                        info("Structure Generated!")
                        self.structure.from_build(self.sbus)
                        return
        
    def bonding_check(self):
        """Evaluate the presence of bonds between existing SBUs"""
        bond_points = [(ind,cp) for ind, sbu in enumerate(self.sbus)
                       for cp in sbu.connect_points]
        for (ind1, cp1), (ind2, cp2) in itertools.combinations(bond_points, 2):
            if self._valid_bond(ind1, cp1, ind2, cp2):
                distance_vector = cp1.origin[:3] - cp2.origin[:3]
                if self.periodic_index == 3:
                    # shift the vector by periodic boundaries
                    distance_vector = self.periodic_shift(distance_vector)
                if np.linalg.norm(distance_vector) < self.options.distance_tolerance:
                    # local bond
                    debug("Bond found between %s, and %s (%i,%i) at bonds (%i,%i)"%(
                        self.sbus[ind1].name, self.sbus[ind2].name, ind1, ind2,
                        cp1.identifier, cp2.identifier))
                    cp1.connected, cp2.connected = True, True
                    cp1.sbu_bond = (ind2, cp2.identifier)
                    cp2.sbu_bond = (ind1, cp1.identifier)
                    
                elif self._valid_periodic_vector(distance_vector):
                    # new periodic boundary
                    debug("Periodic boundary found (%5.3f, %5.3f, %5.3f)"%(
                        tuple(distance_vector)))
                    self.periodic_vectors.add(self.periodic_index, distance_vector)
                    self.periodic_origins[self.periodic_index][:] = cp2.origin[:3].copy()
                    self.periodic_index += 1
                    cp1.connected, cp2.connected = True, True
                    cp1.sbu_bond = (ind2, cp2.identifier)
                    cp2.sbu_bond = (ind1, cp1.identifier)           
    
    def _completed_structure(self):
        return (self.periodic_index == 3 and
                all([cp.connected for sbu in self.sbus for cp in sbu.connect_points]))
                   
    def periodic_shift(self, vector):
        proj_vect = np.dot(vector, self.periodic_vectors.inverse)
        proj_vect = np.rint(proj_vect)
        shift_vector = np.dot(proj_vect, self.periodic_vectors.lattice)
        return (vector - shift_vector)
    
    def init_debug(self):
        write = {'append':'a', 'overwrite':'w', 'write':'w', 'w':'w', 'a':'a'}
        assert self.options.debug_writing in write.keys()
        filename = os.path.join(self.options.job_dir,
                                self.options.jobname + ".debug.xyz")        
        f = open(filename, write[self.options.debug_writing])
        f.close()
        return filename
              
    def debug_xyz(self, sbu):
        filename = os.path.join(self.options.job_dir,
                                self.options.jobname + ".debug.xyz")
        filestream = open(filename, 'a')
        # determine the number of pseudo atoms and atoms
        natoms = (self.periodic_index + sum([len(q.connect_points) for q in self.sbus])
                  + len(sbu.connect_points))
        natoms += sum([len(bb.atoms) for bb in self.sbus + [sbu]])
        lines = "%i\ndebug_file\n"%(natoms)
        for pv in range(self.periodic_index):
            lines += "Na %12.5f %12.5f %12.5f "%(tuple(self.periodic_origins[pv]))
            lines += self.periodic_vectors.to_xyz()[pv]

        for mol in self.sbus + [sbu]:
            lines += str(mol)
        filestream.writelines(lines)
        filestream.close()
    
    def overlap(self, sbu):
        """Just perform local atom-atom overlap distance checks.
        The periodic boundary checks will commence once the structure is
        generated.
        """
        if not self.options.overlap_tolerance:
            return False
        # just ignore the atoms which will likely be in close contact
        # with other SBUs due to bonding.
        coords1 = np.array([atom.coordinates for atom in sbu.atoms
                            if not atom.sbu_bridge])
        for o_sbu in self.sbus:
            coords2 = np.array([atom.coordinates for atom in o_sbu.atoms
                                if not atom.sbu_bridge])
            dist_mat = distance.cdist(coords1, coords2)
            for (atom1, atom2), dist in np.ndenumerate(dist_mat):
                elem1, elem2 = sbu.atoms[atom1].element, o_sbu.atoms[atom2].element
                if (Radii[elem1] + Radii[elem2]) * self.options.overlap_tolerance > dist:
                    return True
        return False
    
    def _valid_periodic_vector(self, vector):
        if self.periodic_index >= 3:
            return False       
        if self.periodic_index == 0:
            return True
        nvec = vector/np.linalg.norm(vector)
        # check if the vector is a linear combination
        # of the existing vectors
        for latt_vec in self.periodic_vectors.nlattice[:self.periodic_index]:
            if np.allclose(np.dot(latt_vec, nvec), 1., atol=0.1):
                return False
        if self.periodic_index == 2:
            # check for co-planarity
            norm_test = np.dot(nvec, np.cross(self.periodic_vectors.nlattice[0],
                                     self.periodic_vectors.nlattice[1]))
            if np.allclose(norm_test, 0., atol=0.1):
                return False
        return True

        
    def _valid_bond(self, ind1, cp1, ind2, cp2):
        """Check if a bond can be formed between two SBUs"""
        sbu1 = self.sbus[ind1]
        sbu2 = self.sbus[ind2]
        # check if either is bonded already
        if cp1.connected or cp2.connected:
            return False
        # check if the two sbus are the same and the cp's are the same
        check_same = (sbu1.is_metal == sbu2.is_metal,
                      cp1.special == None, cp2.special == None)
        if all(check_same):
            return False
        check_same = (sbu1.is_metal == sbu2.is_metal, cp1.special == cp2.special,
                      cp1.identifier == cp2.identifier)
        if all(check_same):
            return False
        # return false if either connect point has a special flag.
        if cp1.special != cp2.special:
            return False
        
        # check if vectors are aligned
        if not np.allclose(np.dot(cp1.z[:3], cp2.z[:3]), -1., atol=0.05):
            return False
        # (anti)parallel alignment vectors
        if not np.allclose(abs(np.dot(cp2.y[:3], cp2.y[:3])), 1., atol=0.05):
            return False
        return True
        
        
    def translation(self, sbu1, cp1, sbu2, cp2):
        vect = cp1.origin[:3] - cp2.origin[:3]
        sbu2.translate(vect)
        
    def rotation_z(self, sbu1, cp1, sbu2, cp2):
        # first rotation
        cp = np.cross(cp1.z[:3], cp2.z[:3])
        axis = cp/np.linalg.norm(cp)
        angle = self.calc_angle(cp1.z, cp2.z)
        R = self.rotation_matrix(axis, angle, point=cp2.origin[:3])
        sbu2.rotate(R)
        
    def rotation_y(self, sbu1, cp1, sbu2, cp2):
        # second
        cp = np.cross(cp1.y[:3], cp2.y[:3])
        axis = cp/np.linalg.norm(cp)
        angle = self.calc_angle(cp1.y, cp2.y)
        R = self.rotation_matrix(axis, angle, point=cp2.origin[:3])
        sbu2.rotate(R)

    def calc_angle(self, vect1, vect2):
        ic1 = vect1[:3] / np.linalg.norm(vect1[:3])
        ic2 = vect2[:3] / np.linalg.norm(vect2[:3])
        a = min(max(np.dot(ic1, ic2), -1.0), 1.0)
        return math.acos(a)

    def rotation_matrix(self, axis, angle, point=None):
        """
        returns a 3x3 rotation matrix based on the
        provided axis and angle
        """
        axis = np.array(axis)
        axis = axis / np.linalg.norm(axis)
        a = np.cos(angle / 2.)
        b, c, d = -axis*np.sin(angle / 2.)

        R = np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
                  [2*(b*c + a*d), a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
                  [2*(b*d - a*c), 2*(c*d + a*b), a*a + d*d - b*b - c*c]])

        M = np.identity(4)
        M[:3,:3] = R
        if point is not None:
            # rotation not around origin
            point = np.array(point[:3], dtype=np.float64, copy=False)
            M[:3,3] = point - np.dot(R, point)
        return M

