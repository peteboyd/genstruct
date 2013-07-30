#!/usr/bin/env python
import numpy as np

class Structure(object):
    """Structure class contains atom info for MOF."""
    
    def __init__(self, name=None):
        self.name = name
        self.cell = Cell()
        self.atoms = []
        self.bonds = []
        
        
class Cell(object):
    """contains periodic vectors for the structure."""
    
    def __init__(self):
        self.basis = 0
        self.lattice = np.identity(3)
        self.nlattice = np.zeros((3,3))
        self.params = np.zeros(6)
        
        
    @property
    def inverse(self):
        try:
            return self._ilattice
        except AttributeError:
            self._ilattice = np.array(np.matrix(self.lattice).I)
            return self._ilattice
    
    def add(self, index, vector):
        """Adds a periodic vector to the lattice."""
        self.lattice[index][:] = vector.copy()
        self.nlattice[index][:] = vector.copy() / np.linalg.norm(vector)
        
    def to_xyz(self):
        """Returns a list of the strings"""
        lines = []
        for vector in self.lattice:
            lines.append("atom_vector %12.5f %12.5f %12.5f\n"%(tuple(vector)))
                         
        return lines
        
        