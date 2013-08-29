#!/usr/bin/env python
import numpy as np

class Structure(object):
    """Structure class contains atom info for MOF."""
    
    def __init__(self, name=None):
        self.name = name
        self.cell = Cell()
        self.atoms = []
        self.bonds = []
        # re-orient the cell
        # shift all to within the periodic boundaries
        # run symmetry finding on it.
        # write to cif file.

    def from_build(self, build_obj):
        """Build structure up from the builder object"""
        # sort out the connectivity information
        # copy over all the Atoms
        
        for sbu in build_obj.sbus:

    

    def check_atomistic_overlap(self):
        """Determine if there is atom overlap, taking into account
        bonding if bonding information is provided, as well as 
        periodic boundary conditions."""
        pass

    def detect_symmetry(self):
        pass

    def re_orient_cell(self):
        """Adjusts cell vectors to lie in the standard directions."""
        pass

    def write_to_cif(self):
        """Write structure information to a cif file."""
        pass


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
   
    def _re_orient(self):


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

    def __mkparam(self):
        """Update the parameters to match the cell."""
        self.params[0:3] = [np.linalg.norm(i) for i in self.lattice][:]
       
 
        alpha = arccos(sum(self.cell[1, :] * self.cell[2, :]) /
                       (cell_b * cell_c)) * 180 / pi
        beta = arccos(sum(self.cell[0, :] * self.cell[2, :]) /
                      (cell_a * cell_c)) * 180 / pi
        gamma = arccos(sum(self.cell[0, :] * self.cell[1, :]) /
                       (cell_a * cell_b)) * 180 / pi
        self._params = (cell_a, cell_b, cell_c, alpha, beta, gamma)
 
    @property
    def a(self):
        """Magnitude of cell a vector."""
        return self.params[0]

    @property
    def b(self):
        """Magnitude of cell b vector."""
        return self.params[1]

    @property
    def c(self):
        """Magnitude of cell c vector."""
        return self.params[2]

    @property
    def alpha(self):
        """Cell angle alpha."""
        return self.params[3]

    @property
    def beta(self):
        """Cell angle beta."""
        return self.params[4]

    @property
    def gamma(self):
        """Cell angle gamma."""
        return self.params[5]

