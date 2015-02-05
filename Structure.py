#!/usr/bin/env python
import numpy as np
import itertools
from scipy.spatial import distance
from LinAlg import LinAlg, RAD2DEG
from element_properties import Radii, ATOMIC_NUMBER
from CIFer import CIF
import os
import sys

class Structure(object):
    """Structure class contains atom info for MOF."""
    
    def __init__(self, options, name=None):
        self.name = name
        self.options = options
        self.cell = Cell()
        self.atoms = []
        self.bonds = {}
        self.fragments = [] 
        self.build_directives = None
        self.charge = 0
        self.space_group_name = "P1"
        self.space_group_number = 1
        self.symmetry_operations = ['x, y, z']
        self.cell_setting = 'monoclinic'

    def from_build(self, build_obj):
        """Build structure up from the builder object"""
        # sort out the connectivity information
        # copy over all the Atoms
        self.cell = build_obj.periodic_vectors
        index_count = 0
        for order, sbu in enumerate(build_obj.sbus):
            sbu.update_atoms(index_count, order)
            self.charge += sbu.charge
            self.fragments.append((sbu.name, order))
            self.atoms += sbu.atoms
            if any([i in self.bonds.keys() for i in sbu.bonds.keys()]):
                warning("Two bonds with the same indices found when forming"+
                        " the bonding table for the structure!")
            self.bonds.update(sbu.bonds)
            index_count += len(sbu)

        for id, sbu in enumerate(build_obj.sbus):
            for cp in sbu.connect_points:
                sbu2 = build_obj.sbus[cp.sbu_bond[0]]
                cp2_id = cp.sbu_bond[1]
                self.compute_inter_sbu_bonding(sbu, cp.identifier, sbu2, cp2_id)

    def compute_inter_sbu_bonding(self, sbu1, cp1_id, sbu2, cp2_id):
        # find out which atoms are involved in inter-sbu bonding
        atoms1 = [i for i in sbu1.atoms if cp1_id in i.sbu_bridge]
        atoms2 = [j for j in sbu2.atoms if cp2_id in j.sbu_bridge]
        # measure minimum distances between atoms to get the 
        # correct bonding.
        base_atoms = atoms1 if len(atoms1) >= len(atoms2) else atoms2
        bond_atoms = atoms2 if len(atoms2) <= len(atoms1) else atoms1
        for atom in base_atoms:
            if bond_atoms:
                shifted_coords = self.min_img(atom, bond_atoms)
                dist = distance.cdist([atom.coordinates[:3]], shifted_coords)
                dist = dist[0].tolist()
                bond_atom = bond_atoms[dist.index(min(dist))]
                if tuple([atom.index, bond_atom.index]) in self.bonds.keys():
                    bond = tuple([bond_atom.index, atom.index])
                elif tuple([bond_atom.index, atom.index]) in self.bonds.keys():
                    bond = tuple([atom.index, bond_atom.index])
                else:
                    bond = tuple(sorted((atom.index, bond_atom.index)))

                self.bonds.update({bond: "S"})

    def _compute_bond_info(self):
        """Update bonds to contain bond type, distances, and min img
        shift."""
        supercells = np.array(list(itertools.product((-1, 0, 1), repeat=3)))
        unit_repr = np.array([5,5,5], dtype=int)
        for (at1, at2), val in self.bonds.items():
            atom1 = self.atoms[at1]
            atom2 = self.atoms[at2]
            fcoords = atom2.scaled_pos(self.cell.inverse) + supercells
            coords = []
            for j in fcoords:
                coords.append(np.dot(j, self.cell.lattice))
            coords = np.array(coords)
            dists = distance.cdist([atom1.coordinates[:3]], coords)
            dists = dists[0].tolist()
            image = dists.index(min(dists))
            dist = min(dists)
            sym = '.' if all([i==0 for i in supercells[image]]) else \
                    "1_%i%i%i"%(tuple(np.array(supercells[image],dtype=int) +
                                      unit_repr))
            self.bonds[(at1, at2)] = (val, dist, sym) 

    def compute_overlap(self):
        """Determines if there is atomistic overlap. Includes periodic
        boundary considerations."""
        for id, atom in enumerate(self.atoms):
            elem1 = atom.element
            non_bonded = [i for i in self.atoms[id:] if 
                    tuple(sorted((atom.index, i.index))) not in self.bonds.keys()]
            indices = [i.index for i in non_bonded] 
            shifted_vectors = self.min_img(atom, non_bonded)
            dist_mat = distance.cdist([atom.coordinates[:3]], shifted_vectors)
            for (atom1, atom2), dist in np.ndenumerate(dist_mat):
                id2 = indices[atom2]
                elem2 = self.atoms[id2].element
                if (Radii[elem1] + Radii[elem2])*self.options.overlap_tolerance > dist:
                    return True
        return False

    def min_img(self, atom, atoms):
        """Orient all atoms to within the minimum image 
        of the provided atom."""
        sc_atom = atom.scaled_pos(self.cell.inverse)
        shifted_coords = []
        for at in atoms:
            scaled = at.scaled_pos(self.cell.inverse)
            shift = np.around(sc_atom - scaled)
            shifted_coords.append(np.dot((scaled+shift), self.cell.lattice))
        return shifted_coords

    def re_orient(self):
        """Adjusts cell vectors to lie in the standard directions."""
        frac_coords = [i.in_cell_scaled(self.cell.inverse) 
                        for i in self.atoms]
        self.cell.reorient_lattice()
        for id, atom in enumerate(self.atoms):
            atom.coordinates[:3] = np.dot(frac_coords[id],
                                          self.cell.lattice)
    
    def compute_symmetry(self):
        """Convert the P1 structure into it's symmetry irreduciple representation"""
        sym = Symmetry(self.options)
        sym.add_structure(self)
        sym.refine_cell()
        self.space_group_name = sym.get_space_group_name()
        self.space_group_number = sym.get_space_group_number()
        self.symmetry_operations = sym.get_space_group_operations()
        self.cell_setting = sym.cell_setting[self.space_group_number]
        symmetry_equiv_atoms = sym.get_equiv_atoms()

        self.cell.lattice = sym._lattice
        self.cell._ilattice = np.array(np.matrix(sym._lattice).I)
        syms = {}
        bonding = {}
        equivs = sym.get_equiv_atoms()
        for ind, i in enumerate(equivs):
            syms.setdefault(i, []).append(ind)
        #ats = [syms[0][0]]
        ats = syms.keys()
        ats_equivs = [0]
        symbonds = {}
        #while len(ats) < len(syms.keys()):
        #    for i in ats:
        #        bonds = [j for j in self.bonds.keys() if i in j]
        #        symbonds.update({j:self.bonds[j] for j in bonds})
        #        neighbours = [j[0] for j in bonds if j[0] != i ]
        #        neighbours += [j[1] for j in bonds if j[1] != i ]
        #        for j in neighbours:
        #            ksym = [k for k in syms.keys() if j in syms[k]][0]
        #            
        #            if ksym in ats_equivs:
        #                pass
        #            else:
        #                ats_equivs.append(ksym)
        #                ats.append(j)
        #                if len(ats) >= len(syms.keys()):
        #                    break
        for i in ats:
            bonds = [j for j in self.bonds.keys() if i in j]
            symbonds.update({j:self.bonds[j] for j in bonds})

        atom_equivs = {}
        for atat in ats:
            for j in syms.keys():
                if atat in syms[j]:
                    break
            for atatat in syms[j]:
                atom_equivs[atatat] = atat


        for j in symbonds.keys():
            i0 = atom_equivs[j[0]]
            i1 = atom_equivs[j[1]]
            
            n1 = ats.index(i0)
            n2 = ats.index(i1)
            t, d, p = symbonds[j]
            symbonds[j] = (n1, n2, t, d, p)
            
        self.bonds = symbonds
         
        atoms = [self.atoms[i] for i in ats]
        for ind, atom in enumerate(atoms):

            c = sym._scaled_coords[ats[ind]]
            atom.coordinates = np.dot(c, self.cell.lattice)

        #for j in self.bonds.keys():
        #    i0 = ats.index(j[0])
        #    i1 = ats.index(j[1])
        #    val = self.bonds.pop(j)
        #    self.bonds[(i0, i1)] = val
        self.cell.reparam()
        self.atoms = atoms

    def write_cif(self):
        """Write structure information to a cif file."""
        self._compute_bond_info()
        if self.options.find_symmetry:
            self.compute_symmetry()

        c = CIF(name=self.name)
        if self.space_group_name == "P1":
            c.insert_block_order("fragment", 4)

        labels = []
        # data block
        c.add_data("data", data_=self.name)
        c.add_data("data", _audit_creation_date=
                            CIF.label(c.get_time()))
        c.add_data("data", _audit_creation_method=
                            CIF.label("Genstruct v.%4.3f"%(
                                    self.options.version)))
        if self.charge:
            c.add_data("data", _chemical_properties_physical=
                               "net charge is %i"%(self.charge))

        # sym block
        c.add_data("sym", _symmetry_space_group_name_H_M=
                            CIF.label(self.space_group_name))
        c.add_data("sym", _symmetry_Int_Tables_number=
                            CIF.label(str(self.space_group_number)))
        c.add_data("sym", _symmetry_cell_setting=
                            CIF.label(self.cell_setting))

        # sym loop block
        for i in self.symmetry_operations: 
            c.add_data("sym_loop", _symmetry_equiv_pos_as_xyz=
                                   CIF.label("'%s'"%i))

        # cell block
        c.add_data("cell", _cell_length_a=CIF.cell_length_a(self.cell.a))
        c.add_data("cell", _cell_length_b=CIF.cell_length_b(self.cell.b))
        c.add_data("cell", _cell_length_c=CIF.cell_length_c(self.cell.c))
        c.add_data("cell", _cell_angle_alpha=CIF.cell_angle_alpha(self.cell.alpha))
        c.add_data("cell", _cell_angle_beta=CIF.cell_angle_beta(self.cell.beta))
        c.add_data("cell", _cell_angle_gamma=CIF.cell_angle_gamma(self.cell.gamma))

        if self.space_group_name == "P1":
            for name, order in self.fragments:
                c.add_data("fragment", _chemical_identifier=CIF.label(order),
                                   _chemical_name=CIF.label(name))
        # atom block
        element_counter = {}
        for atom in self.atoms:
            label = c.get_element_label(atom.element)
            labels.append(label)
            c.add_data("atoms", _atom_site_label=
                                    CIF.atom_site_label(label))
            c.add_data("atoms", _atom_site_type_symbol=
                                    CIF.atom_site_type_symbol(atom.element))
            c.add_data("atoms", _atom_site_description=
                                    CIF.atom_site_description(atom.force_field_type))
            if self.space_group_name == "P1":
                c.add_data("atoms", _atom_site_fragment=CIF.atom_site_fragment(atom.sbu_order))
            fc = atom.scaled_pos(self.cell.inverse)
            c.add_data("atoms", _atom_site_fract_x=
                                    CIF.atom_site_fract_x(fc[0]))
            c.add_data("atoms", _atom_site_fract_y=
                                    CIF.atom_site_fract_y(fc[1]))
            c.add_data("atoms", _atom_site_fract_z=
                                    CIF.atom_site_fract_z(fc[2]))

        # bond block
        if self.options.find_symmetry:
            for (dump1, dump2), (at1, at2, type, dist, sym) in self.bonds.items():
                label1 = labels[at1]
                label2 = labels[at2]
                c.add_data("bonds", _geom_bond_atom_site_label_1=
                                            CIF.geom_bond_atom_site_label_1(label1))
                c.add_data("bonds", _geom_bond_atom_site_label_2=
                                            CIF.geom_bond_atom_site_label_2(label2))
                c.add_data("bonds", _geom_bond_distance=
                                            CIF.geom_bond_distance(dist))
                c.add_data("bonds", _geom_bond_site_symmetry_2=
                                            CIF.geom_bond_site_symmetry_2(sym))
                c.add_data("bonds", _ccdc_geom_bond_type=
                                            CIF.ccdc_geom_bond_type(type))
        else:
            for (at1, at2), (type, dist, sym) in self.bonds.items():
                label1 = labels[at1]
                label2 = labels[at2]
                c.add_data("bonds", _geom_bond_atom_site_label_1=
                                            CIF.geom_bond_atom_site_label_1(label1))
                c.add_data("bonds", _geom_bond_atom_site_label_2=
                                            CIF.geom_bond_atom_site_label_2(label2))
                c.add_data("bonds", _geom_bond_distance=
                                            CIF.geom_bond_distance(dist))
                c.add_data("bonds", _geom_bond_site_symmetry_2=
                                            CIF.geom_bond_site_symmetry_2(sym))
                c.add_data("bonds", _ccdc_geom_bond_type=
                                            CIF.ccdc_geom_bond_type(type))

        file = open("%s.cif"%self.name, "w")
        file.writelines(str(c))
        file.close()

class Cell(object):
    """contains periodic vectors for the structure."""
    
    def __init__(self):
        self.basis = 0
        self.lattice = np.identity(3)
        self.nlattice = np.zeros((3,3))
        
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
       
    def remove(self, index):
        if index == 0:
            self.lattice[:2] = self.lattice[1:]
        elif index == 1:
            self.lattice[1] = self.lattice[2]
        self.lattice[2] = np.identity(3)[2]

    def to_xyz(self):
        """Returns a list of the strings"""
        lines = []
        for vector in self.lattice:
            lines.append("atom_vector %12.5f %12.5f %12.5f\n"%(tuple(vector)))
                         
        return lines

    def reparam(self):
        self.__mkparam()

    def __mkparam(self):
        """Update the parameters to match the cell."""
        self._params = np.zeros(6)
        # cell lengths
        self._params[0:3] = [np.linalg.norm(i) for i in self.lattice][:]
        # angles in rad
        self._params[3:6] = [LinAlg.calc_angle(i, j) for i, j in
                            reversed(list(itertools.combinations(self.lattice, 2)))]

    def reorient_lattice(self):
        self.__mkparam()
        a, b, c = self._params[:3]
        al, be, ga = self._params[3:]
        cos_be = np.cos(be)
        cos_ga = np.cos(ga)
        sin_ga = np.sin(ga)
        cos_al = np.cos(al)
        c_x = c*cos_be
        c_y = c*(cos_al - cos_ga*cos_be)/sin_ga
        c_z = np.sqrt(c**2 - c_x**2 - c_y**2)
        self.lattice = np.array([[a, 0., 0.],
                                [b*cos_ga, b*sin_ga, 0.],
                                [c_x, c_y, c_z]])
        del self._ilattice # re-compute the inverse

    @property
    def a(self):
        """Magnitude of cell a vector."""
        return self._params[0]

    @property
    def b(self):
        """Magnitude of cell b vector."""
        return self._params[1]

    @property
    def c(self):
        """Magnitude of cell c vector."""
        return self._params[2]

    @property
    def alpha(self):
        """Cell angle alpha."""
        return self._params[3]*RAD2DEG

    @property
    def beta(self):
        """Cell angle beta."""
        return self._params[4]*RAD2DEG

    @property
    def gamma(self):
        """Cell angle gamma."""
        return self._params[5]*RAD2DEG

class Symmetry(object):
    def __init__(self, options):
        self.options = options
        assert os.path.isdir(options.symmetry_dir)
        sys.path.append(options.symmetry_dir)
        self.spg = __import__('pyspglib._spglib')._spglib
        #import pyspglib._spglib as spg
        self._symprec = options.symmetry_precision
        self._lattice = None
        self._inv_latt = None
        self._scaled_coords = None
        self._element_symbols = None
        self.dataset = {}

    def add_structure(self, structure):
        self._lattice = structure.cell.lattice.copy()
        self._inv_latt = structure.cell.inverse.copy()
        self._scaled_coords = np.array([atom.in_cell_scaled(self._inv_latt) for
                                        atom in structure.atoms])
        self._angle_tol = -1.0
        self._element_symbols = [atom.element for atom in structure.atoms]
        self._numbers = np.array([ATOMIC_NUMBER.index(i) for i in 
                                    self._element_symbols])

    def refine_cell(self):
        """
        get refined data from symmetry finding
        """
        # Temporary storage of structure info
        _lattice = self._lattice.T.copy()
        _scaled_coords = self._scaled_coords.copy()
        _symprec = self._symprec
        _angle_tol = self._angle_tol
        _numbers = self._numbers.copy()
        
        keys = ('number',
                'international',
                'hall',
                'transformation_matrix',
                'origin_shift',
                'rotations',
                'translations',
                'wyckoffs',
                'equivalent_atoms')
        dataset = {}

        dataset['number'] = 0
        while dataset['number'] == 0:

            # refine cell
            num_atom = len(_scaled_coords)
            ref_lattice = _lattice.copy()
            ref_pos = np.zeros((num_atom * 4, 3), dtype=float)
            ref_pos[:num_atom] = _scaled_coords.copy()
            ref_numbers = np.zeros(num_atom * 4, dtype=int)
            ref_numbers[:num_atom] = _numbers.copy()
            num_atom_bravais = self.spg.refine_cell(ref_lattice,
                                       ref_pos,
                                       ref_numbers,
                                       num_atom,
                                       _symprec,
                                       _angle_tol)
            for key, data in zip(keys, self.spg.dataset(ref_lattice.copy(),
                                    ref_pos[:num_atom_bravais].copy(),
                                ref_numbers[:num_atom_bravais].copy(),
                                            _symprec,
                                            _angle_tol)):
                dataset[key] = data

            _symprec = _symprec * 0.5

        # an error occured with met9, org1, org9 whereby no
        # symmetry info was being printed for some reason.
        # thus a check is done after refining the structure.

        if dataset['number'] == 0:
            warning("WARNING - Bad Symmetry found!")
        else:

            self.dataset['number'] = dataset['number']
            self.dataset['international'] = dataset['international'].strip()
            self.dataset['hall'] = dataset['hall'].strip()
            self.dataset['transformation_matrix'] = np.array(dataset['transformation_matrix'])
            self.dataset['origin_shift'] = np.array(dataset['origin_shift'])
            self.dataset['rotations'] = np.array(dataset['rotations'])
            self.dataset['translations'] = np.array(dataset['translations'])
            letters = "0abcdefghijklmnopqrstuvwxyz"
            try:
                self.dataset['wyckoffs'] = [letters[x] for x in dataset['wyckoffs']]
            except IndexError:
                print dataset['wyckoffs']
            self.dataset['equivalent_atoms'] = np.array(dataset['equivalent_atoms'])
            self._lattice = ref_lattice.T.copy()
            self._scaled_coords = ref_pos[:num_atom_bravais].copy()
            self._numbers = ref_numbers[:num_atom_bravais].copy()
            self._element_symbols = [ATOMIC_NUMBER[i] for 
                i in ref_numbers[:num_atom_bravais]]

    def get_space_group_name(self):
        return self.dataset["international"] 

    def get_space_group_operations(self):
        return [self.convert_to_string((r, t)) 
                for r, t in zip(self.dataset['rotations'], 
                    self.dataset['translations'])]

    def get_space_group_number(self):
        return self.dataset["number"]

    def get_equiv_atoms(self):
        """Returs a list where each entry represents the index to the
        asymmetric atom. If P1 is assumed, then it just returns a list
        of the range of the atoms."""
        return self.dataset["equivalent_atoms"]

    def get_equivalent_hydrogens(self):
        at_equiv = self.get_equiv_atoms()
        h_equiv = {}
        h_id = list(set([i for id, i in enumerate(at_equiv) 
                    if self._element_symbols[id] == "H"]))
        for id, i in enumerate(self._element_symbols):
            if i == "H":
                h_equiv[id] = h_id.index(at_equiv[id])
        return h_equiv
    
    def convert_to_string(self, operation):
        """ takes a rotation matrix and translation vector and
        converts it to string of the format "x, y, z" """
        def tofrac(x, largest_denom = 32):
        
            negfrac = False
            if not x >= 0:
                negfrac = True
                x = abs(x)
        
            scaled = int(round(x * largest_denom))
            whole, leftover = divmod(scaled, largest_denom)
            if leftover:
                while leftover % 2 == 0:
                    leftover >>= 1
                    largest_denom >>= 1
            if negfrac:
                return -1*whole, leftover, largest_denom
        
            else:
                return whole, leftover, largest_denom
       
        def to_x(val):
            """ assumes integer value returned to x """
            if val == 0:
                return ""
            elif val == 1:
                return "x"
            elif val == -1:
                return "-x"
            else:
                return "%ix"%val
        
        def to_y(val):
            if val == 0:
                return ""
            elif val == 1:
                return "y"
            elif val == -1:
                return "-y"
            else:
                return "%iy"%val
        
        def to_z(val):
            if val == 0:
                return ""
            elif val == 1:
                return "z"
            elif val == -1:
                return "-z"
            else:
                return "%iz"%val
        # operation[0][0] is the first entry,
        # operation[0][1] is the second entry,
        # operation[0][2] is the third entry,
        # operation[1][1, 2, 3] are the translations
        fracs = [tofrac(i) for i in operation[1]]
        string = ""
        for idx, op in enumerate(operation[0]):
            x,y,z = op
            # note, this assumes 1, 0 entries in the rotation operation
            str_conv = (to_x(x), "+"*abs(x)*(y) + to_y(y),
                        "+"*max(abs(x),abs(y))*(z) + to_z(z))
            # determine if translation needs to be included
            for ind, val in enumerate((x,y,z)):
                frac = ""
                if val and fracs[ind][1] != 0:
                    frac = "+%i/%i"%(fracs[ind][1], fracs[ind][2])
                string += str_conv[ind] + frac

            # function to add comma delimiter if not the last entry
            f = lambda p: p < 2 and ", " or ""
            string += f(idx) 

        return string
    
    @property
    def cell_setting(self):
        return {
                0   :   "triclinic",
                1   :   "triclinic",       
                2   :   "triclinic",       
                3   :   "monoclinic",      
                4   :   "monoclinic",      
                5   :   "monoclinic",      
                6   :   "monoclinic",      
                7   :   "monoclinic",      
                8   :   "monoclinic",      
                9   :   "monoclinic",      
               10   :   "monoclinic",      
               11   :   "monoclinic",      
               12   :   "monoclinic",      
               13   :   "monoclinic",      
               14   :   "monoclinic",      
               15   :   "monoclinic",      
               16   :   "orthorhombic",      
               17   :   "orthorhombic",      
               18   :   "orthorhombic",      
               19   :   "orthorhombic",      
               20   :   "orthorhombic",      
               21   :   "orthorhombic",      
               22   :   "orthorhombic",      
               23   :   "orthorhombic",      
               24   :   "orthorhombic",      
               25   :   "orthorhombic",      
               26   :   "orthorhombic",      
               27   :   "orthorhombic",      
               28   :   "orthorhombic",      
               29   :   "orthorhombic",      
               30   :   "orthorhombic",      
               31   :   "orthorhombic",      
               32   :   "orthorhombic",      
               33   :   "orthorhombic",      
               34   :   "orthorhombic",      
               35   :   "orthorhombic",      
               36   :   "orthorhombic",      
               37   :   "orthorhombic",      
               38   :   "orthorhombic",      
               39   :   "orthorhombic",      
               40   :   "orthorhombic",      
               41   :   "orthorhombic",      
               42   :   "orthorhombic",      
               43   :   "orthorhombic",      
               44   :   "orthorhombic",      
               45   :   "orthorhombic",      
               46   :   "orthorhombic",      
               47   :   "orthorhombic",      
               48   :   "orthorhombic",      
               49   :   "orthorhombic",      
               50   :   "orthorhombic",      
               51   :   "orthorhombic",      
               52   :   "orthorhombic",      
               53   :   "orthorhombic",      
               54   :   "orthorhombic",      
               55   :   "orthorhombic",      
               56   :   "orthorhombic",      
               57   :   "orthorhombic",      
               58   :   "orthorhombic",      
               59   :   "orthorhombic",      
               60   :   "orthorhombic",      
               61   :   "orthorhombic",      
               62   :   "orthorhombic",      
               63   :   "orthorhombic",      
               64   :   "orthorhombic",      
               65   :   "orthorhombic",      
               66   :   "orthorhombic",      
               67   :   "orthorhombic",      
               68   :   "orthorhombic",      
               69   :   "orthorhombic",      
               70   :   "orthorhombic",      
               71   :   "orthorhombic",      
               72   :   "orthorhombic",      
               73   :   "orthorhombic",      
               74   :   "orthorhombic",      
               75   :   "tetragonal",        
               76   :   "tetragonal",        
               77   :   "tetragonal",        
               78   :   "tetragonal",        
               79   :   "tetragonal",        
               80   :   "tetragonal",        
               81   :   "tetragonal",        
               82   :   "tetragonal",        
               83   :   "tetragonal",        
               84   :   "tetragonal",        
               85   :   "tetragonal",        
               86   :   "tetragonal",        
               86   :   "tetragonal",        
               87   :   "tetragonal",        
               88   :   "tetragonal",        
               89   :   "tetragonal",        
               90   :   "tetragonal",        
               91   :   "tetragonal",        
               92   :   "tetragonal",        
               93   :   "tetragonal",        
               94   :   "tetragonal",        
               95   :   "tetragonal",        
               96   :   "tetragonal",        
               97   :   "tetragonal",        
               98   :   "tetragonal",        
               99   :   "tetragonal",        
              100   :   "tetragonal",        
              101   :   "tetragonal",        
              102   :   "tetragonal",        
              103   :   "tetragonal",        
              104   :   "tetragonal",        
              105   :   "tetragonal",        
              106   :   "tetragonal",        
              107   :   "tetragonal",        
              108   :   "tetragonal",        
              109   :   "tetragonal",        
              110   :   "tetragonal",        
              111   :   "tetragonal",        
              112   :   "tetragonal",        
              113   :   "tetragonal",        
              114   :   "tetragonal",        
              115   :   "tetragonal",        
              116   :   "tetragonal",        
              117   :   "tetragonal",        
              118   :   "tetragonal",        
              119   :   "tetragonal",        
              120   :   "tetragonal",        
              121   :   "tetragonal",        
              122   :   "tetragonal",        
              123   :   "tetragonal",        
              124   :   "tetragonal",        
              125   :   "tetragonal",        
              126   :   "tetragonal",        
              127   :   "tetragonal",        
              128   :   "tetragonal",        
              129   :   "tetragonal",        
              130   :   "tetragonal",        
              131   :   "tetragonal",        
              132   :   "tetragonal",        
              133   :   "tetragonal",        
              134   :   "tetragonal",        
              135   :   "tetragonal",        
              136   :   "tetragonal",        
              137   :   "tetragonal",        
              138   :   "tetragonal",        
              139   :   "tetragonal",        
              140   :   "tetragonal",        
              141   :   "tetragonal",        
              142   :   "tetragonal",        
              143   :   "trigonal",          
              144   :   "trigonal",          
              145   :   "trigonal",          
              146   :   "rhombohedral",   
              147   :   "trigonal",       
              148   :   "rhombohedral",   
              149   :   "trigonal",       
              150   :   "trigonal",       
              151   :   "trigonal",       
              152   :   "trigonal",       
              153   :   "trigonal",       
              154   :   "trigonal",       
              155   :   "rhombohedral",   
              156   :   "trigonal",       
              157   :   "trigonal",       
              158   :   "trigonal",       
              159   :   "trigonal",       
              160   :   "rhombohedral",   
              161   :   "rhombohedral",   
              162   :   "trigonal",       
              163   :   "trigonal",       
              164   :   "trigonal",       
              165   :   "trigonal",       
              166   :   "rhombohedral",   
              167   :   "rhombohedral",   
              168   :   "hexagonal",      
              169   :   "hexagonal",      
              170   :   "hexagonal",      
              171   :   "hexagonal",      
              172   :   "hexagonal",      
              173   :   "hexagonal",      
              174   :   "hexagonal",      
              175   :   "hexagonal",      
              176   :   "hexagonal",      
              177   :   "hexagonal",      
              178   :   "hexagonal",      
              179   :   "hexagonal",      
              180   :   "hexagonal",      
              181   :   "hexagonal",      
              182   :   "hexagonal",      
              183   :   "hexagonal",      
              184   :   "hexagonal",      
              185   :   "hexagonal",      
              186   :   "hexagonal",      
              187   :   "hexagonal",      
              188   :   "hexagonal",      
              189   :   "hexagonal",      
              190   :   "hexagonal",      
              191   :   "hexagonal",      
              192   :   "hexagonal",      
              193   :   "hexagonal",      
              194   :   "hexagonal",      
              195   :   "cubic",          
              196   :   "cubic",          
              197   :   "cubic",          
              198   :   "cubic",          
              199   :   "cubic",          
              200   :   "cubic",          
              201   :   "cubic",          
              202   :   "cubic",          
              203   :   "cubic",          
              204   :   "cubic",          
              205   :   "cubic",          
              206   :   "cubic",          
              207   :   "cubic",          
              208   :   "cubic",          
              209   :   "cubic",          
              210   :   "cubic",          
              211   :   "cubic",          
              212   :   "cubic",          
              213   :   "cubic",          
              214   :   "cubic",          
              215   :   "cubic",          
              216   :   "cubic",          
              217   :   "cubic",          
              218   :   "cubic",          
              219   :   "cubic",          
              220   :   "cubic",          
              221   :   "cubic",          
              222   :   "cubic",          
              223   :   "cubic",          
              224   :   "cubic",          
              225   :   "cubic",          
              226   :   "cubic",          
              227   :   "cubic",          
              228   :   "cubic",          
              229   :   "cubic",          
              230   :   "cubic",          
              } 
    
