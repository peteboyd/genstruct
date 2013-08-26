#!/usr/bin/env python
import sys
import os
import openbabel as ob
import pybel
import numpy
from element_properties import ATOMIC_NUMBER 

intype = "mol"
toptype = "eta"
genstrtype = "default"
metal = "False"
index = "?"
path = "."
filelist = os.listdir(path)
mollist = []

def clean(name):
    if name[-4:] == ".mol":
        return name[:-4]
    return name

class InputSBU(object):
    """Contains the necessary information to produce an input for
    Genstruct. This input file is a necessary step in case bonding
    flags or symmetry are incorrect."""
    def __init__(self, filename):
        self.data = {'index':'', 'metal':'', 'topology':'', 'parent':''
                'atomic_info':'', 'bond_table':'', 'connectivity':'',
                'connect_flag':'', 'connect_sym':''}
        self.name = clean(filename) 
        self.mol = pybel.readfile('mol', filename).next()
        self._reset_formal_charges()

    def get_metal(self):
        if "m" in self.name[-2:]:
            self.update(metal="True")
        else:
            self.update(metal="False")

    def special(self):
        """If the .mol file ends with an 's', this will interpret
        it as a child SBU, the parent will be the mol name before the 's'"""
        if "s" in self.name[-1:]:
            self.upate(parent=self.name[:-1])

    def add_data(self, **kwargs):
        self.data.update(kwargs)

    def update(self, **kwargs):
        for key, val in kwargs.items()
            self.data[key] += val

    def _reset_formal_charges(self):
        """Set all formal charges to zero, this is how special
        information will be passed to oBMol objects."""
        for atom in self.mol:
            atom.OBAtom.SetFormalCharge(0)

    def _remove_atoms(self, *args):
        for obatom in args:
            self.mol.OBMol.DeleteAtom(obatom)

    def get_connect_info(self):
        """Grab all the atoms which are flagged by this program to be
        connectivity points. Namely, Xe, Y, and Rn. Ac series
        elements are replacement Xe atoms for special bonding purposes.
        """
        special = []
        connect_index = 0
        for ind, atom in enumerate(self.mol):
            N = atom.atomicnum
            if N == 54 or (N >= 89 and N <= 102):
                connect_index += 1
                X = "%12.4f %8.4f %8.4f"%(atom.coords)
                if (N >= 89 and N <= 102):
                    special.append((connect_index-1, N%89+1))
                for neighbour in ob.OBAtomAtomIter(atom.OBAtom):
                    x = neighbour.GetX() - atom.coords[0]
                    y = neighbour.GetY() - atom.coords[1]
                    z = neighbour.GetZ() - atom.coords[2]
                    if neighbour.GetAtomicNum() == 39:
                        net_atom = neighbour
                        net_vector = "%12.4f %8.4f %8.4f"%(x, y, z)
                    elif neighbour.GetAtomicNum() == 86:
                        bond_atom = neighbour
                        bond_vector = "%12.4f %8.4f %8.4f"%(x, y, z)
                    else:
                        neighbour.SetFormalCharge(connect_index)

                con_line = "".join([X, bond_vector, net_vector, "\n"])
                self.update(connectivity=con_line)
                self._remove_atoms(atom.OBAtom, net_atom, bond_atom)

        # include special considerations
        for (i, spec) in special:
            if spec == 89:
                bond_partner = 1
            elif spec == 90:
                bond_partner = 2
            else:
                bond_partner = 0
            const_line = '%5i%5i%5i\n'%(i, spec, bond_partner)
            self.update(connect_flag = const_line)

    def get_atom_info(self):
        for atom in self.mol:
            N = atom.OBAtom.GetAtomicNum()
            element = ATOMIC_NUMBER[N]
            coordlines = "%4s  %-6s %8.4f %8.4f %8.4f\n"%(
                    element, self._get_ff_type(atom), atom.coords[0],
                    atom.coords[1], atom.coords[2])
            self.update(atomic_info=coordlines)
            if atom.formalcharge != 0:
                conn_atom = str(atom.formalcharge) + "C"
                order = "S" # currently set to a single bond.
                tableline = "%4i%4s%4s\n"%(atom.idx-1, conn_atom, order)
                self.update(bond_table=tableline)

    def get_bond_info(self):
        for bond in ob.OBMolBondIter(self.mol.OBMol):
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            type = self.return_bondtype(bond)
            line = "%4i%4i%4s\n"%(start_idx-1, end_idx-1, type)

    def return_bondtype(self, bond):
        start_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if bond.IsSingle():
            return "S"
        elif bond.IsDouble():
            return "D"
        elif bond.IsTriple():
            return "T"
        elif bond.IsAromatic():
            return "A"
        elif start_atom.GetType()[-1] == "R" and end_atom.GetType()[-1] == "R"\
                and start_atom.ExplicitHydrogenCount() == 1 and\
                end_atom.ExplicitHydrogenCount() == 1:
            return "A"
        elif bond.IsAmide():
            return "Am"

    def _set_uff(self):
        """Adds UFF atomtyping to the openbabel molecule description"""
        uff = ob.OBForceField_FindForceField('uff')
        uff.Setup(self.mol.OBMol)
        uff.GetAtomTypes(self.mol.OBMol)
   def _get_ff_type(self, pyatom):
       return pyatom.OBAtom.GetData("FFAtomType").GetValue()

for ind, filename in enumerate(mollist):
    index = ind + 1
    pybel_mol = pybel.readfile(intype, filename).next()
    # Note, there seems to be a shift when the coordinates are reported.
    # this may be unacceptable when adding to an existing database

    uff = ob.OBForceField_FindForceField('uff')
    uff.Setup(pybel_mol.OBMol)
    uff.GetAtomTypes(pybel_mol.OBMol)
    file_name = filename.split(".mol")[0]
    index = file_name.lstrip('index')
    headlines = "[" + file_name + "]\n"
    headlines += "index = " + str(index) + "\n"
    headlines += "metal = " + metal + "\n"
    headlines += "topology = " + toptype + "\n"
    coordlines = "coordinates:\n"
    toplines = "connectivity:\n"
    tablelines = "table:\n"
    connect_lines = "connect_flag:\n"
    connect_sym_lines = "connect_sym:\n"

    # Initial scan for the connectivity atoms
    conn_count = 0
    connectpoints = []
    special = []
    connect_dic = {}
    for ob_atom in pybel_mol:
        # use formal charge to keep track of the connectivity bonding
        ob_atom.OBAtom.SetFormalCharge(0)

        atnum = ob_atom.atomicnum
        # Need to include some connectivity info (unique bonds etc..)
        if atnum == 54 or (atnum >= 89 and atnum <= 102):
            conn_count += 1
            # need to change this to input reading for molecules with
            # different bonding types.
            btype = "%5i"%(conn_count)
            if atnum >= 89 and atnum <= 102:
                special.append([conn_count-1, atnum])
            X = "%12.3f%8.3f%8.3f"%(ob_atom.coords)
            for neighbour in ob.OBAtomAtomIter(ob_atom.OBAtom):
                xcoord = neighbour.GetX() - ob_atom.coords[0]
                ycoord = neighbour.GetY() - ob_atom.coords[1]
                zcoord = neighbour.GetZ() - ob_atom.coords[2]
                if neighbour.GetAtomicNum() == 39:
                    Y = "%12.3f%8.3f%8.3f"%(xcoord, ycoord, zcoord)
                    Yatom = neighbour
                    
                elif neighbour.GetAtomicNum() == 86:
                    Z = "%12.3f%8.3f%8.3f"%(xcoord, ycoord, zcoord)
                    Zatom = neighbour

                else:
                    # flag neighbour as a bonding atom to this
                    # connectivity point with btype
                    neighbour.SetFormalCharge(conn_count)

            connectpoints.append([btype, X, Z, Y, "%5i"%(1), "%8s"%("None"),"\n"])
            pybel_mol.OBMol.DeleteAtom(ob_atom.OBAtom)
            pybel_mol.OBMol.DeleteAtom(Zatom)
            pybel_mol.OBMol.DeleteAtom(Yatom)
       
    for ind, i in enumerate(special):
        connectpoints[i[0]][5] = "%8i"%((i[1]%89)+1)
        # intermetal constraints
        if i[1] == 89:
            connect_lines += '%5i%5i\n'%(i[0]+1, 2)
        elif i[1] == 90:
            bond_const += '%5i%5i\n'%(i[0]+1, 1)

    buttsoup = ["".join(i) for i in connectpoints]
    toplines += "".join(buttsoup)
        
    # Now scan for the coordinates
    for ob_atom in pybel_mol:
        atomic_number = ob_atom.OBAtom.GetAtomicNum()
        # Note, change these atomic numbers to something else.  Zn (30) is commonly
        # used as an SBU
        if atomic_number == 39 or atomic_number == 54 or atomic_number == 86:
            pass
        else:
        
            element = ATOMIC_NUMBER[ob_atom.atomicnum]
        
            #filelines += "%4i"%(ob_atom.OBAtom.GetIdx())
            coordlines += "%4s"%(element)
            # two spaces
            coordlines += "%2s"%("")
            coordlines += "%-6s"%(ob_atom.OBAtom.GetData("FFAtomType").GetValue())
            coordlines += "%8.3f%8.3f%8.3f\n"%(ob_atom.coords)
            if ob_atom.formalcharge != 0:
                soup = str(ob_atom.formalcharge) + 'C'
                bond_order = 'S'
                tablelines += "%4i%4s%4s\n"%(ob_atom.idx-1, soup, bond_order)

    for bond in ob.OBMolBondIter(pybel_mol.OBMol):
        start_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        start_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        bond_num = bond.GetBondOrder()
        # Test for carboxylate

        if bond.IsSingle():
            bond_order = "S"
        if bond.IsDouble():
            bond_order = "D"
        if bond.IsTriple():
            bond_order = "T"
        if bond.IsAromatic():
            bond_order = "A"
        if start_atom.GetType()[-1] == "R" and end_atom.GetType()[-1] == "R" \
           and start_atom.ExplicitHydrogenCount() == 1 and \
           end_atom.ExplicitHydrogenCount() == 1:
            bond_order = "A";
        if bond.IsAmide():
            bond_order = "Am";
            

        tablelines += "%4i%4i%4s\n"%(start_idx-1, end_idx-1, bond_order)

    filelines += headlines + coordlines + tablelines + toplines + bond_const

outfile = open("fordatabase.out", "w")
outfile.writelines(filelines)
outfile.close()

