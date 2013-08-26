#!/usr/bin/env python
import sys
import os
import openbabel as ob
import pybel
import numpy
from element_properties import ATOMIC_NUMBER 
#=============================================================================
# Program Begins
#=============================================================================
order = {1:"S", 2:"D", 3:"T"}
intype = "mol"
toptype = "eta"
genstrtype = "default"
metal = "False"
index = "?"
path = "."
filelist = os.listdir(path)
mollist = []
for file in filelist:
    if file[-3:] == intype:
        mollist.append(file)
filelines = ""

# change this so that the headings are not [] but :
for ind, filename in enumerate(mollist):
    index = ind + 1
    headlines = ""; coordlines = ""; tablelines = ""; toplines = ""
    bond_const = ""
    pybel_mol = pybel.readfile(intype, filename).next()
    # Note, there seems to be a shift when the coordinates are reported.
    # this may be unacceptable when adding to an existing database

    uff = ob.OBForceField_FindForceField('uff')
    uff.Setup(pybel_mol.OBMol)
    uff.GetAtomTypes(pybel_mol.OBMol)
    file_name = filename.split(".mol")[0]
    index = file_name.lstrip('index')
    headlines += "[" + file_name + "]\n"
    headlines += "index = " + str(index) + "\n"
    headlines += "metal = " + metal + "\n"
    headlines += "topology = " + toptype + "\n"
    coordlines += "coordinates:\n"
    toplines += "connectivity:\n"
    tablelines += "table:\n"
    bond_const += "bond_constraints:\n"

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
            bond_const += '%5i%5i\n'%(i[0]+1, 2)
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

