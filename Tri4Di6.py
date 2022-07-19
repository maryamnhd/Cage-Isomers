from utils import *
import os
from ase import Atoms, Atom
from ase.io import read, write
from ase.neighborlist import *
from itertools import product


def transformation_library(FGs, Vs):
    tv = list(product([1,-1], repeat=4))
    TL = []
    for itv in tv: 
        tFGs = FGs.copy()
        for r in range(2):
            if r ==1: tFGs = reflect(tFGs,[1,0,0])
            for d in range(3):
                for iv in range(len(Vs)):  
                    tFGs.rotate(120,v=itv[iv]*Vs[iv],center=[0,0,0]) 
                    TL += [tuple(project(FGs,tFGs))]  
    return list(zip(*set(TL)))

# ================================================================================= #
#                            Generate cage/base structure
# ================================================================================= #

import stk

#Zr as metal centre
Zr_atom = stk.BuildingBlock(
    smiles='[Zr]',
    functional_groups=(stk.SingleAtom(stk.Zr(0)) for i in range(3)),
    position_matrix=[[0., 0., 0.]],
)

# Benzene as linker 
bb1 = stk.BuildingBlock('Nc1ccc(N)cc1',[stk.PrimaryAminoFactory()]) #c1cnccn1

# Construct cage structure
cage = stk.ConstructedMolecule(
    topology_graph=stk.cage.FourPlusSix(building_blocks=(Zr_atom, bb1),optimizer=stk.Collapser(step_size=0.01, distance_threshold=5)),   
)

# Write to file
writer = stk.PdbWriter()
writer.write(molecule=cage, path='base_structure/Tri4Di6.pdb')

# ================================================================================= #
#                       Functionalisation and Preparation
# ================================================================================= #

# Applying full functionalisation to the base structure
cage = read("base_structure/Tri4Di6.pdb")
FGs, FG_indices = full_functionalisation(cage.copy())

# Number of functional group slot
nFGs = len(FGs)

# Number of linkers
nLinkers = int(nFGs/4)

# Unique flag identifier list 
uniqueflag = [None]*(4**nLinkers)

# Create directory for isomer structure files
path = 'Tri4Di6_isomers'
if not os.path.exists(path): os.mkdir(path)

# ================================================================================= #
#                       Calculate number of unique isomers
# ================================================================================= #

# Linker indexing 
linkers = [[i+(j*4) for i in range(4)] for j in range(nLinkers)]

# Transformation library
TL = transformation_library(FGs.copy(), cage.positions[:4])

# Enumeration and identification of unique isomers
iIs = 0
while iIs < (4**nLinkers):
    if uniqueflag[iIs] == None:
        uniqueflag[iIs] = True
        Is = convert10to4(iIs, nLinkers)
        generate_isomer_structure_file(cage, FGs, Is, FG_indices, path)
        eq_isomers = [TL[i] for i in Is]
        eq_isomers = list(zip(*eq_isomers))
        eq_isomers = [sorted(ei) for ei in eq_isomers]
        for ei in eq_isomers:
            if higherthan(ei,Is):
                iei = convert4to10(ei, nLinkers)
                uniqueflag[iei] = False
    iIs += 1

nisomers = uniqueflag.count(True)
print("Number of unique isomers = ",nisomers)
print("Isomer structure files are generated to folder %s" %path)
   


