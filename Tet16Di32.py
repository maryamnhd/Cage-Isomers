from utils import *
import os
from ase import Atoms, Atom
from ase.io import read, write
from ase.neighborlist import *


def transformation_library(FGs):
    TL = []

    tFGs = FGs.copy()
    for deg in range(4):
        if deg!=0:
            tFGs.rotate(90,'z') 
        TL += [project(FGs,tFGs)] 

    tFGs = reflect(FGs.copy(),[0,0,1])
    for deg in range(4):
        if deg!=0:
            tFGs.rotate(90,'z') 
        TL += [project(FGs,tFGs)]    

    tFGs = reflect(FGs.copy(),[1,0,0])
    for deg in range(4):
        if deg!=0:
            tFGs.rotate(90,'z') 
        TL += [project(FGs,tFGs)]  
  
    tFGs = reflect(FGs.copy(),[0,0,1])
    tFGs = reflect(tFGs.copy(),[1,0,0])
    for deg in range(4):
        if deg!=0:
            tFGs.rotate(90,'z') 
        TL += [project(FGs,tFGs)]      
  
    TL.pop(0)
    return TL


# ================================================================================= #
#                            Generate cage/base structure
# ================================================================================= #

import stk

#Zr as metal centre
Zr_atom = stk.BuildingBlock(
    smiles='[Zr]',
    functional_groups=(stk.SingleAtom(stk.Zr(0)) for i in range(4)),
    position_matrix=[[0., 0., 0.]],
)

# Benzene as linker 
bb1 = stk.BuildingBlock('Nc1ccc(N)cc1',[stk.PrimaryAminoFactory()]) #c1cnccn1

# Construct cage structure
cage = stk.ConstructedMolecule(
    topology_graph=stk.cage.M12L24(building_blocks=(Zr_atom, bb1),optimizer=stk.Collapser(step_size=1.5, distance_threshold=6)),   #, distance_threshold=5
)

# Write to file
writer = stk.PdbWriter()
writer.write(molecule=cage, path='base_structure/Tet16Di32.pdb')



# ================================================================================= #
#                       Functionalisation and Preparation
# ================================================================================= #

# Applying full functionalisation to the base structure
cage = read("base_structure/Tet16Di32.pdb")
FGs, FG_indices = full_functionalisation(cage.copy())

# Number of functional group slot
nFGs = len(FGs)

# Number of linkers
nLinkers = int(nFGs/4)


# ================================================================================= #
#                       Calculate number of unique isomers
# ================================================================================= #

# Linker indexing 
linkers = [[i+(j*4) for i in range(4)] for j in range(nLinkers)]

# Transformation library
TL = transformation_library(FGs.copy())

# Number of unique isomers identified by Cauchy's theorem (group theory)
nisomers = count_unique_isomers(TL, nLinkers, nFGs)
print("Number of unique isomers = ", nisomers)

# Enumeration and identification of unique isomers
# this approach is avoided because the number of iterations is too large
if False:
    # Unique flag identifier list 
    uniqueflag = [None]*(4**nLinkers)

    # Create directory for isomer structure files
    path = 'Tet16Di32_isomers'
    if not os.path.exists(path): os.mkdir(path)

    TL = list(zip(*TL))

    iIs = 0
    while iIs < (4**nLinkers):
        if uniqueflag[iIs] == None:
            uniqueflag[iIs] = True
            Is = convert10to4(iIs, nLinkers)
            #generate_isomer_structure_file(cage, FGs, Is, FG_indices, path)
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
       
       



