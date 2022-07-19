from utils import *
import os
from ase import Atoms, Atom
from ase.io import read, write
from ase.neighborlist import *


def transformation_library(FGs):
    TL = []
    Vs = [ 'z', '-z', 'y', '-y', 'x', '-x' ]
    for V in Vs:
        tFGs = FGs.copy()
        tFGs.rotate(V, 'z',center=[0,0,0])
        for deg in range(4):
            if deg!=0:
                tFGs.rotate(90,'z') 
            TL += [project(FGs,tFGs)] 

    for V in Vs:
        tFGs = reflect(FGs.copy(),[1,0,0])
        tFGs.rotate(V, 'z',center=[0,0,0])
        
        for deg in range(4):
            if deg!=0:
                tFGs.rotate(90,'z') 
            TL += [project(FGs,tFGs)] 

    TL.pop(0)
    return list(zip(*TL))




# ================================================================================= #
#                       Functionalisation and Preparation
# ================================================================================= #

# Applying full functionalisation to the base structure
cage = read("base_structure/Tet6Di12.pdb")
#cage.rotate(45,'z')
FGs, FG_indices = full_functionalisation(cage.copy())

# Number of functional group slot
nFGs = len(FGs)

# Number of linkers
nLinkers = int(nFGs/4)


# Unique flag identifier list 
uniqueflag = [None]*(4**nLinkers)

# Create directory for isomer structure files
path = 'Tet6Di12_isomers'
if not os.path.exists(path): os.mkdir(path)



# ================================================================================= #
#                       Calculate number of unique isomers
# ================================================================================= #

# Linker indexing 
linkers = [[i+(j*4) for i in range(4)] for j in range(nLinkers)]

# Transformation library
TL = transformation_library(FGs.copy())

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
   


