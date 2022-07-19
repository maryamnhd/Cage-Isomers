from utils import *
import os
from ase import Atoms, Atom
from ase.io import read, write
from ase.neighborlist import *


def transformation_library(FGs):
    TL = []

    tFGs = FGs.copy()
    for deg in range(3):
        if deg!=0:
            tFGs.rotate(120,'z') 
        TL += [project(FGs,tFGs)] 

    tFGs = reflect(FGs.copy(),[0,0,1])
    for deg in range(3):
        if deg!=0:
            tFGs.rotate(120,'z') 
        TL += [project(FGs,tFGs)]   

    tFGs = reflect(FGs.copy(),[1,0,0])
    for deg in range(3):
        if deg!=0:
            tFGs.rotate(120,'z') 
        TL += [project(FGs,tFGs)]  
  
    tFGs = reflect(FGs.copy(),[0,0,1])
    tFGs = reflect(tFGs.copy(),[1,0,0])
    for deg in range(3):
        if deg!=0:
            tFGs.rotate(120,'z') 
        TL += [project(FGs,tFGs)]      

    TL.pop(0)
    return list(zip(*TL))

# ================================================================================= #
#                            Generate cage/base structure
# ================================================================================= #

import stk
from stk.molecular.topology_graphs.cage import Cage
from stk.molecular.topology_graphs.cage.vertices import LinearVertex, NonLinearVertex
from stk.molecular.topology_graphs.topology_graph import Edge

class TwoPlusThree(Cage):
    _vertex_prototypes = (
        NonLinearVertex(0, [0, 0, 1]),
        NonLinearVertex(1, [0, 0, -1]),
        LinearVertex(2, [np.sin(np.pi*2/3), np.cos(np.pi*2/3), 0], False),
        LinearVertex(3, [np.sin(np.pi*4/3), np.cos(np.pi*4/3), 0], False),
        LinearVertex(4, [0, 1, 0], False),
    )
    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[2]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(3, _vertex_prototypes[1], _vertex_prototypes[2]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[3]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[4]),
    )
    _num_windows = 3
    _num_window_types = 1

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
    topology_graph=TwoPlusThree(building_blocks=(Zr_atom, bb1),optimizer=stk.Collapser(step_size=0.01, distance_threshold=5)),   
)
 
# Write to file
writer = stk.PdbWriter()
writer.write(molecule=cage, path='base_structure/Tri2Di3.pdb')

# ================================================================================= #
#                       Functionalisation and Preparation
# ================================================================================= #

# Applying full functionalisation to the base structure
cage = read("base_structure/Tri2Di3.pdb")
FGs, FG_indices = full_functionalisation(cage.copy())

# Number of functional group slot
nFGs = len(FGs)

# Number of linkers
nLinkers = int(nFGs/4)


# Unique flag identifier list 
uniqueflag = [None]*(4**nLinkers)

# Create directory for isomer structure files
path = 'Tri2Di3_isomers'
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
        generate_isomer_structure_file(cage.copy(), FGs, Is, FG_indices, path)
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
   


