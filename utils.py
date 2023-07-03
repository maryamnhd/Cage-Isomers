import numpy
from ase import Atoms, Atom
from ase.io import read, write
import os
from ase.neighborlist import *

def FG-FG_distance_count(cage, Is):

    pairds = list(itertools.combinations(Is, 2))

    inNs = []

    for ilinker in range(int(len(cage) / 4)):
        ds = []
        for islot in range(4):
            d = round(numpy.linalg.norm(cage[ilinker * 4 + islot].position), 1)
            ds += [d]
        ds = numpy.array(ds)
        inNs += [ilinker * 4 + islot for slots in numpy.where(ds == ds.min()) for islot in list(slots)]

    dskeys = []
    for paird in list(itertools.combinations(range(len(cage)), 2)): #pairds
        if int(paird[0]/4) == int(paird[1]/4):continue
        dskeys.append(round(cage.get_distance(paird[0], paird[1]), 2))
    dskeys = sorted(list(set(dskeys)))

    dds = dict((i, []) for i in dskeys)
    for paird in pairds:
        k = round(cage.get_distance(paird[0], paird[1]), 2)
        dds[k].append(paird)

    countds = dict((i, 0) for i in dskeys)

    ds = dict((i, 1) for i in pairds)

    for d in dskeys:
        for iis in dds[d]:
            countds[d] += ds[iis]

    # ==================================================================================
    count3ds = dict((i, [0, 0, 0]) for i in dskeys)

    for d in dskeys:
        for iis in dds[d]:
            in0 = int(iis[0] not in inNs)
            in1 = int(iis[1] not in inNs)
            count3ds[d][in0 + in1] += ds[iis]

    return dskeys,[countds[k] for k in dskeys], [count3ds[k][0] for k in dskeys], [count3ds[k][1] for k in dskeys], [count3ds[k][2] for k in dskeys]


def reflect(cage,v):
    norm = numpy.linalg.norm
    v /= norm(v)
    vp = cage.positions
    cage.positions[:] -= numpy.outer(2*numpy.dot(vp, v), v)
    return cage

def project(FGs, tFGs):
    rot = []
    nFGs = len(FGs)
    merged = tFGs+FGs
    for ia in range(0,nFGs):
        for ib in range(nFGs,nFGs*2):
            if merged.get_distance(ia,ib)<0.3:
                rot+=[ib-nFGs]
                break    
    return rot

def convert4to10(Is,nLinkers):
    cc = 0
    for i in range(nLinkers):
        rr = Is[i]%4
        if rr!=0:
            cc+=(rr*(4**(nLinkers-1-i)))
    return cc

def convert10to4(i,nLinkers):
    Is = []; rr = []
    while i > 0:
        rr.append(i % 4)
        i = i // 4
    rr = rr+[0]*(nLinkers-len(rr))
    for ir in range(len(rr)):
        Is += [ir*4+rr[nLinkers-1-ir]] 
    return Is

def higherthan(list1,list2):
    for i in range(len(list1)):
        if list1[i]>list2[i]:
            return True
        elif list1[i]<list2[i]:
            return False
    return False

def full_functionalisation(cage):
    Cs = neighbor_list('i', cage, {('C', 'H'): 1.3})
    Hs = neighbor_list('j', cage, {('C', 'H'): 1.3})
    FG_indices = []
    for i,j in zip(Cs,Hs):
        if cage[i].symbol == 'C':
            cage.set_distance(i,j,1.47,fix=0)
            cage[j].symbol = 'N'
            cage[j].position = [round(p,3) for p in cage[j].position]
            FG_indices += [j]
    return cage[FG_indices], FG_indices
    
               
def writexyz(path,atoms):
    sym = atoms.symbols
    pos = atoms.positions
    lines = ['%i\n'%len(atoms),'\n']
    for iline in range(len(atoms)):
        lines += ['{0:3}{1:8}{2:8}{3:8}\n'.format(sym[iline],pos[iline][0],pos[iline][1],pos[iline][2])]
    file = open(path,'w')
    file.writelines(lines)
    file.close()

def generate_isomer_structure_file(cage,FGs,Is,FG_indices,path ):
    isomer = cage.copy()
    isomer.symbols = ['X' if s == 'N' else s for s in isomer.symbols]
    fname = ""
    for i in Is:
        isomer[FG_indices[i]].symbol = "N"
        isomer[FG_indices[i]].position = FGs[i].position
        fname += str(i)+'-'
    writexyz("{}/{}.xyz".format(path,fname[:-1]), isomer)

def count_unique_isomers(TL, nlinkers, nFGs):

    # identity symmetry element
    cycsymm = [4**nlinkers]
    
    # iterate over other symmetry elements  
    for rot in TL:
        flag = []
        ncyc = {}

        # iterate all functional group slots to find cyclic group
        for i in range(nFGs):
            if i in flag:
                continue
            flag.append(i)
            j=i
            n=1

            # finding and listing the the cyclic group
            while (rot[j] != i ):
                n+=1
                j= rot[j]
                flag.append(j)

            # n is the size of the cyclic group, ncyc[n] is the frequency of the cyclic group
            if n not in ncyc.keys():
                ncyc[n]=1
            else:
                ncyc[n]+=1

            if n==2:
                if int(rot[i]/4) == int(i/4):
                    cycsymm.append(0)
                    ncyc = {}
                    break
        if ncyc:
            mult = 0
            for n in ncyc.keys():
                mult+= ncyc[n]*n
            mult/=sum(list(ncyc.values()))
            cycsymm.append(int(4**(nlinkers/mult)))
    nsymm = len(TL)+1
    return sum(cycsymm)/nsymm



