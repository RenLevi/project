from interface import Classifier
import numpy as np
import numba as nb
from numba.core import types
import ase
from typing import Iterator

class CosineSimilarity(Classifier):
    def __init__(self,inputStructures: list[ase.Atoms],paraThreads: int):
        self._inputStructures=inputStructures
        nb.config.NUMBA_DEFAULT_NUM_THREADS=paraThreads

    def setDescriptor(self,descriptors: np.ndarray[np.ndarray[np.float32]]) -> None:
        self._descriptors=descriptors

    def classify(self,unsimilarIdxs: list[int],diffLimit: float) -> Iterator[np.ndarray[np.int32]]:
        unsimilarIdxs=np.array(unsimilarIdxs,dtype=np.int32)
        limitCosDis=1-diffLimit

        while len(unsimilarIdxs)>0:
            centerIdx=np.random.choice(unsimilarIdxs)
            similarIdxs,unsimilarIdxs=paraSplitIdxs(unsimilarIdxs,centerIdx,limitCosDis,self._descriptors)
            yield similarIdxs


#test by 3000 configurations with each 1000 dim random descriptor for 20 loops(diffLimit=0.005):
#with jit: 41.1s
#  no jit: 208.3s
#with jit and  2 cores parallelization: 49.4s
#with jit and  4 cores parallelization: 30.6s
#with jit and  8 cores parallelization: 19.0s
#with jit and 16 cores parallelization: 12.7s
@nb.njit(types.UniTuple(nb.int32[:],2)(nb.int32[:],nb.int32,nb.float64,nb.float64[:,::1]),nogil=True,parallel=True)
def paraSplitIdxs(calcIdxs: np.ndarray[np.int32],centerIdx: int,limitCosDis: float,
            descriptors: np.ndarray[np.ndarray[np.float64]]) -> tuple[np.ndarray[np.int32],np.ndarray[np.int32]]:
    centerVec=descriptors[centerIdx]
    centerVecLen=np.linalg.norm(centerVec,ord=2)

    THREADS=nb.config.NUMBA_DEFAULT_NUM_THREADS
    blockSize=(len(calcIdxs)+THREADS-1)//THREADS
    paraSimilarIdxs=np.empty((THREADS,blockSize),dtype=np.int32)
    paraUnsimilarIdxs=np.empty((THREADS,blockSize),dtype=np.int32)
    paraSips=[0]*THREADS
    paraUips=[0]*THREADS

    #for pid in nb.prange(THREADS):
    #   if pid*blockSize>=len(calcIdxs): continue
    #   for i in range(pid*blockSize,min((pid+1)*blockSize,len(calcIdxs))):
    #       idx=calcIdxs[i]
    #       currVec=descriptors[idx]
    #       currVecLen=np.linalg.norm(currVec,ord=2)
    #       cosDis=centerVec@currVec/(centerVecLen*currVecLen)
    #       if cosDis>=limitCosDis:
    #           sip=paraSips[pid]
    #           paraSimilarIdxs[pid][sip]=idx
    #           paraSips[pid]=sip+1
    #       else:
    #           uip=paraUips[pid]
    #           paraUnsimilarIdxs[pid][uip]=idx
    #           paraUips[pid]=uip+1

    for pid in nb.prange(THREADS):
        #9 calcIdxs with 4 threads: [0..3) [3..6) [6..9) [9..12), range [9..12) should be skipped
        if pid*blockSize>=len(calcIdxs): continue

        sip,uip=0,0
        for i in range(pid*blockSize,min((pid+1)*blockSize,len(calcIdxs))):
            idx=calcIdxs[i]
            currVec=descriptors[idx]
            currVecLen=np.linalg.norm(currVec,ord=2)
            cosDis=centerVec@currVec/(centerVecLen*currVecLen)
            if cosDis>=limitCosDis:
                paraSimilarIdxs[pid][sip]=idx
                sip+=1
            else:
                paraUnsimilarIdxs[pid][uip]=idx
                uip+=1
        paraSips[pid]=sip
        paraUips[pid]=uip

    similarIdxs=paraSimilarIdxs[0][:paraSips[0]]
    unsimilarIdxs=paraUnsimilarIdxs[0][:paraUips[0]]
    for pid in range(1,THREADS):
        similarIdxs=np.append(similarIdxs,paraSimilarIdxs[pid][:paraSips[pid]])
        unsimilarIdxs=np.append(unsimilarIdxs,paraUnsimilarIdxs[pid][:paraUips[pid]])

    #similarIdxsLen=np.sum(paraSips)
    #similarIdxs=np.empty(similarIdxsLen,dtype=np.int32)
    #unsimilarIdxs=np.empty(len(calcIdxs)-similarIdxsLen,dtype=np.int32)
    #msip,muip=0,0  #idx for final merged similar idx and unsimilar idx

    #for pid in range(THREADS):
    #   #9 calcIdxs with 4 threads: [0..3) [3..6) [6..9) [9..12), range [9..12) should be skipped
    #   if pid*blockSize>=len(calcIdxs): continue
    #   sis,uis=paraSips[pid],paraUips[pid]  #current similar idx size & current unsimilar idx size
    #   similarIdxs[msip:msip+sis]=paraSimilarIdxs[pid][0:sis]
    #   unsimilarIdxs[muip:muip+uis]=paraUnsimilarIdxs[pid][0:uis]
    #   msip+=sis
    #   muip+=uis

    return similarIdxs,unsimilarIdxs

#test by 3000 configurations with each 1000 dim random descriptor for 20 loops(diffLimit=0.005):
#with jit: 41.2s
#  no jit: 203.1s
@nb.njit(types.UniTuple(nb.int32[:],2)(nb.int32[:],nb.int32,nb.float32,nb.float32[:,::1]),nogil=True)
def splitIdxs(calcIdxs: np.ndarray[np.int32],centerIdx: int,limitCosDis: float,
            descriptors: np.ndarray[np.ndarray[np.float32]]) -> tuple[np.ndarray[np.int32],np.ndarray[np.int32]]:
    centerVec=descriptors[centerIdx]
    centerVecLen=np.linalg.norm(centerVec,ord=2)

    similarIdxs=np.empty(len(calcIdxs),dtype=np.int32)
    unsimilarIdxs=np.empty(len(calcIdxs),dtype=np.int32)
    sip,uip=0,0

    for idx in calcIdxs:
        currVec=descriptors[idx]
        currVecLen=np.linalg.norm(currVec,ord=2)
        cosDis=centerVec@currVec/(centerVecLen*currVecLen)
        if cosDis>=limitCosDis:
            similarIdxs[sip]=idx
            sip+=1
        else:
            unsimilarIdxs[uip]=idx
            uip+=1
    return similarIdxs[:sip],unsimilarIdxs[:uip]

#def splitIdxs(calcIdxs: np.ndarray[np.int32],centerIdx: int,limitCosDis: float,
#           descriptors: np.ndarray[np.ndarray[np.float32]]) -> list[np.ndarray[np.int32],np.ndarray[np.int32]]:
#   allVecs=descriptors[calcIdxs]
#   allVecLens=np.linalg.norm(allVecs,ord=2,axis=1)
#   centerVec=descriptors[centerIdx]
#   centerVecLen=np.linalg.norm(centerVec,ord=2)
#
#   cosDiss=allVecs@centerVec/allVecLens/centerVecLen  #opt for parallelization
#   similarIdxs=np.where(cosDiss>=limitCosDis)[0]
#   unsimilarIdxs=np.where(cosDiss<limitCosDis)[0]
#   return [similarIdxs,unsimilarIdxs]

if __name__=='__main__':
    import ase.io,time
    ips=ase.io.read('Pt_CO_O.xyz',':')
    startTime=time.time()
    cst=CosineSimilarity(ips)
    cst.setDescriptor(np.array(list([random.random() for _ in range(1000)] for _ in range(len(ips)))))
    for _ in range(20):
        for idxs in cst.classify(list(range(len(ips))),0.005):
            pass
            #print(f'Received similar idxs: {idxs}')
    endTime=time.time()
    print(f'Elapsed time: {endTime-startTime:.2f}s')
