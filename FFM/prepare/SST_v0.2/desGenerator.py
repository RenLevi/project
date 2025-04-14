from interface import DescriptorGenerator
from dscribe.descriptors import SOAP
import numpy as np
import ase

class StructSOAPGenerator(DescriptorGenerator):
    def __init__(self,elements: list[str],paraThreads: int):
        self._soapGenerator=SOAP(species=elements,r_cut=6,n_max=8,l_max=6,average='inner',
                                    sparse=False,periodic=True,dtype='float64')
        self._paraThreads=paraThreads

    def getDescriptor(self,inputStructures: list[ase.Atoms]) -> np.ndarray[np.ndarray[np.float64]]:
        descriptors=self._soapGenerator.create(inputStructures,n_jobs=self._paraThreads)

        #计算当前组分所有结构SOAP向量的模长
        #import matplotlib.pyplot as plt
        #for i,desc in enumerate(descriptors,start=1):
        #   print(f'{composition}:{i:>4}/{len(descriptors):>4} descriptor dimensions = {len(desc)}, length = {np.linalg.norm(desc,ord=2):.2f}')
        #fig=plt.figure(figsize=(7,5))
        #plt.plot(range(1,len(descriptors)+1),[np.linalg.norm(desc,ord=2) for desc in descriptors])
        #plt.show()

        return descriptors

class AtomSOAPGenerator(DescriptorGenerator):
    def __init__(self,elements: list[str]):
        self._soapGenerator=SOAP(species=elements,r_cut=6,n_max=8,l_max=6,average='off',
                                    sparse=False,periodic=True,dtype='float32')

    def getDescriptor(self,inputStructure: ase.Atoms,centerIdxs: list[int]) -> np.ndarray[np.ndarray[np.float32]]:
        descriptors=self._soapGenerator.create(inputStructure,centers=centerIdxs)

        #计算所有原子SOAP向量的模长
        #import matplotlib.pyplot as plt
        #for i,desc in enumerate(descriptors,start=1):
        #   print(f'{composition}:{i:>4}/{len(descriptors):>4} descriptor dimensions = {len(desc)}, length = {np.linalg.norm(desc,ord=2):.2f}')
        #fig=plt.figure(figsize=(7,5))
        #plt.plot(range(1,len(descriptors)+1),[np.linalg.norm(desc,ord=2) for desc in descriptors])
        #plt.show()

        return descriptors

if __name__=='__main__':
    import ase.io
    from dscribe.kernels import AverageKernel
    des=AtomSOAPGenerator(['Pt','C','O']).getDescriptor(ase.io.read('Pt_CO_O.xyz','0'),[1,2,3])
    print(des.shape)
    #re=AverageKernel(metric="rbf", gamma=1)
    #reK=re.create([des1,des2])
