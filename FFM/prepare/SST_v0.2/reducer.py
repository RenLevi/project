from interface import Reducer
import numpy as np
from sklearn.cluster import KMeans
import ase


class KMeansCluster(Reducer):
    def __init__(self,inputStructures: list[ase.Atoms]):
        self._inputStructures=inputStructures

    def setDescriptor(self,descriptors: np.ndarray[np.ndarray[np.float32]]) -> None:
        self._descriptors=descriptors

    def reduce(self,similarIdxs: list[int],kCount: int) -> list[int]:
        if len(similarIdxs)<=kCount: return list(similarIdxs)

        kMeans=KMeans(n_clusters=kCount,init='k-means++',max_iter=1000,n_init=20)
        clusterModel=kMeans.fit(self._descriptors[similarIdxs])
        clusterCenters=clusterModel.cluster_centers_
        labels=clusterModel.labels_  #聚类成k组时，标签为range(0,k)

        #按聚类标签分组
        similarIdxGroupByLabel=[[] for _ in range(kCount)]
        for label,idx in zip(labels,similarIdxs):
            similarIdxGroupByLabel[label].append(idx)

        #在每个分好的组内找出离质心最近的结构
        selectedStructuresIdx=[-1]*kCount
        for label,center in enumerate(clusterCenters):
            minEuclidDis=0x3f3f3f3f
            for structId in similarIdxGroupByLabel[label]:
                if (currEuclidDis:=np.linalg.norm(center-self._descriptors[structId],ord=2))<minEuclidDis:
                    minEuclidDis=currEuclidDis
                    selectedStructuresIdx[label]=structId

        return selectedStructuresIdx

class StructureCluster:

    def __init__(self, inputStructures: list[ase.Atoms]):
        self.inputStructures=inputStructures if isinstance(inputStructures,list) else ase.io.read(inputStructures,':')
        self.descriptors=self.getDescriptors(self.inputStructures)

    # modify this method if you have better descriptors
    def getDescriptors(self, inputStructures: list[ase.Atoms]):
        composition=inputStructures[0].symbols
        elements=set(composition)
        soapGenerator=SOAP(species=elements,r_cut=6,n_max=8,l_max=6,average='inner',sparse=False,periodic=True)
        descriptors=[soapGenerator.create(struct) for struct in inputStructures]
        return descriptors

    def findBestKCount(self, kCountTestRange: list[int]) -> None:
        shScoreRecords,chScoreRecords=list(),list()
        print('Starting find the best K from given values: ')
        for kCount in kCountTestRange:
            model=KMeans(n_clusters=kCount,init='k-means++',max_iter=1000,n_init=20)
            model.fit(self.descriptors)

            #clusterCenters=model.cluster_centers_
            labels=model.labels_
            shScore=silhouette_score(self.descriptors,labels,metric='euclidean')
            chScore=calinski_harabasz_score(self.descriptors,labels)

            shScoreRecords.append(shScore)
            chScoreRecords.append(chScore)
            print(f'Current K: {kCount:>3}, silhouette score: {shScore:.3f}, calinski harabasz score: {chScore:.3f}')

        fig,ax1=plt.subplots()
        ax1.plot(kCountTestRange,shScoreRecords,color='k',label='Silhouette')
        ax1.set_ylabel('Silhouette Score',color='k')
        plt.legend()
        ax2=ax1.twinx()
        ax2.plot(kCountTestRange,chScoreRecords,color='b',label='Calinski-Harabasz')
        ax2.set_ylabel('Calinski Harabasz Score',color='b')
        plt.legend()
        plt.show()

    def runCluster(self, kCount: int) -> list[ase.Atoms]:
        model=KMeans(n_clusters=kCount,init='k-means++',max_iter=1000,n_init=20).fit(self.descriptors)

        #自动计算分成簇的数量
        #model=MeanShift(bandwidth=0.5).fit(self.descriptors)
        #labels=model.labels_
        #shScore=silhouette_score(self.descriptors,labels,metric='euclidean')
        #chScore=calinski_harabasz_score(self.descriptors,labels)
        #print(f'current K: {len(np.unique(labels))}, silhouette score: {shScore:.3f}, calinski harabasz score: {chScore:.3f}')

        #自动计算分成簇的数量
        #model=AffinityPropagation(damping=0.5,max_iter=500).fit(self.descriptors)
        #labels=model.labels_
        #shScore=silhouette_score(self.descriptors,labels,metric='euclidean')
        #chScore=calinski_harabasz_score(self.descriptors,labels)
        #print(f'current K: {len(np.unique(labels))}, silhouette score: {shScore:.3f}, calinski harabasz score: {chScore:.3f}')

        clusterCenters=model.cluster_centers_
        labels=model.labels_

        #按聚类标签分组
        structIdGroupByCluster={label:[] for label in range(kCount)}
        for i,label in enumerate(labels):
            structIdGroupByCluster[label].append(i)

        #在每个分好的组内找出离质心最近的结构
        selectedStructures=[-1]*kCount
        selectedStructuresIdx=[-1]*kCount
        for label,center in enumerate(clusterCenters):
            minEuclidDis=0x3f3f3f3f
            for structId in structIdGroupByCluster[label]:
                if (currEuclidDis:=np.linalg.norm(center-self.descriptors[structId],ord=2))<minEuclidDis:
                    minEuclidDis=currEuclidDis
                    selectedStructures[label]=self.inputStructures[structId]
                    selectedStructuresIdx[label]=structId

        with open('idxFromSOAPCls_25k.txt','a+') as idxFile:
            idxFile.write(f"{self.inputStructures[0].get_chemical_formula(mode='hill',empirical=False)}\n")
            idxFile.write(','.join(str(idx) for idx in selectedStructuresIdx))
            idxFile.write('\n')

        if outputFilePath: ase.io.write(outputFilePath,selectedStructures)
        return selectedStructures

if __name__ == '__main__':
    import ase.io
    print(f'Reading initial dataset ...')
    ips=ase.io.read('Pt_CO_O.xyz',':1000')

    import random
    print('Generating descriptors of input structures ...')
    desp=np.array(list([random.random() for _ in range(1000)] for _ in range(len(ips))))

    from classifier import CosineSimilarity
    cst=CosineSimilarity(ips)
    cst.setDescriptor(desp)
    kmc=KMeansCluster(ips)
    kmc.setDescriptor(desp)

    print('Starting classification and reduction ...')
    for similarIdxs in cst.classify(0.25):
        print(f'similarIdxs = {similarIdxs}')
        reducedIdxs=kmc.reduce(similarIdxs,2)
        print(f'reducedIdxs = {reducedIdxs}')
        print('\n------------------\n')
