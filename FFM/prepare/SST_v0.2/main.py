from desGenerator import StructSOAPGenerator,AtomSOAPGenerator
from classifier import CosineSimilarity
from reducer import KMeansCluster
from collections import defaultdict
from functools import reduce
import logging,os,sys,time,random
import ase,ase.data,ase.io
import numpy as np

logging.basicConfig(format='%(levelname)s %(asctime)s %(message)s',datefmt='[%Y/%m/%d %H:%M:%S]',level=logging.INFO)

class SST:
    def __init__(self,dataFile: str,resultFile: str,elements: list[str]):
        #show version info
        logging.info('*** WELCOME TO SST(Samples Select Twice) v0.2 ***')
        logging.info('*** Author: dlsyty,  Release date: 2024/03/26 ***')
        logging.warning('*** SST v0.2 is a developing version and may be unstable ***')

        #check IO file errors
        if not os.path.isfile(dataFile):
            raise IOError(f'DataFile \'{dataFile}\' seems invalid')
        if not resultFile.endswith('.xyz'):
            raise RuntimeError('Results should be dumped to xyz format file')
        self._dataFile=dataFile
        self._resultFile=resultFile

        #check elements list errors
        for element in elements:
            if not isinstance(element,str) or element not in ase.data.chemical_symbols:
                raise ValueError(f'Element \'{element}\' seems invalid')
        self._elements=elements

    def loadDataset(self) -> list[ase.Atoms]:
        #load initial dataset from `self._dataFile`
        logging.info(f'Loading initial dataset from {self._dataFile} ...')
        inputStructures=ase.io.read(self._dataFile,':')
        logging.info(f'{len(inputStructures)} structures loaded from {self._dataFile}')

        logging.info('Purging highly unstable and twisted structures from initial dataset ...')
        R_KSPACING=25  # see https://www.vasp.at/wiki/index.php/KPOINTS for more information
        stableStructures=[]
        for struct in inputStructures:  # purge highly unstable and twisted structures
            if struct.get_potential_energy()/len(struct)<-1.2 and \
               max(np.linalg.norm(struct.get_forces(),axis=1))<40 and \
               min(struct.cell.angles())>25 and max(struct.cell.angles())<155 and \
               max(struct.cell.lengths())<40 and \
               reduce(int.__mul__,[int(max(1,R_KSPACING*np.linalg.norm(recVector)+0.5)) for recVector in struct.cell.reciprocal()])<=350:
                stableStructures.append(struct)
        logging.info(f'{len(stableStructures)} structures selected from initial dataset, with {len(inputStructures)-len(stableStructures)} highly unstable and twisted structures purged')
        return stableStructures

    def runCompositionSST(self,clsLimit: float=0.05,redCount: int=1) -> None:
        #load initial dataset from `self._dataFile`
        inputStructures=self.loadDataset()

        #group initial dataset by composition
        logging.info('Grouping initial dataset by composition ...')
        compos2Idxs=defaultdict(list)
        for i,struct in enumerate(inputStructures):
            compos2Idxs[struct.get_chemical_formula(mode='hill',empirical=False)].append(i)
        logging.info(f'{len(compos2Idxs)} group(s) created from initial dataset')

        #generate SOAP descriptors
        logging.info(f'Elements participating in SOAP: {"  ".join(self._elements)}')
        logging.info('Generating SOAP descriptors for structures ...')
        structDescriptors=StructSOAPGenerator(self._elements,paraThreads=8).getDescriptor(inputStructures)
        logging.info('Generation of structure descriptors completed')

        #build selection pipeline
        logging.info('Building selection pipeline ...')
        clser=CosineSimilarity(inputStructures,paraThreads=8)
        clser.setDescriptor(structDescriptors)
        reder=KMeansCluster(inputStructures)
        reder.setDescriptor(structDescriptors)
        logging.info(f'Selection pipeline: {clser}, clsLimit = {clsLimit}')
        logging.info(f'                    {reder}, redCount = {redCount}')

        #enter main loop
        logging.info('Entering selection loop ...')
        selectedIdxs=[]
        for cpsId,(compos,unclassifiedIdxs) in enumerate(compos2Idxs.items()):
            #process current composition
            currComposSelectedIdxs=[]
            for similarIdxs in clser.classify(unclassifiedIdxs,diffLimit=clsLimit):
                currComposSelectedIdxs.extend(reder.reduce(similarIdxs,kCount=redCount))
            logging.info(f'{compos:>30}: {len(currComposSelectedIdxs):>5} of {len(unclassifiedIdxs):>7} selected [{(cpsId+1)*100//len(compos2Idxs):>3}% OK]')
            selectedIdxs+=currComposSelectedIdxs
        logging.info(f'{"Total":>30}: {len(selectedIdxs):>5} of {len(inputStructures):>7} selected')

        #dump selected structures and indexs
        self.dump(inputStructures,selectedIdxs)
        logging.info('Composition SST task finished successfully')

    def runAtomSST(self,clsLimit: float=0.05,redCount: int=1) -> None:
        #load initial dataset from `self._dataFile`
        inputStructures=self.loadDataset()

        #build selection pipeline
        logging.info('Building selection pipeline ...')
        clser=CosineSimilarity(inputStructures,paraThreads=8)
        reder=KMeansCluster(inputStructures)
        logging.info(f'Selection pipeline: {clser}, clsLimit = {clsLimit}')
        logging.info(f'                    {reder}, redCount = {redCount}')

        #build atom SOAP generator
        logging.info(f'Elements participating in SOAP: {"  ".join(self._elements)}')
        soapGene=AtomSOAPGenerator(self._elements)

        selectedStructIds=dict()  #map of elements -> selected structIds by this element
        for element in self._elements:
            #generate SOAP descriptors
            logging.info(f'Generating SOAP descriptors for element {element} ...')

            #employ a 2-loop method to decrease the use of memory
            #first loop of enumerating input structures: get the count of current element atoms,
            #so the size of var `atomDescriptors` can be known
            eleAtomIdOfStructs=[None]*len(inputStructures)  #map of structId -> current element atom ids
            idx2StructId=[]  #map of atom descriptor index -> structId
            for structId,struct in enumerate(inputStructures):
                atomIds=[i for i,symbol in enumerate(struct.symbols) if symbol==element]  #find atomIds of current element
                if not atomIds: continue  #skip, then eleAtomIdOfStructs[structId]==None
                eleAtomIdOfStructs[structId]=atomIds
                idx2StructId+=[structId]*len(atomIds)  #establish map of atom descriptor index -> structId

            #second loop of enumerating input structures: generate SOAP descriptors
            atomDescriptors=None  #wait for allocating memory
            atmDesSize=0
            for structId,struct in enumerate(inputStructures):
                atomIds=eleAtomIdOfStructs[structId]
                if not atomIds: continue
                des=soapGene.getDescriptor(struct,atomIds)
                if atomDescriptors is None:  #initialize and allocate memory for var `atomDescriptors`
                    atomDescriptors=np.empty((len(idx2StructId),len(des[0])),dtype=np.float32)
                atomDescriptors[atmDesSize:atmDesSize+len(des)]=des  #copy var `des` to var `atomDescriptors`
                atmDesSize+=len(des)
                del des  #del var `des` to deallocate memory

            logging.info(f'Generation of descriptors of {len(idx2StructId)} atom(s) of element {element} completed')

            #set descriptors for pipeline
            clser.setDescriptor(atomDescriptors)
            reder.setDescriptor(atomDescriptors)

            #enter main loop
            logging.info(f'Entering selection loop of element {element} ...')
            selectedIdxs=[]  #record the index of selected atom descriptors
            for similarIdxs in clser.classify(list(range(len(atomDescriptors))),diffLimit=clsLimit):
                selectedIdxs.extend(reder.reduce(similarIdxs,kCount=redCount))
            logging.info(f'{len(selectedIdxs)} atom(s) of element {element} selected')

            #free the memory used by atomDescriptors
            clser.setDescriptor(None)
            reder.setDescriptor(None)

            #map selected atom descriptors to structures
            selectedStructIds[element]={idx2StructId[idx] for idx in selectedIdxs}
            logging.info(f'{len(selectedStructIds[element])} structure(s) selected by selection of element {element}')

        #merge selected structure ids, remove those repeated and dump the final selected structures
        finalStructIds=set()
        for sids in selectedStructIds.values():
            finalStructIds.update(sids)
        logging.info(f'{len(finalStructIds)} of {len(inputStructures)} structures finally selected')
        self.dump(inputStructures,finalStructIds)
        logging.info('Atom SST task finished successfully')

    #dump selected structures, unselected structures and indexs
    def dump(self,inputStructures: list[ase.Atoms],selectedIdxs: list[int]) -> None:
        logging.info('Dumping selected structures for training, unselected structures for validation and indexs of selected structures ...')
        dumpTime=time.strftime('%Y%m%d_%H%M%S')
        xyzTrainFile=f'{self._resultFile[:-4]}_train_{dumpTime}.xyz'
        xyzValidFile=f'{self._resultFile[:-4]}_val_{dumpTime}.xyz'

        ase.io.write(xyzTrainFile,[inputStructures[idx] for idx in selectedIdxs])
        logging.info(f'{len(selectedIdxs)} selected structures for training dumped to file: {xyzTrainFile}')

        logging.info('Ratio of training data size to validation data size set to 10:1')
        validDataSize=len(selectedIdxs)//10  # Setting validation data size manually is also allowed

        assert len(selectedIdxs)+validDataSize<=len(inputStructures),'Too large validation data size'
        validIdxs=random.sample(list(set(range(len(inputStructures)))-set(selectedIdxs)),validDataSize)
        ase.io.write(xyzValidFile,[inputStructures[idx] for idx in validIdxs])
        logging.info(f'{validDataSize} unselected structures for validation dumped to file: {xyzValidFile}')

        idxTrainFile=f'{self._resultFile[:-4]}_idx_{dumpTime}.txt'
        with open(idxTrainFile,'w') as file:
            file.write(','.join(map(str,sorted(selectedIdxs))))
        logging.info(f'Indexs of selected structures dumped to file: {idxTrainFile}')

def runStructureSST(dataFile: str,resultFile: str,elements: list[str]):
    logging.info('*** WELCOME TO SST(Samples Select Twice) v0.1 ***')
    #check errors
    if not os.path.isfile(dataFile):
        raise RuntimeError('Invalid dataset file')
    if not resultFile.endswith('xyz'):
        raise RuntimeError('Results should be dumped to xyz format file')

    #load initial dataset from `dataFile`
    logging.info(f'Loading initial dataset from {dataFile} ...')
    inputStructures=ase.io.read(dataFile,':')
    logging.info(f'{len(inputStructures)} structures loaded from {dataFile}')

    #generate SOAP descriptors
    logging.info('Generating SOAP descriptors for structures ...')
    structDescriptors=StructSOAPGenerator(elements).getDescriptor(inputStructures)
    logging.info('SOAP descriptors Generation completed')

    #build selection pipeline
    logging.info('Building selection pipeline ...')
    clser=CosineSimilarity(inputStructures)
    clser.setDescriptor(structDescriptors)
    reder=KMeansCluster(inputStructures)
    reder.setDescriptor(structDescriptors)
    logging.info(f'Selection pipeline: {clser}')
    logging.info(f'                    {reder}')

    #enter main loop
    logging.info('Entering selection loop ...')
    selectedIdxs=[]
    for similarIdxs in clser.classify(diffLimit=0.05):
        selectedIdxs.extend(reder.reduce(similarIdxs,kCount=1))
    logging.info(f'Selected {len(selectedIdxs)} structures by SST')

    #dump selected structures and indexs
    logging.info('Dumping selected structures and indexs ...')
    dumpTime=time.strftime('%Y%m%d_%H%M%S')
    xyzResultFile=f'{resultFile[:-4]}_{dumpTime}.xyz'
    for idx in selectedIdxs:
        ase.io.write(xyzResultFile,inputStructures[idx],append=True)
    logging.info(f'Dumped structures to file: {xyzResultFile}')

    idxResultFile=f'{resultFile[:-4]}Idx_{dumpTime}.txt'
    with open(idxResultFile,'a') as file:
        file.write(','.join(str(idx) for idx in selectedIdxs))
    logging.info(f'Dumped indexs to file: {idxResultFile}')
    logging.info('SST task finished successfully')

if __name__=='__main__':
    sst=SST(dataFile=sys.argv[1],resultFile='try.xyz',elements=['Ru','C','H','O'])
    sst.runCompositionSST(clsLimit=0.00014,redCount=1)
