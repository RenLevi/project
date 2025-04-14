from abc import ABC,abstractclassmethod
from typing import Union,Iterator
import numpy as np
import ase

class DescriptorGenerator(ABC):
    @abstractclassmethod
    def getDescriptor(self, **structuresAndConfig) -> Union[list,Iterator]:
        ...

    def __repr__(self) -> str:
        return f'<DescriptorGenerator> {self.__class__.__name__}'

class Classifier(ABC):
    @abstractclassmethod
    def __init__(self,inputStructures: list[ase.Atoms]):
        ...

    @abstractclassmethod
    def setDescriptor(self,descriptors: np.ndarray[np.ndarray[np.float32]]) -> None:
        ...

    @abstractclassmethod
    def classify(self, **classificationConfig) -> Union[list[list[int]],Iterator[list[int]]]:
        ...

    def __repr__(self) -> str:
        return f'<Classifier> {self.__class__.__name__}'

class Reducer(ABC):
    @abstractclassmethod
    def __init__(self, inputStructures: list[ase.Atoms]):
        ...

    @abstractclassmethod
    def setDescriptor(self,descriptors: np.ndarray[np.ndarray[np.float32]]) -> None:
        ...

    @abstractclassmethod
    def reduce(self, **reducerConfig) -> list[int]:
        ...

    def __repr__(self) -> str:
        return f'<Reducer>    {self.__class__.__name__}'
