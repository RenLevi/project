from pymatgen.core.periodic_table import Element
class bond():
    def __init__(self, ele1,ele2,dis):
        self.ele1 = ele1
        self.ele2 = ele2
        self.length = dis
    def judge_bondorder(self):
        r_i = Element(self.ele1).atomic_radius
        r_j = Element(self.ele2).atomic_radius
        delta = 0.45
        if self.length <= r_i+r_j+delta:
            return 1
        else:
            return 0