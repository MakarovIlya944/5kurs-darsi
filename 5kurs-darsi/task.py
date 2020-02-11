from matrix import MatrixFabric
from net import NetFabric
from solver import Solver

def RightPart(r,z):
    return r*r+z*z

class TaskFabric():
    def CreateTask(**kwargs):
        r = kwargs.get('r') or [0,1]
        z = kwargs.get('z') or [0,1]
        f = kwargs.get('f') or RightPart
        b = kwargs.get('b') or RightPart
        return Task()

class Task():

    Solver = 0
    Matrix = 0
    Net = 0
    Border = 0

    Result = 0

    def __init__(self):
        super().__init__()

    def Run(self):
        return 'result'
