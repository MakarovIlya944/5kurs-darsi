
class NetFabric():
    def CreateNet(self):
        return Net()

class Net():

    rMin = 0
    rMax = 1
    R = []

    zMin = 0
    zMax = 1
    Z = []

    NVTR = []
    NVCAT = []

    def __init__(self, **kwargs):
        self.RMin = kwargs['r'][0]
        self.RMax = kwargs['r'][1]

        self.ZMin = kwargs['z'][0]
        self.ZMax = kwargs['z'][1]

        n = kwargs['n']
        h = (self.ZMax - self.ZMin) / n