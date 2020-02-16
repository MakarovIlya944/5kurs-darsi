from cpp_solver import NetGenerator, Net

class NetFabric():
    def CreateNetDefault(netFile="..\\solver\\net.txt",borderFile="..\\solver\\border.txt",timeFile="..\\solver\\time.txt"):
        a = NetGenerator()
        return a.GenerateFromFiles(netFile,borderFile,timeFile)

    def CreateNet():
      with open(netFile, 'r') as f:
            f.readline()
            r = f.readline().split(' ')
            r = [float(e) for e in r]
            f.readline()
            z = f.readline().split(' ')
            z = [float(e) for e in z]