from decouple import config
LINUX = config('LINUX',cast=bool) or False
if not LINUX:
    from cpp_solver import NetGenerator, Net
else:
    import cppyy
    import os
    solver_path = os.path.dirname(os.path.realpath(__file__) )+ '/../solver'
    cppyy.include(solver_path + '/net.hpp')
    cppyy.include(solver_path + '/net.cpp')
    from cppyy.gbl import NetGenerator, Net

class NetFabric():
    def CreateNetDefault(netFile="..\\solver\\net.txt",borderFile="..\\solver\\border.txt",timeFile="..\\solver\\time.txt"):
        a = NetGenerator()
        a.isLogging = False;
        a.test(netFile,borderFile,timeFile)
        return a.GenerateFromFiles(netFile,borderFile,timeFile)

    def CreateNet():
      with open(netFile, 'r') as f:
            f.readline()
            r = f.readline().split(' ')
            r = [float(e) for e in r]
            f.readline()
            z = f.readline().split(' ')
            z = [float(e) for e in z]