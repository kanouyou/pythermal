
import numpy as np

class IMakeNode:

    def __init__(self, nx, ny):
        # initialization
        self.fMsh   = {"x":nx, "y":ny}
        self.fSpace = {"x0":0., "y0":0., "lx":1., "ly":1.}
        self.fNode = {"id":np.array([], dtype="int"), "x":np.array([], dtype="float"),  \
                       "y":np.array([], dtype="float")}


    def id(self, i, j):
        return j + self.fMsh["x"]*i


    def GetMesh(self):
        return self.fMsh


    def SetSpace(self, x0=0., y0=0., lx=1., ly=1.):
        self.fSpace = {"x0":x0, "y0":y0, "lx":lx, "ly":ly}

        dx = self.fSpace["lx"] / self.fMsh["x"]
        dy = self.fSpace["ly"] / self.fMsh["y"]

        # meshing and set node id
        for i in range(self.fMsh["y"]):
            y = self.fSpace["y0"] + i*dy
            for j in range(self.fMsh["x"]):
                x = self.fSpace["x0"] + j*dx
                self.fNode["id"] = np.append(self.fNode["id"], self.id(i,j))
                self.fNode[ "y"] = np.append(self.fNode[ "y"], y)
                self.fNode[ "x"] = np.append(self.fNode[ "x"], x)


    def GetCellSize(self):
        dist = {"dx": self.fSpace["lx"]/self.fMsh["x"], "dy": self.fSpace["ly"]/self.fMsh["y"]}
        return dist


    def MakeRectangle(self, x0, y0, lx, ly):
        node = np.array([], dtype="int")

        for i in range(len(self.fNode["id"])):
            if self.fNode["x"][i]>=x0 and self.fNode["x"][i]<x0+lx and \
               self.fNode["y"][i]>=y0 and self.fNode["y"][i]<y0+ly:
                node = np.append(node, self.fNode["id"][i])

        return node


    def MakeCircle(self, x0=0., y0=0., R=1.):
        node = np.array([], dtype="int")

        for i in range(len(self.fNode["id"])):
            if (self.fNode["x"][i]-x0)**2 + (self.fNode["y"][i]-y0)**2 <= R**2:
                node = np.append(node, self.fNode["id"][i])

        return node


    def GetVertexBoundary(self):
        node = np.array([], dtype="int")

        node = np.append(node, self.fNode["id"][self.id(0,0)])
        node = np.append(node, self.fNode["id"][self.id(0,self.fMsh["x"]-1)])
        node = np.append(node, self.fNode["id"][self.id(self.fMsh["y"]-1,0)])
        node = np.append(node, self.fNode["id"][self.id(self.fMsh["y"]-1,self.fMsh["x"]-1)])

        return node


    def GetBoundaryLine(self, opt, lmin=0., lmax=1000.):
        node = np.array([], dtype="int")

        if opt=="t" or opt=="Top":
            for j in range(1, self.fMsh["x"]-1):
                if self.fNode["x"][self.id(self.fMsh["y"]-1,j)]>=lmin and \
                   self.fNode["x"][self.id(self.fMsh["y"]-1,j)]<=lmax:
                    node = np.append(node, self.fNode["id"][self.id(self.fMsh["y"]-1,j)])

        if opt=="b" or opt=="Bottom":
            for j in range(1, self.fMsh["x"]-1):
                if self.fNode["x"][self.id(0,j)]>=lmin and self.fNode["x"][self.id(0,j)]<=lmax:
                    node = np.append(node, self.fNode["id"][self.id(0,j)])

        if opt=="l" or opt=="Left":
            for i in range(1, self.fMsh["y"]-1):
                if self.fNode["y"][self.id(i,0)]>=lmin and self.fNode["y"][self.id(i,0)]<=lmax:
                    node = np.append(node, self.fNode["id"][self.id(i,0)])

        if opt=="r" or opt=="Right":
            for i in range(1, self.fMsh["y"]-1):
                if self.fNode["y"][self.id(i,self.fMsh["x"]-1)]>=lmin and \
                   self.fNode["y"][self.id(i,self.fMsh["x"]-1)]<=lmax:
                    node = np.append(node, self.fNode["id"][self.id(i,self.fMsh["x"]-1)])

        if opt=="all" or opt=="All":
            for i in range(1, self.fMsh["y"]-1):
                node = np.append(node, self.fNode["id"][self.id(i,0)])
                node = np.append(node, self.fNode["id"][self.id(i,self.fMsh["x"]-1)])
            for j in range(1, self.fMsh["x"]-1):
                node = np.append(node, self.fNode["id"][self.id(0,j)])
                node = np.append(node, self.fNode["id"][self.id(self.fMsh["y"]-1,j)])

        return node


    def GetNodeInfo(self):
        return self.fNode


    def GetBoundaryX(self, x0):
        node = np.array([], dtype="int")

        return node
