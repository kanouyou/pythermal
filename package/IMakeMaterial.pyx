
from IMakeNode import IMakeNode
import matplotlib.pyplot as plt
import numpy as np

class IMakeMaterial(IMakeNode):

    def __init__(self, nx, ny):
        IMakeNode.__init__(self,nx,ny)

        self.fMat = {"rho":np.array([], dtype="float"),    "C":np.array([], dtype="float"), \
                      "kx":np.array([], dtype="float"),   "ky":np.array([], dtype="float"), \
                       "q":np.array([], dtype="float"), "name":np.array([], dtype="str")}

        # initialize
        for i in range(ny):
            for j in range(nx):
                self.fMat[ "rho"] = np.append(self.fMat[ "rho"], 1.)
                self.fMat[   "C"] = np.append(self.fMat[   "C"], 1.)
                self.fMat[  "kx"] = np.append(self.fMat[  "kx"], 0.)
                self.fMat[  "ky"] = np.append(self.fMat[  "ky"], 0.)
                self.fMat[   "q"] = np.append(self.fMat[   "q"], 0.)
                self.fMat["name"] = np.append(self.fMat["name"], "Air")


    def SetMaterial(self, node, rho=1., C=1., kx=0., ky=0., q=0., name="Air"):
        self.fMat[ "rho"][node] = rho
        self.fMat[   "C"][node] = C
        self.fMat[  "kx"][node] = kx
        self.fMat[  "ky"][node] = ky
        self.fMat[   "q"][node] = q
        self.fMat["name"][node] = str(name)


    def GetMaterial(self):
        return self.fMat


    def CheckMaterial(self):
        mat = [self.fMat["name"][0]]

        # count
        for i in range(len(self.fMat["name"])):
            cnt = 0
            for j in range(len(mat)):
                if self.fMat["name"][i]==mat[j]:
                    cnt += 1
            if cnt==0:
                mat.append(self.fMat["name"][i])

        return mat


    def Plot(self, with_num=False, save=False):
        mat = self.CheckMaterial()
        posz = [ [] for i in range(len(mat)) ]
        posr = [ [] for i in range(len(mat)) ]

        for i in range(len(mat)):
            for j in range(len(self.fMat["name"])):
                if self.fMat["name"][j]==mat[i]:
                    posz[i].append(self.fNode["x"][j])
                    posr[i].append(self.fNode["y"][j])

            plt.plot(posz[i], posr[i], marker="o", ls="none", markeredgewidth=0, label=mat[i])
        plt.legend(loc=1, numpoints=1)
        plt.title("Node")
        plt.xlabel(" Z [m]")
        plt.ylabel(" R [m]")

        # plot the node id
        if with_num==True:
            for i in range(len(self.fNode["id"])):
                plt.text(self.fNode["x"][i], self.fNode["y"][i], "%i" %self.fNode["id"][i])

        if save==True:
            plt.savefig("node.pdf")
        plt.show()


