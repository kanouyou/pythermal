
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from petsc4py import PETSc

class IPostprocess:

    def __init__(self, mat):
        self.fMsh   = mat.GetMesh()
        self.fSpace = mat.GetCellSize()
        self.fNode  = mat.GetNodeInfo()
        self.fRes   = np.array([])
        self.fMat   = mat.GetMaterial()
        self.fSize  = mat.GetCellSize()
        self.dt     = 1.


    def SetTimeStep(self, dt):
        self.dt = dt


    def SetResult(self, T):
        begin, end = T.getOwnershipRange()

        for i in range(begin, end):
            self.fRes = np.append(self.fRes, T.getValues(i))


    def getmax(self):
        return max(self.fRes)


    def getmin(self):
        return min(self.fRes)


    def PlotMatTogether(self, save=False):
        fig = plt.figure(figsize=(12,8))
        plt.subplot(3,2,1)
        self.PlotMat("rx")

        plt.subplot(3,2,2)
        self.PlotMat("ry")

        plt.subplot(3,2,3)
        self.PlotMat("kx")

        plt.subplot(3,2,4)
        self.PlotMat("ky")

        plt.subplot(3,2,5)
        self.PlotMat("q")

        plt.subplot(3,2,6)
        self.plot2d()

        plt.subplots_adjust(hspace=0.7)
        plt.subplots_adjust(left=0.1)

        if save==True:
            plt.savefig("result.pdf")


    def PlotMat(self, opt):
        nx = self.fMsh["x"]
        ny = self.fMsh["y"]

        x = np.array([])
        y = np.array([])
        z = np.array([])

        vmin = 0.
        vmax = 0.

        for i in range(len(self.fMat["rho"])):
            x = np.append(x, self.fNode["x"][i])
            y = np.append(y, self.fNode["y"][i])
            if opt=="rho" or opt=="Density":
                z = np.append(z, self.fMat ["rho"][i])
            elif opt=="C" or opt=="Capacity":
                z = np.append(z, self.fMat ["C"][i])
            elif opt=="kx":
                z = np.append(z, self.fMat ["kx"][i])
            elif opt=="ky":
                z = np.append(z, self.fMat ["ky"][i])
            elif opt=="q":
                z = np.append(z, self.fMat ["q"][i])
            elif opt=="rx":
                alp = self.fMat["kx"][i]/self.fMat["rho"][i]/self.fMat["C"][i]
                rx  = alp * self.dt / self.fSize["dx"]**2
                z = np.append(z, rx)
            elif opt=="ry":
                alp = self.fMat["ky"][i]/self.fMat["rho"][i]/self.fMat["C"][i]
                ry  = alp * self.dt / self.fSize["dy"]**2
                z = np.append(z, ry)

        vmin = min(z)
        vmax = max(z)

        x = x.reshape((ny,nx))
        y = y.reshape((ny,nx))
        z = z.reshape((ny,nx))

        #if opt=="kx" or opt=="ky":
        plt.contourf(x,y,z,80, cmap="Blues", vmin=vmin, vmax=vmax)
        #else:
            #plt.contourf(x,y,z,60,cmap="RdBu")
        #plt.contourf(x,y,z,90,cmap=plt.cm.jet,vmin=z.min(),vmax=z.max())
        plt.title(opt)
        plt.xlabel("Z [m]")
        plt.ylabel("R [m]")
        cb = plt.colorbar(extend="both")
        if opt=="rho" or opt=="Density":
            cb.set_label("Density [kg/m$^3$]")
        elif opt=="C" or opt=="Capacity":
            cb.set_label("Heat Capacity [J/kg/K]")
        elif opt=="kx":
            cb.set_label("Thermal Conductivity [W/m/K]")
        elif opt=="ky":
            cb.set_label("Thermal Conductivity [W/m/K]")
        elif opt=="q":
            cb.set_label("Heat Generation [W/m$^3$]")
        elif opt=="rx" or opt=="ry":
            cb.set_label("r")



    def plot2d(self, save=False):
        nx = self.fMsh["x"]
        ny = self.fMsh["y"]

        x = np.array([])
        y = np.array([])
        z = np.array([])

        for i in range(len(self.fRes)):
            x = np.append(x, self.fNode["x"][i])
            y = np.append(y, self.fNode["y"][i])
            z = np.append(z, self.fRes[i])

        x = x.reshape((ny,nx))
        y = y.reshape((ny,nx))
        z = z.reshape((ny,nx))

        #plt.contourf(x,y,z,50,cmap=plt.cm.rainbow)
        plt.contourf(x,y,z,90,cmap=plt.cm.jet)
        #plt.contourf(x,y,z,90,cmap=plt.cm.jet,vmin=z.min(),vmax=z.max())
        plt.title("Temperature")
        plt.xlabel("Z [m]")
        plt.ylabel("R [m]")
        cb = plt.colorbar()
        cb.set_label("Temperature [K]")
        if save==True:
            plt.savefig("temp2d.pdf")

