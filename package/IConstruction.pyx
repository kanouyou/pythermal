
import petsc4py
import sys

from petsc4py import PETSc
from IMakeMaterial import IMakeMaterial
import numpy as np

class IConstruction:

    def __init__(self):
        self.fPara = {"rx":np.array([],dtype="float"), "ry":np.array([],dtype="float"), \
                      "Td":np.array([],dtype="float"),  "q":np.array([],dtype="float")}


    def SetTimeStep(self, dt):
        self.fTimeStep = dt


    def SetInitialTemp(self, T0):
        self.fTi = T0


    def SetBoundaryTemp(self, Tb):
        self.fTb = Tb


    def SetMaterial(self, makemat):
        mat = makemat.GetMaterial()
        self.fMsh  = makemat.GetMesh()
        self.fCell = makemat.GetCellSize()

        for i in range(len(mat["name"])):
            rho = mat["rho"][i]
            C   = mat[  "C"][i]
            kx  = mat[ "kx"][i]
            ky  = mat[ "ky"][i]
            q   = mat[  "q"][i]

            alpx = kx / rho / C
            alpy = ky / rho / C

            dx = self.fCell["dx"]
            dy = self.fCell["dy"]

            w = q / rho / C

            self.fPara["rx"] = np.append(self.fPara["rx"], alpx * self.fTimeStep / dx**2)
            self.fPara["ry"] = np.append(self.fPara["ry"], alpy * self.fTimeStep / dy**2)
            self.fPara["Td"] = np.append(self.fPara["Td"], 0.)
            self.fPara[ "q"] = np.append(self.fPara[ "q"], w)

        node = makemat.GetVertexBoundary()
        for i in node:
            self.fPara["q"][node] /= 4.

        node = makemat.GetBoundaryLine("all")
        for i in node:
            self.fPara["q"][node] /= 2.


    def MakeMatrix(self):
        nx = self.fMsh["x"]
        ny = self.fMsh["y"]
        n  = nx * ny

        # make matrix
        A = PETSc.Mat()
        A.create(PETSc.COMM_WORLD)
        A.setSizes((n, n))
        A.setUp()

        begin, end = A.getOwnershipRange()

        for node in range(begin,end):
            i = node // nx
            j = node - i*nx
            # fill the matrix except the boundary
            if i!=0 and j!=0 and i!=ny-1 and j!=nx-1:
                A.setValue(node, node-nx, -self.fPara["ry"][node])
                A.setValue(node,  node-1, -self.fPara["rx"][node])
                A.setValue(node,    node, 1. + 2.*self.fPara["rx"][node] + 2.*self.fPara["ry"][node])
                A.setValue(node,  node+1, -self.fPara["rx"][node])
                A.setValue(node, node+nx, -self.fPara["ry"][node])

        #A.assemble()
        return A


    def MakeVector(self):
        nx = self.fMsh["x"]
        ny = self.fMsh["y"]
        n  = nx * ny

        # make vector
        b = PETSc.Vec()
        b.create(PETSc.COMM_WORLD)
        b.setSizes(n)
        b.setUp()

        begin, end = b.getOwnershipRange()

        for node in range(begin,end):
            b.setValue(node, 0. + self.fPara["q"][node])

        return b


    def GetParameters(self):
        para = {"rx":np.array([]), "ry":np.array([])}

        for i in range(len(self.fPara["rx"])):
            para["rx"] = np.append(para["rx"], self.fPara["rx"][i])
            para["ry"] = np.append(para["ry"], self.fPara["ry"][i])

        return para
