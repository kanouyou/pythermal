
import numpy as np


class IBoundary:

    def __init__(self, nx, ny):
        self.fMsh = {"x":nx, "y":ny}
        self.fNode = 0
        self.fBound = {"x":"Neumann", "y":"Neumann"}
        self.fBc = {"x":0., "y":0.}


    def id(self,i,j):
        return j + self.fMsh["x"] * i


    def SetBoundary(self, node, xb="Neumann", xbc=0., yb="Neumann", ybc=0.):
        self.fNode = node
        self.fBound = {"x":xb, "y":yb}
        self.fBc = {"x":xbc, "y":ybc}


    def SetLineBoundary(self, node, b="Neumann", bc=0.):
        self.fNode = node
        self.fBound["x"] = b
        self.fBc["x"] = bc


    def SetParameters(self, rx, ry):
        self.rx = rx
        self.ry = ry


    def MakeVertexMat(self, A, b):
        node = self.fNode
        nx = self.fMsh["x"]
        ny = self.fMsh["y"]

        A.setValue(node, node, 1.+self.rx[node]/2.+self.ry[node]/2.)

        # both z and r direction are the neumann b.c.
        if (self.fBound["x"]=="Neumann" or self.fBound["x"]=="N") and \
           (self.fBound["y"]=="Neumann" or self.fBound["y"]=="N"):
            b[node] += 0.
            if node==0:
                A.setValue(node, node+nx, -self.ry[node]/2.)
                A.setValue(node,  node+1, -self.rx[node]/2.)
            elif node==nx-1:
                A.setValue(node, node+nx, -self.ry[node]/2.)
                A.setValue(node,  node-1, -self.rx[node]/2.)
            elif node==nx*(ny-1):
                A.setValue(node, node-nx, -self.ry[node]/2.)
                A.setValue(node,  node+1, -self.rx[node]/2.)
            elif node==nx*ny-1:
                A.setValue(node, node-nx, -self.ry[node]/2.)
                A.setValue(node,  node-1, -self.rx[node]/2.)

        # z: dirichlet , r: neumann b.c.
        elif (self.fBound["x"]=="Dirichlet" or self.fBound["x"]=="D") and \
             (self.fBound["y"]=="Neumann"   or self.fBound["y"]=="N"):
            if node==0:
                A.setValue(node, node+nx, -self.ry[node]/2.)
                A.setValue(node,  node+1, -self.rx[node]/4.)
                b.setValue(node, b.getValues(node) + self.rx[node]/4. * self.fBc["x"])
                #b.setValue(node, b.getValues(node) + self.ry[node]/4.*self.fBc["x"])
            elif node==nx-1:
                A.setValue(node, node+nx, -self.ry[node]/2.)
                A.setValue(node,  node-1, -self.rx[node]/4.)
                b.setValue(node, b.getValues(node) + self.rx[node]/4. * self.fBc["x"])
                #b.setValue(node, b.getValues(node) + self.ry[node]/4.*self.fBc["x"])
            elif node==nx*(ny-1):
                A.setValue(node, node-nx, -self.ry[node]/2.)
                A.setValue(node,  node+1, -self.rx[node]/4.)
                b.setValue(node, b.getValues(node) + self.rx[node]/4. * self.fBc["x"])
                #b.setValue(node, b.getValues(node) + self.ry[node]/4.*self.fBc["x"])
            elif node==nx*ny-1:
                A.setValue(node, node-nx, -self.ry[node]/2.)
                A.setValue(node,  node-1, -self.rx[node]/4.)
                b.setValue(node, b.getValues(node) + self.rx[node]/4. * self.fBc["x"])
                #b.setValue(node, b.getValues(node) + self.ry[node]/4.*self.fBc["x"])

        # z: neumann b.c , r: dirichlet b.c
        elif (self.fBound["x"]=="Neumann"   or self.fBound["x"]=="N") and \
             (self.fBound["y"]=="Dirichlet" or self.fBound["y"]=="D"):
            if node==0:
                A.setValue(node, node+nx, -self.ry[node]/4.)
                A.setValue(node,  node+1, -self.rx[node]/2.)
                b.setValue(node, b.getValues(node) + self.ry[node]/4.*self.fBc["y"])
                #b.setValue(node, b.getValues(node) + self.rx[node]/4.*self.fBc["y"])
            elif node==nx-1:
                A.setValue(node, node+nx, -self.ry[node]/4.)
                A.setValue(node,  node-1, -self.rx[node]/2.)
                b.setValue(node, b.getValues(node) + self.ry[node]/4. * self.fBc["y"])
                #b.setValue(node, b.getValues(node) + self.rx[node]/4. * self.fBc["y"])
            elif node==nx*(ny-1):
                A.setValue(node, node-nx, -self.ry[node]/4.)
                A.setValue(node,  node+1, -self.rx[node]/2.)
                b.setValue(node, b.getValues(node) + self.ry[node]/4. * self.fBc["y"])
                #b.setValue(node, b.getValues(node) + self.rx[node]/4. * self.fBc["y"])
            elif node==nx*ny-1:
                A.setValue(node, node-nx, -self.ry[node]/4.)
                A.setValue(node,  node-1, -self.rx[node]/2.)
                b.setValue(node, b.getValues(node) + self.ry[node]/4. * self.fBc["y"])
                #b.setValue(node, b.getValues(node) + self.rx[node]/4. * self.fBc["y"])

        # z: dirichlet b.c. , r: dirichlet b.c.
        elif (self.fBound["x"]=="Dirichlet" or self.fBound["x"]=="D")  and \
             (self.fBound["y"]=="Dirichlet" or self.fBound["y"]=="D"):
            if node==0:
                A.setValue(node, node+nx, -self.ry[node]/4.)
                A.setValue(node,  node+1, -self.rx[node]/4.)
                b.setValue(node, b.getValues(node) + self.ry[node]/4.*self.fBc["y"] + self.rx[node]/4.*self.fBc["x"])
            elif node==nx-1:
                A.setValue(node, node+nx, -self.ry[node]/4.)
                A.setValue(node,  node-1, -self.rx[node]/4.)
                b.setValue(node, b.getValues(node) + self.ry[node]/4.*self.fBc["y"] + self.rx[node]/4.*self.fBc["x"])
            elif node==nx*(ny-1):
                A.setValue(node, node-nx, -self.ry[node]/4.)
                A.setValue(node,  node+1, -self.rx[node]/4.)
                b.setValue(node, b.getValues(node) + self.ry[node]/4.*self.fBc["y"] + self.rx[node]/4.*self.fBc["x"])
            elif node==nx*ny-1:
                A.setValue(node, node-nx, -self.ry[node]/4.)
                A.setValue(node,  node-1, -self.rx[node]/4.)
                b.setValue(node, b.getValues(node) + self.ry[node]/4.*self.fBc["y"] + self.rx[node]/4.*self.fBc["x"])



    def MakeMatrix(self, A, b):
        node = self.fNode
        nx = self.fMsh["x"]
        ny = self.fMsh["y"]

        A.setValue(node, node, 1.+2.*self.rx[node]+2.*self.ry[node])

        # neumann b.c.
        if self.fBound["x"]=="N" or self.fBound["x"]=="Neumann":
            # bottom
            if node-nx < 0:
                A.setValue(node, node+nx, -2.*self.ry[node])
                A.setValue(node,  node+1, -self.rx[node])
                A.setValue(node,  node-1, -self.rx[node])
            # top
            elif node+nx > nx*ny-1:
                A.setValue(node, node-nx, -2.*self.ry[node])
                A.setValue(node,  node+1, -self.rx[node])
                A.setValue(node,  node-1, -self.rx[node])

            for i in range(1,ny-1):
                # left
                if node == self.id(i,0):
                    A.setValue(node, node-nx, -self.ry[node])
                    A.setValue(node, node+nx, -self.ry[node])
                    A.setValue(node,  node+1, -2.*self.rx[node])
                # right
                elif node == self.id(i,nx-1):
                    A.setValue(node, node-nx, -self.ry[node])
                    A.setValue(node, node+nx, -self.ry[node])
                    A.setValue(node,  node-1, -2.*self.rx[node])
        # dirichlet
        elif self.fBound["x"]=="D" or self.fBound["x"]=="Dirichlet":
            # bottom
            if node-nx < 0:
                A.setValue(node, node+nx, -self.ry[node])
                A.setValue(node,  node+1, -self.rx[node])
                A.setValue(node,  node-1, -self.rx[node])
                b.setValue(node, b.getValues(node) + self.ry[node] * self.fBc["x"])
                #b.setValue(node, b.getValues(node) + self.rx[node] * self.fBc["x"])
            # top
            elif node+nx > nx*ny-1:
                A.setValue(node, node-nx, -self.ry[node])
                A.setValue(node,  node+1, -self.rx[node])
                A.setValue(node,  node-1, -self.rx[node])
                b.setValue(node, b.getValues(node) + self.ry[node] * self.fBc["x"])
                #b.setValue(node, b.getValues(node) + self.rx[node] * self.fBc["x"])

            for i in range(1,ny-1):
                # left
                if node == self.id(i,0):
                    A.setValue(node, node-nx, -self.ry[node])
                    A.setValue(node, node+nx, -self.ry[node])
                    A.setValue(node,  node+1, -self.rx[node])
                    b.setValue(node, b.getValues(node) + self.rx[node] * self.fBc["x"])
                    #b.setValue(node, b.getValues(node) + self.ry[node]*self.fBc["x"])
                # right
                elif node == self.id(i,nx-1):
                    A.setValue(node, node-nx, -self.ry[node])
                    A.setValue(node, node+nx, -self.ry[node])
                    A.setValue(node,  node-1, -self.rx[node])
                    b.setValue(node, b.getValues(node) + self.rx[node] * self.fBc["x"])
                    #b.setValue(node, b.getValues(node) + self.ry[node]*self.fBc["x"])

