
from IMakeMaterial import IMakeMaterial
from IConstruction import IConstruction
from IBoundary     import IBoundary
from IPostprocess  import IPostprocess

from petsc4py import PETSc

import matplotlib.pyplot as plt

mm = 1.0e-3
cm = 1.0e-2
m  = 1.

class IConstructCOMET:

    def __init__(self, nx, ny):
        self.fMsh  = {"x":nx, "y":ny}
        self.fPrep = IMakeMaterial(nx, ny)
        self.dt = 1.
        self.tf = 60.
        self.T0 = 4.5


    def MakeStrip(self, r0, lr):
        l = 1.35*m
        #node = self.fPrep.MakeRectangle(0.,r0,l+40.*cm,lr)
        node = self.fPrep.MakeRectangle(0.,r0,l,lr)
        for i in node:
            self.fPrep.SetMaterial(i, rho=2700., C=0.1, kx=500., ky=1., q=0.03*2700, name="Strip")


    def MakeConductor(self,r0):
        l = 1.35*m
        #node = self.fPrep.MakeRectangle(20.*cm,r0,l,15.*mm)
        node = self.fPrep.MakeRectangle(0.*cm,r0,l,15.*mm)
        for i in node:
            self.fPrep.SetMaterial(i, rho=4000., C=0.2, kx=0.011, ky=0.01, q=0.03*4000., name="Conductor")


    def SetGeometry(self, plot=True):
        r0 = 673.*mm
        r1 = r0 + 1.*mm*9 + 15.*mm*9
        l  = 1.35 * m
        lr = 1.*mm
        lcdt = 15.*mm
        #self.fPrep.SetSpace(0.*cm, r0, l+40.*cm, r1-r0)
        self.fPrep.SetSpace(0.*cm, r0, l, r1-r0)

        # 1st layer:
        self.MakeStrip(r0, lr)
        self.MakeConductor(r0+lr)
        # 2nd layer:
        self.MakeStrip(r0+lr+lcdt, lr)
        self.MakeConductor(r0+lr*2+lcdt)
        # 3rd layer:
        self.MakeStrip(r0+lr*2+lcdt*2, lr)
        self.MakeConductor(r0+lr*3+lcdt*2)
        # 4th layer:
        self.MakeStrip(r0+lr*3+lcdt*3, lr)
        self.MakeConductor(r0+lr*4+lcdt*3)
        # 5th layer:
        self.MakeStrip(r0+lr*4+lcdt*4, lr)
        self.MakeConductor(r0+lr*5+lcdt*4)
        # 6th layer:
        self.MakeStrip(r0+lr*5+lcdt*5, lr)
        self.MakeConductor(r0+lr*6+lcdt*5)
        # 7th layer:
        self.MakeStrip(r0+lr*6+lcdt*6, lr)
        self.MakeConductor(r0+lr*7+lcdt*6)
        # 8th layer:
        self.MakeStrip(r0+lr*7+lcdt*7, lr)
        self.MakeConductor(r0+lr*8+lcdt*7)
        # 9th layer:
        self.MakeStrip(r0+lr*8+lcdt*8, lr)
        self.MakeConductor(r0+lr*9+lcdt*8)

        self.fPrep.Plot()
        #node = self.fPrep.GetNodeInfo()["id"]
        #for i in node:
            #print self.fPrep.GetMaterial()["kx"][i], self.fPrep.GetMaterial()["ky"][i]


    def SetTime(self, tf=60., dt=1.):
        self.tf = tf
        self.dt = dt


    def SetInitialTemp(self, T0):
        self.T0 = T0


    def SetMethod(self, name, A):
        ksp = PETSc.KSP()
        ksp.create(PETSc.COMM_WORLD)
        ksp.setOperators(A)
        ksp.getPC().setType(name)
        ksp.setFromOptions()

        return ksp


    def Solve(self, A, T, b):
        begin, end = T.getOwnershipRange()

        q = PETSc.Vec()
        q.create(PETSc.COMM_WORLD)
        q.setSizes(end)
        q.setUp()

        for i in range(begin,end):
            q.setValue(i, b.getValues(i))

        ksp = self.SetMethod("lu", A)

        print "***************************"
        print "start time loop!"
        print "***************************"

        time = 0.

        while (time<=self.tf):

            for i in range(begin, end):
                preT = T.getValues(i)
                heat = q.getValues(i)
                b.setValue(i, preT + heat)

            ksp.solve(b, T)

            print time, T.getValues(5), ksp.getIterationNumber(), ksp.getResidualNorm()

            time += self.dt

        ksp.view()
        ksp.destroy()


    def SetTwoSide(self, bound, A, b, lmin, lmax):
        line = self.fPrep.GetBoundaryLine("Left",lmin=lmin,lmax=lmax)
        for i in line:
            bound.SetBoundary(i, xb="D", xbc=self.T0)
            bound.MakeMatrix(A,b)

        line = self.fPrep.GetBoundaryLine("Right",lmin=lmin,lmax=lmax)
        for i in line:
            bound.SetBoundary(i, xb="D", xbc=self.T0)
            bound.MakeMatrix(A,b)


    def Run(self):
        cst = IConstruction()
        cst.SetTimeStep(self.dt)
        cst.SetInitialTemp(self.T0)
        cst.SetBoundaryTemp(self.T0)
        cst.SetMaterial(self.fPrep)
        A = cst.MakeMatrix()
        b = cst.MakeVector()

        # set boundary condition
        r0   = 673.*mm
        lAl  = 1.*mm
        lCdt = 15.*mm
        bdy = IBoundary(self.fMsh["x"], self.fMsh["y"])
        # get vertex node number
        vet = self.fPrep.GetVertexBoundary()
        # get material property
        r = cst.GetParameters()
        bdy.SetParameters(r["rx"], r["ry"])

        bdy.SetBoundary(vet[0], xb="D", xbc=self.T0, yb="N", ybc=0.)
        bdy.MakeVertexMat(A,b)
        bdy.SetBoundary(vet[1], xb="D", xbc=self.T0, yb="N", ybc=0.)
        bdy.MakeVertexMat(A,b)
        bdy.SetBoundary(vet[2], xb="N", xbc=0., yb="D", ybc=self.T0)
        bdy.MakeVertexMat(A,b)
        bdy.SetBoundary(vet[3], xb="N", xbc=0., yb="D", ybc=self.T0)
        bdy.MakeVertexMat(A,b)

        # set line boundary
        # get boundary line
        line = self.fPrep.GetBoundaryLine("Top")
        for i in line:
            bdy.SetBoundary(i, xb="D", xbc=self.T0)
            bdy.MakeMatrix(A,b)

        line = self.fPrep.GetBoundaryLine("Bottom")
        for i in line:
            bdy.SetBoundary(i, xb="N", xbc=0.)
            bdy.MakeMatrix(A,b)

        # left
        line = self.fPrep.GetBoundaryLine("Left")
        for i in line:
            bdy.SetBoundary(i, xb="N", xbc=0.)
            bdy.MakeMatrix(A,b)

        # right
        line = self.fPrep.GetBoundaryLine("Right")
        for i in line:
            bdy.SetBoundary(i, xb="N", xbc=0.)
            bdy.MakeMatrix(A,b)

        # two side
        self.SetTwoSide(bdy, A, b, r0, r0+lAl)
        self.SetTwoSide(bdy, A, b, r0+lAl+lCdt, r0+lAl*2+lCdt)
        self.SetTwoSide(bdy, A, b, r0+lAl*2+lCdt*2, r0+lAl*3+lCdt*2)
        self.SetTwoSide(bdy, A, b, r0+lAl*3+lCdt*3, r0+lAl*4+lCdt*3)
        self.SetTwoSide(bdy, A, b, r0+lAl*4+lCdt*4, r0+lAl*5+lCdt*4)
        self.SetTwoSide(bdy, A, b, r0+lAl*5+lCdt*5, r0+lAl*6+lCdt*5)
        self.SetTwoSide(bdy, A, b, r0+lAl*6+lCdt*6, r0+lAl*7+lCdt*6)
        self.SetTwoSide(bdy, A, b, r0+lAl*7+lCdt*7, r0+lAl*8+lCdt*7)
        self.SetTwoSide(bdy, A, b, r0+lAl*8+lCdt*8, r0+lAl*9+lCdt*8)

        A.assemble()
        b.assemble()

        # create vector for initial temperature
        T = PETSc.Vec()
        T.create(PETSc.COMM_WORLD)
        T.setSizes(self.fMsh["x"]*self.fMsh["y"])
        T.setUp()

        begin, end = T.getOwnershipRange()

        for node in range(begin,end):
            T.setValue(node, self.T0)

        T.assemble()

        self.Solve(A,T,b)

        return T


    def Postprocess(self, T):
        post = IPostprocess(self.fPrep)
        post.SetTimeStep(self.dt)
        post.SetResult(T)

        print "min: %.2f K, max: %.2f K" %(post.getmin(), post.getmax())
        post.PlotMatTogether(save=True)
        #plt.subplot_tool()
        plt.show()


