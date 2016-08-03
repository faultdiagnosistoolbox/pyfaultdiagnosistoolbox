class VarIdGen:
    xlabel = 'x'
    zlabel = 'z'
    flabel = 'f'
    elabel = 'e'

    def __init__(self):
        self.x = 1
        self.z = 1
        self.f = 1
        self.e = 1
        
    def setXLabel(self,l):
        self.xlabel = l

    def setFLabel(self,l):
        self.flabel = l

    def setZLabel(self,l):
        self.zlabel = l

    def setELabel(self,l):
        self.elabel = l

    def reset(self):
        self.resetX()
        self.resetF()
        self.resetZ()
        self.resetE()

    def state(self):
        return (self.x, self.z, self.f, self.e)
    
    def setState( self, s ):
        self.x = s[0];
        self.z = s[1];
        self.f = s[2];
        self.e = s[3];

    def resetX(self):
        self.x = 1
    def resetF(self):
        self.f = 1
    def resetZ(self):
        self.z = 1
    def resetE(self):
        self.e = 1

    def NewX(self):
        r = self.xlabel + str(self.x)
        self.x = self.x + 1
        return r
    def NewZ(self):
        r = self.zlabel + str(self.z)
        self.z = self.z + 1
        return r
    def NewF(self):
        r = self.flabel + str(self.f)
        self.f = self.f + 1
        return r
    def NewE(self):
        r = self.elabel + str(self.e)
        self.e = self.e + 1
        return r
