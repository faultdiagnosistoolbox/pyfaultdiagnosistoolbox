"""Class for generating unique variable names."""


class VarIdGen:
    """Base class for variable name generator."""

    xlabel = 'x'
    zlabel = 'z'
    flabel = 'f'
    elabel = 'e'

    def __init__(self):
        """Initializer."""
        self.x = 1
        self.z = 1
        self.f = 1
        self.e = 1

    def setXLabel(self, lab):
        """Set x label."""
        self.xlabel = lab

    def setFLabel(self, lab):
        """Set f label."""
        self.flabel = lab

    def setZLabel(self, lab):
        """Set z label."""
        self.zlabel = lab

    def setELabel(self, lab):
        """Set e label."""
        self.elabel = lab

    def reset(self):
        """Reset all name generators."""
        self.resetX()
        self.resetF()
        self.resetZ()
        self.resetE()

    def state(self):
        """Return full state of name generators."""
        return self.x, self.z, self.f, self.e

    def setState(self, s):
        """Set full state of name generators."""
        self.x = s[0]
        self.z = s[1]
        self.f = s[2]
        self.e = s[3]

    def resetX(self):
        """Reset x name generator."""
        self.x = 1

    def resetF(self):
        """Reset f name generator."""
        self.f = 1

    def resetZ(self):
        """Reset z name generator."""
        self.z = 1

    def resetE(self):
        """Reset e name generator."""
        self.e = 1

    def NewX(self):
        """Generate new x name."""
        r = self.xlabel + str(self.x)
        self.x = self.x + 1
        return r

    def NewZ(self):
        """Generate new z name."""
        r = self.zlabel + str(self.z)
        self.z = self.z + 1
        return r

    def NewF(self):
        """Generate new f name."""
        r = self.flabel + str(self.f)
        self.f = self.f + 1
        return r

    def NewE(self):
        """Generate new e name."""
        r = self.elabel + str(self.e)
        self.e = self.e + 1
        return r
