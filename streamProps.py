from CoolProp.CoolProp import PropsSI

class ThermoState():
    def __init__(self, pair, fluid, value1, value2, value3=None, name=None):
        self.fluid = fluid
        self.name  = name
        self.M     = value3 #kg/s
        if   pair == "TQ":
            self.T = value1 #K
            self.Q = value2 #K
            self.P = PropsSI("P", "T", self.T, "Q", value2, self.fluid) #Pa
            self.H = PropsSI("H", "T", self.T, "Q", value2, self.fluid) #J/kg
            self.S = PropsSI("S", "T", self.T, "Q", value2, self.fluid) #J/kgK
        elif pair == "TP":
            self.T = value1 #K
            self.P = value2 #Pa
            self.H = PropsSI("H", "T", self.T, "P", value2, self.fluid) #J/kg
            self.S = PropsSI("S", "T", self.T, "P", value2, self.fluid) #J/kgK
            self.Q = PropsSI("Q", "T", self.T, "P", value2, self.fluid)
        elif pair == "PS":
            self.P = value1 #Pa
            self.S = value2 #J/kgK
            self.H = PropsSI("H", "P", self.P, "S", self.S, self.fluid) #J/kg
            self.T = PropsSI("T", "P", self.P, "S", self.S, self.fluid) #K
            self.Q = PropsSI("Q", "P", self.P, "S", self.S, self.fluid)
        elif pair == "PH":
            self.P = value1 #Pa
            self.H = value2 #J/kg
            self.T = PropsSI("T", "P", self.P, "H", self.H, self.fluid) #K
            self.S = PropsSI("S", "P", self.P, "H", self.H, self.fluid) #J/kgK
            self.Q = PropsSI("Q", "P", self.P, "H", self.H, self.fluid)
        elif pair == "HS":
            self.H = value1 #J/kg
            self.S = value2 #J/kgK
            self.P = PropsSI("P", "H", self.H, "S", self.S, self.fluid) #bar
            self.T = PropsSI("T", "H", self.H, "S", self.S, self.fluid) #K
            self.Q = PropsSI("Q", "P", self.P, "S", self.S, self.fluid) 
        else:
            print("ERROR: pair = ", pair, " not implemented.")

        self.Z = PropsSI("Z", "P", self.P, "H", self.H, self.fluid)
        self.A = PropsSI("A", "P", self.P, "H", self.H, self.fluid) #m/s
        self.D = PropsSI("D", "P", self.P, "H", self.H, self.fluid) #kg/m3

    def printState(self):
        if self.name != None:
            print("Stream " + self.name)
        if self.M != None:
            print("M [kg/s]    = ", round(self.M,4))
        print("P [bar]     = ", round(self.P/1e5,3))
        print("T [Kelvin]  = ", round(self.T,2))
        print("H [kJ/kg]   = ", round(self.H/1e3,3))
        print("S [kJ/kg.K] = ", round(self.S/1e3,3))
        print("\n")
