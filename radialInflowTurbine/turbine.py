## Creating a class for radial inflow turbine design
import matplotlib.pyplot as plt
import numpy as np
import math
import pprint

import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

import streamProps

class Format:
    end = '\033[0m'
    underline = '\033[4m'

class Turbine:

    def __init__(self, fluid, T_in, P_in, P_out, massFlow):

        self.fluid    = fluid
        self.T_in     = T_in        #K
        self.P_in     = P_in*1e5    #Pa
        self.P_out    = P_out*1e5   #Pa
        self.massFlow = massFlow    #kg/s

        self.E_ratio  = self.P_out / self.P_in

        #Inlet thermo props
        self.inlet    = streamProps.ThermoState("TP", self.fluid, self.T_in, self.P_in)
        self.inletV   = massFlow / self.inlet.D   #m3/s

        #Isentropic compressor
        self.isen     = self.Isentropic("PS", self.inlet, self.P_out, self.massFlow)

    def initialSpec(self):

        print(Format.underline + "TURBO-EXPANDER INITIAL SPECIFICATION INITIAL SPECIFICATION" + Format.end)
        print("Fluid: " + self.fluid[9:])
        print("Mass flow [kg/s] : " + str(round(self.massFlow,4)) + "\n")
        print(Format.underline + "Inlet Condition" + Format.end)
        self.inlet.printState()
        print(Format.underline + "Isentropic Outlet Condition" + Format.end)
        self.isen.outlet.printState()
        self.isen.printPower()

    def noOfStages(self, u1_max):

        #Maximum tip velocity is a constraint that needs to be defined by the designer
        self.u1_max   = u1_max                            #m/s
        self.maxDH    = (u1_max/self.isen.nu_s)**2 / 2    #J/kg

        self.minZ     = round(self.isen.DH / self.maxDH)
        if self.minZ == 0:
            self.minZ = 1

        self.maxZ     = self.minZ+2

        print(Format.underline + "Mininum and Maximum Number of Stage" + Format.end)
        print("Wheel Tip Velocity [m/s]  = ", self.u1_max)
        print("\u0394h_max/stage [kJ/kg] = ", round(self.maxDH/1e3,3))
        print("Min. number of stage = ", self.minZ)
        print("Max. number of stage = ", self.maxZ)

    def rotationalSpeedRange(self):

        self.optOmega_s = np.array([0.4,0.8])

        outletV       = self.massFlow / self.isen.outlet.D    #m3/s
        firstStage    = streamProps.ThermoState("HS", self.fluid, self.inlet.H-(self.isen.DH/self.maxZ), self.inlet.S)
        firstStageV   = self.massFlow / firstStage.D         #m3/s

        n_minZ        = self.optOmega_s * (60 / (2*math.pi)) * ((self.isen.DH/self.minZ)**0.75) / outletV**0.5 #rpm
        n_maxZ        = self.optOmega_s * (60 / (2*math.pi)) * ((self.isen.DH/self.maxZ)**0.75) / firstStageV**0.5 #rpm

        nMin          = np.min(np.append(n_minZ,n_maxZ))      #rpm
        nMax          = np.max(np.append(n_minZ,n_maxZ))      #rpm

        self.n_range  = np.append(nMin,nMax)

        print(Format.underline + "Minimum & Maximum Rotational Speed" + Format.end)
        print("Min. RPM = ", round(nMin,3))
        print("Max. RPM = ", round(nMax,3))

    def nzSelectionPlot(self):

        z             = np.arange(self.minZ, self.maxZ+1,1)
        ns            = np.linspace(min(self.n_range), max(self.n_range), 4)

        omegaFirst    = np.empty(shape=(len(ns),len(z)))
        omegaLast     = np.empty(shape=(len(ns),len(z)))

        print(Format.underline + "Rotational Speed Selection Plot \n" + Format.end)
        print("Specific speed range for radial compressor is between 0.4-0.8")

        for i in range(len(ns)):

            for j in range(len(z)):

                firstStage     = streamProps.ThermoState("HS", self.fluid, self.inlet.H - (self.isen.DH/z[j]), self.inlet.S)
                firstStageV    = self.massFlow / firstStage.D
                lastStageV     = self.massFlow / self.isen.outlet.D    #m3/s

                omegaFirst[i,j] = (2*math.pi*ns[i] / 60) * (firstStageV**0.5 / (self.isen.DH/z[j])**0.75)
                omegaLast[i,j]  = (2*math.pi*ns[i] / 60) * (lastStageV**0.5 / (self.isen.DH/z[j])**0.75)

            firstStage = plt.plot(z, omegaFirst[i,:], label=str(round(ns[i]))+" RPM")
            plt.plot(z, omegaLast[i,:], "--", color = firstStage[0].get_color())

        plt.fill_between(z, self.optOmega_s[0], self.optOmega_s[1], alpha=0.1)
        z_int = range(math.floor(min(z)), math.ceil(max(z))+1)
        plt.xticks(z_int)
        plt.xlabel("Number of Stages", fontsize=14)
        plt.ylabel("Specific Speed", fontsize=14)
        plt.ylim(0.2,1.4)
        plt.grid(which='major', color='#DDDDDD', linewidth=1)
        plt.grid(which='minor', color='#EEEEEE', linewidth=0.5)
        plt.legend(fontsize=12)
        plt.yticks(fontsize=13)
        plt.xticks(fontsize=13)
        plt.tight_layout()
        plt.savefig('nz.png')
        plt.show()

    def stagesCalculation(self, n, z, warning=False):

        self.n          = n #rpm
        self.z          = z #number of stage
        stageInletV     = self.inletV
        DHperStage      = self.isen.DH/self.z

        self.isenStages = np.array([])
        self.stages     = np.array([])
        self.stageInV   = np.array([])

        #Assumptions
        self.deltah     = 0.185
        self.alpha2     = 0

        self.thermodynamicDict = {}
        self.kinematicDict     = {}
        self.geometricDict     = {}

        for i in range(z):

            self.stageInV    = np.append(self.stageInV, stageInletV)

            self.isenStages  = np.append(self.isenStages,self.Isentropic("HS", self.inlet, \
                                                                    self.inlet.H - DHperStage*(i+1), self.massFlow))
            stageOutletP     = self.isenStages[i].outlet.P
            stageOutletV_Is  = self.massFlow / self.isenStages[i].outlet.D

            if i == 0:
                stageInlet              = self.inlet

                omega_s, nu_s, alpha1,  \
                isPsi, u1, D1, D0,      \
                eta_ts, eta_is, eta_S,  \
                Psi, c1m, b1, xi, phi,  \
                nextH                   = self.designParameters(n, stageInletV, stageOutletV_Is, DHperStage, stageInlet.H, stageInlet.S, stageOutletP)

            else:
                stageInlet              = self.stages[i-1]

                omega_s, nu_s, alpha1,  \
                isPsi, u1, D1, D0,      \
                eta_ts, eta_is, eta_S,  \
                Psi, c1m, b1, xi, phi,  \
                nextH                   = self.designParameters(n, stageInletV, stageOutletV_Is, DHperStage, stageInlet.H, stageInlet.S, stageOutletP)

            if u1 > self.u1_max:
                print("Tip velocity at stage", i+1, ", exceeds the maximum value: ", \
                      u1, "m/s (>", self.u1_max, "m/s)")
                break

            ER              = stageInlet.P/stageOutletP

            self.stages     = np.append(self.stages, streamProps.ThermoState("PH", self.fluid, stageOutletP, nextH))

            stageInletV     = self.massFlow / self.stages[i].D
            stageDH         = stageInlet.H - self.stages[i].H
            stagePower      = stageDH*self.massFlow

            #omega_s         = self.specificSpeed(n, stageInletV, DHperStage)
            deltat          = self.rotorTipDiameterRatio(omega_s, stageInletV, phi, u1, D1)
            b2              = self.bladeHeightRotorOutlet(D1, deltat)
            R               = self.degreeReaction(Psi, phi, xi, deltat)
            beta1, beta2    = self.relativeFlowAngles(deltat, phi, Psi, xi)

            if R > (1+(phi**2)/(2*Psi)-(Psi/2)):

                print("Degree of reaction at stage", i+1, ", reached the maximum value", \
                     R, " (>", (1+(phi**2)/(2*Psi)-(Psi/2)))
                break

            c1m, c1u, c2m,  \
            c2u, w1u, w2u,  \
            c1, c2, c0, w1, \
            w2, u2          = self.velocityKinematics(u1, phi, xi, alpha1, beta1, beta2, deltat)

            statorOutlet    = self.statorOutlet(stageInlet.H, stageInlet.S, c0, c1, eta_S)

            Ma0, Ma1S, Ma1R, Ma2 = self.MachNumbers(stageInlet.A, statorOutlet.A, self.stages[i].A, c0, c1, w1, w2)

            self.thermodynamicDict['Stage '+str(i+1)] = {
                "Stage inlet volumetric [m3/h]"       : round(self.stageInV[i]*3600,2),
                "Stator outlet pressure [bar]"        : round(statorOutlet.P/1e5,2),
                "Stator outlet temperature [\u2103]"  : round(statorOutlet.T-273.15,2),
                "Rotor outlet pressure [bar]"         : round(self.stages[i].P/1e5,2),
                "Rotor outlet temperature [\u2103]"   : round(self.stages[i].T-273.15,2),
                "Total-to-static efficiency [%]"      : round(eta_ts*100,2),
                "Isentropic efficiency [%]"           : round(eta_is*100,2)
            }

            self.geometricDict['Stage '+str(i+1)] = {
                "Rotor tip diameter ratio"          : round(deltat,3),
                "Rotor hub diameter ratio"          : round(self.deltah,3),
                "Stator inlet diameter [cm]"        : round(D0*100,2),
                "Rotor inlet diameter [cm]"         : round(D1*100,2),
                "Rotor outlet diameter [cm]"        : round(D1*deltat*100,2),
                "Blade height at rotor inlet  [cm]" : round(b1*100,2),
                "Blade height at rotor outlet  [cm]": round(b2*100,2)
            }

            self.kinematicDict['Stage '+str(i+1)] = {
                "Velocities [m/s]" : {
                    "Absolute meridional rotor inlet"    : round(c1m,1),
                    "Absolute tangential rotor inlet"    : round(c1u,1),
                    "Absolute rotor inlet"               : round(c1,1),
                    "Absolute tangential rotor outlet"   : round(c2m,1),
                    "Absolute meridional rotor outlet"   : round(c2m,1),
                    "Absolute tangential rotor outlet"   : round(c2u,1),
                    "Absolute rotor outlet"              : round(c2,1),
                    "Absolute stator inlet"              : round(c0,1),
                    "Relative tangential rotor inlet"    : round(w1u,1),
                    "Relative rotor inlet"               : round(w1,1),
                    "Relative tangential rotor outlet"   : round(w2u,1),
                    "Relative rotor outlet"              : round(w2,1),
                    "Blade speed rotor inlet"            : round(u1,1),
                    "Blade speed rotor outlet"           : round(u2,1)
                },
                "Mach Numbers" : {
                    "Ma at stator inlet"                 : round(Ma0,3),
                    "Ma at stator outlet"                : round(Ma1S,3),
                    "Ma at rotor inlet"                  : round(Ma1R,3),
                    "Ma at rotor outlet"                 : round(Ma2,3)
                },
                "Non-dimensional" : {
                    "Specific speed"                     : round(omega_s,3),
                    "Flow coefficient"                   : round(phi,3),
                    "Work coefficient"                   : round(Psi,3),
                    "Isentropic velocity ratio"          : round(math.sqrt(1/(2*Psi)),3),
                    "Isentropic work coefficient"        : round(isPsi,3),
                    "Rotor meridional velocity ratio"    : round(xi,3),
                    "Degree of reaction"                 : round(R,3)
                },
                "Flow angles [\u00b0]" : {
                    "Absolute rotor inlet"               : round(alpha1,2),
                    "Absolute rotor outlet"              : round(self.alpha2,2),
                    "Relative rotor inlet"               : round(beta1,2),
                    "Relative rotor outlet"              : round(beta2,2)
                }
            }

            print(Format.underline + 'Stage '+str(i+1) + Format.end)
            print("Delta enthalpy [kJ/kg]:" + str(round(stageDH/1000,2)))
            print("Expansion ratio:" + str(round(ER,2)))
            print("Power [kW]:" + str(round(stagePower/1000,2)))
            print("\033[1mThermodynamic Parameters:\033[0m")
            pprint.pprint(self.thermodynamicDict['Stage '+str(i+1)], sort_dicts=False)
            print("\033[1mGeometric Parameters:\033[0m")
            pprint.pprint(self.geometricDict['Stage '+str(i+1)], sort_dicts=False)
            print("\033[1mKinematic Parameters:\033[0m")
            pprint.pprint(self.kinematicDict['Stage '+str(i+1)], sort_dicts=False)

            if warning == True:
                print("\033[1mWarnings:\033[0m")
                if (omega_s < 0.4) or (omega_s > 0.8):
                    print("Non-optimal specific speed (0.4-0.8): ", round(omega_s,3))
                if (phi < 0.2) or (phi > 0.3):
                    print("Non-optimal flow coefficient (0.2-0.3): ", round(phi,3))
                if (Psi < 0.8) or (Psi > 1.0):
                    print("Non-optimal work coefficient (0.8-1.0): ", round(Psi,3))
                if (R < 0.45) or (R > 0.65):
                    print("Non-optimal degree of reaction (0.45-0.65): ", round(R,3))
                if (xi < 0.65) or (xi > 1.00):
                    print("Non-optimal degree of reaction (0.65-1.00): ", round(xi,3))
                if (deltat > 0.7):
                    print("Non-optimal rotor tip diameter ratio (<0.7): ", round(deltat,3))
                if ((self.deltah/deltat) < 0.4):
                    print("Non-optimal rotor hub/tip ratio (>0.4): ", round(self.deltah/deltat,3))
                if (alpha1 < 66) or (alpha1 > 78):
                    print("Non-optimal rotor inlet absolute flow angle (66\u00b0-78\u00b0): ", round(alpha1,2))
                if (beta1 < 20) or (beta1 > 60):
                    print("Non-optimal rotor inlet relative flow angle (20\u00b0-40\u00b0): ", round(beta1,2))

    def turbDictionary(self):

        self.outlet     = streamProps.ThermoState("PH", self.fluid, self.stages[-1].P, self.stages[-1].H)

        self.DH         = self.inlet.H - self.outlet.H
        self.Power      = self.DH * self.massFlow
        self.turbEta_is = self.Power / self.isen.Power
        self.turbEta_p  = self.polytropicEfficiency()

        self.turbineDict = {
            "Fluid"                   : self.fluid[9:],
            "Mass flow [kg/s]"        : self.massFlow,
            "Turbine type"            : "Radial inflow turbine",
            "Rotational speed [RPM]"  : self.n,
            "Number of stages"        : self.z,
            "Max. tip velocity [m/s]" : self.u1_max,
            "Inlet pressure [bar]"    : round(self.inlet.P/1e5,2),
            "Inlet temperature [C]"   : round(self.inlet.T-273.15,2),
            "Inlet volumetric [m3/h]" : round(self.inletV*3600,2),
            "Delta enthalpy [kJ/kg]"  : round(self.DH/1000,2),
            "Expansion ratio"         : round(self.E_ratio,2),
            "Turbine power [kW]"      : round(self.Power/1000,2),
            "Turb. Isentropic eff. [%]" : round(self.turbEta_is*100,2),
            "Turb. Polytropic eff. [%]" : round(self.turbEta_p*100,2)
        }

        print(Format.underline + "RADIAL INFLOW TURBINE FINAL SPECIFICATION" + Format.end)
        pprint.pprint(self.turbineDict, sort_dicts=False)

    def MachNumbers(self, stageInletA, statorOutletStateA, stageOutletA, c0, c1, w1, w2):

        Ma0       = c0 / stageInletA
        Ma1S      = c1 / statorOutletStateA
        Ma1R      = w1 / statorOutletStateA
        Ma2       = w2 / stageOutletA

        return Ma0, Ma1S, Ma1R, Ma2

    def statorOutlet(self, stageInletH, stageInletS, c0, c1, eta_S):

        statorInletTotalH    = stageInletH + (c0**2)/2
        statorOutletTotalH   = statorInletTotalH
        statorOutletH        = statorOutletTotalH - (c1**2)/2
        statorOutletH_is     = statorInletTotalH - (statorInletTotalH-statorOutletH)/eta_S
        statorOutlet_is      = streamProps.ThermoState("HS", self.fluid, statorOutletH_is, stageInletS)
        statorOutlet         = streamProps.ThermoState("PH", self.fluid, statorOutlet_is.P, statorOutletH)

        return statorOutlet

    def velocityKinematics(self, u1, phi, xi, alpha1, beta1, beta2, deltat):

        tanAlpha1 = math.tan(alpha1*(math.pi/180))
        tanAlpha2 = math.tan(self.alpha2*(math.pi/180))
        tanBeta1  = math.tan(beta1*(math.pi/180))
        tanBeta2  = math.tan(beta2*(math.pi/180))

        c1m       = u1*xi*phi
        c1u       = c1m*tanAlpha1
        c2m       = u1*phi
        c2u       = c2m*tanAlpha2
        w1u       = u1*xi*phi*tanBeta1
        w2u       = u1*phi*tanBeta2
        c1        = u1*xi*phi*math.sqrt(1+(tanAlpha1**2))
        c2        = u1*phi*math.sqrt(1+(tanAlpha2**2))
        c0        = c2 #Assumption
        w1        = u1*xi*phi*math.sqrt(1+(tanBeta1**2))
        w2        = u1*phi*math.sqrt(1+(tanBeta2**2))
        u2        = u1*deltat

        return c1m, c1u, c2m, c2u, w1u, w2u, c1, c2, c0, w1, w2, u2

    def designParameters(self, n, inletV, outletV, DHperStage, inletH, inletS, stageOutletP):

        omega_s   = self.specificSpeed(n, outletV, DHperStage)
        nu_s      = 0.737*(omega_s**0.2)
        alpha1    = 79.2 - 14.2*(omega_s**2)
        isPsi     = 1/ (2*(nu_s**2))
        u1        = math.sqrt(DHperStage/isPsi)
        D1        = 60*u1 / (math.pi*n)
        eta_ts    = self.total2staticEfficiency(omega_s)

        tanAlpha1 = math.tan(alpha1*(math.pi/180))
        tanAlpha2 = math.tan(self.alpha2*(math.pi/180))

        eta_is    = 0.7 #initial guess
        statorOutletV = inletV #initial guess
        error     = 10  #initial value

        while error > 0.001:

            Psi         = isPsi*eta_is
            c1m         = u1*(Psi/tanAlpha1)
            b1          = statorOutletV/(math.pi*D1*c1m)
            xi          = (1 + 20*((b1/D1)**2))**-1
            phi         = c1m / (xi*u1)
            c1          = u1*xi*phi*math.sqrt(1+(tanAlpha1**2))
            c2          = u1*phi*math.sqrt(1+(tanAlpha2**2))
            c0          = c2 #Assumption

            eta_S       = eta_is #Assumption
            statorOutlet    = self.statorOutlet(inletH, inletS, c0, c1, eta_S)
            statorOutletV   = self.massFlow / statorOutlet.D

            isenExpansion   = streamProps.ThermoState("PS", self.fluid, stageOutletP, inletS)

            nextH       = inletH - ( ((c0**2)/2) + (inletH - isenExpansion.H) )*eta_ts

            guessEta_is = (inletH - nextH) / (inletH - isenExpansion.H)
            error       = abs(eta_is - guessEta_is)

            eta_is     += 0.001

        eta_S     = eta_is
        D0        = 1.25*D1*(1 + 4*(b1/D1)*math.cos(alpha1*(math.pi/180)))

        return omega_s, nu_s, alpha1, isPsi, u1, D1, D0, eta_ts, eta_is, \
                eta_S, Psi, c1m, b1, xi, phi, nextH

    def polytropicEfficiency(self):

        isobarH    = streamProps.ThermoState("PH", self.fluid, self.inlet.P, self.outlet.H)

        turbEta_p  = (self.inlet.S-isobarH.S) / (self.outlet.S - isobarH.S)

        return turbEta_p

    def specificSpeed(self, n, outletV, DHperStage):

        omega_s   = (2*math.pi*n / 60) * (outletV**0.5 / (DHperStage)**0.75)

        return omega_s

    def total2staticEfficiency(self, omega_s):

        w                        = omega_s-0.55
        eta_ts                   = 0.87 - 1.07*(w**2) - 0.5*(w**3)

        return eta_ts

    def rotorTipDiameterRatio(self, omega_s, outletV, phi, u1, D1):

        A         = (4*outletV) / (math.pi*phi*u1*(D1**2))
        deltat    = math.sqrt(self.deltah**2 + A)

        return deltat

    def bladeHeightRotorOutlet(self, D1, deltat):

        b2        = (D1/2) * (deltat - self.deltah)

        return b2

    def degreeReaction(self, Psi, phi, xi, deltat):

        tanAlpha2 = math.tan(self.alpha2*(math.pi/180))
        A         = (1-(xi**2)) + (tanAlpha2**2)*(1-(deltat**2))
        R         = 1 - Psi/2 + (phi**2)*A/(2*Psi) + phi*deltat*tanAlpha2

        return R

    def relativeFlowAngles(self, deltat, phi, Psi, xi):

        tanAlpha2 = math.tan(self.alpha2*(math.pi/180))
        A         = (Psi-1) / (phi*xi)
        B         = deltat*tanAlpha2 / xi
        beta1     = math.atan(A-B) * 180/math.pi
        beta2     = math.atan((deltat/phi) + tanAlpha2) * 180/math.pi

        return beta1, beta2

    class Isentropic:

        def __init__(self, pair, inlet, outlet_value, massFlow):

            #outlet isentropic thermo props
            if   pair == "PS":
                self.outlet =  streamProps.ThermoState("PS", inlet.fluid, outlet_value, inlet.S)

            elif pair == "HS":
                self.outlet =  streamProps.ThermoState("HS", inlet.fluid, outlet_value, inlet.S, inlet.fluid)

            #isentropic power
            self.DH     = inlet.H - self.outlet.H    #J/kg
            self.Power  = self.DH * massFlow         #W

            #isentropic velocity ratio
            self.nu_s    = 0.7

        def printPower(self):
            print("\u0394h_is [kJ/kg]  = ", round(self.DH/1e3,3))
            print("Power_is [kW]  =", round(self.Power/1e3,3))
