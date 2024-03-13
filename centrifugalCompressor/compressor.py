## Creating a class for radial compressor design
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

class Compressor:

    def __init__(self, fluid, T_in, P_in, P_out, massFlow):

        self.fluid    = fluid
        self.T_in     = T_in        #K
        self.P_in     = P_in*1e5    #Pa
        self.P_out    = P_out*1e5   #Pa
        self.massFlow = massFlow    #kg/s

        self.C_ratio  = self.P_out / self.P_in

        #Inlet thermo props
        self.inlet    = streamProps.ThermoState("TP", self.fluid, self.T_in, self.P_in)
        self.inletV   = massFlow / self.inlet.D   #m3/s

        #Isentropic compressor
        self.isen     = self.Isentropic("PS", self.inlet, self.P_out, self.massFlow)

    def initialSpec(self):

        print(Format.underline + "RADIAL COMPRESSOR INITIAL SPECIFICATION" + Format.end)
        print("Fluid: " + self.fluid[9:])
        print("Mass flow [kg/s] : " + str(round(self.massFlow,4)) + "\n")
        print(Format.underline + "Inlet Condition" + Format.end)
        self.inlet.printState()
        print(Format.underline + "Isentropic Outlet Condition" + Format.end)
        self.isen.outlet.printState()
        self.isen.printPower()

    def noOfStages(self, u2_max):

        #Maximum tip velocity is a constraint that needs to be defined by the designer
        self.u2_max   = u2_max                    #m/s
        self.maxDH    = self.isen.Psi*(u2_max**2) #J/kg

        self.minZ     = round(self.isen.DH / self.maxDH)
        if self.minZ == 0:
            self.minZ = 1

        self.maxZ     = self.minZ+2

        print(Format.underline + "Mininum and Maximum Number of Stage" + Format.end)
        print("Wheel Tip Velocity [m/s]  = ", self.u2_max)
        print("\u0394h_max/stage [kJ/kg] = ", round(self.maxDH/1e3,3))
        print("Min. number of stage = ", self.minZ)
        print("Max. number of stage = ", self.maxZ)

    def rotationalSpeedRange(self):

        optOmega_s    = np.array([0.4,1.0])
        
        lastStage      = streamProps.ThermoState("HS", self.fluid, self.inlet.H + (self.isen.DH*(self.maxZ-1)/self.maxZ),\
                                                         self.inlet.S)
        lastStageV     = self.massFlow / lastStage.D
        
        n_minZ        = optOmega_s * (60 / (2*math.pi)) * ((self.isen.DH/self.minZ)**0.75) / self.inletV**0.5 #rpm
        n_maxZ        = optOmega_s * (60 / (2*math.pi)) * ((self.isen.DH/self.maxZ)**0.75) / lastStageV**0.5 #rpm

        nMin          = np.min(n_maxZ) #rpm
        nMax          = np.max(n_minZ) #rpm

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
        print("Specific speed range for radial compressor is between 0.4-1.0")

        for i in range(len(ns)):

            for j in range(len(z)):

                lastStage      = streamProps.ThermoState("HS", self.fluid, self.inlet.H + (self.isen.DH*(z[j]-1)/z[j]),\
                                                         self.inlet.S)
                lastStageV     = self.massFlow / lastStage.D

                omegaFirst[i,j] = (2*math.pi*ns[i] / 60) * (self.inletV**0.5 / (self.isen.DH/z[j])**0.75)
                omegaLast[i,j]  = (2*math.pi*ns[i] / 60) * (lastStageV**0.5 / (self.isen.DH/z[j])**0.75)

            firstStage = plt.plot(z, omegaFirst[i,:], label=str(round(ns[i]))+" RPM")
            plt.plot(z, omegaLast[i,:], "--", color = firstStage[0].get_color())

        plt.fill_between(z, 0.4, 1.0, alpha=0.1)
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
        self.deltah     = 0.35
        self.alpha1     = 0

        self.thermodynamicDict = {}
        self.kinematicDict     = {}
        self.geometricDict     = {}

        for i in range(z):

            self.stageInV    = np.append(self.stageInV, stageInletV)
            self.isenStages  = np.append(self.isenStages,self.Isentropic("HS", self.inlet, \
                                                                    self.inlet.H + DHperStage*(i+1), self.massFlow))
            stageOutletP     = self.isenStages[i].outlet.P

            if i == 0:
                stageInlet      = self.inlet

                omega_s, isPsi, \
                deltat, u2, D2, \
                b1, c1m, phi,   \
                phi_prime, D3,  \
                alpha2, eta_p,  \
                eta_R, eta_is,  \
                Psi, xi, R,     \
                beta1, beta2    = self.designParameters(n, stageInletV, DHperStage, stageInlet.H, stageInlet.P, stageOutletP)

                nextH           = stageInlet.H + DHperStage/eta_is

            else:
                stageInlet      = self.stages[i-1]

                omega_s, isPsi, \
                deltat, u2, D2, \
                b1, c1m, phi,   \
                phi_prime, D3,  \
                alpha2, eta_p,  \
                eta_R, eta_is,  \
                Psi, xi, R,     \
                beta1, beta2    = self.designParameters(n, stageInletV, DHperStage, stageInlet.H, stageInlet.P, stageOutletP)

                isenCompression = streamProps.ThermoState("PS", self.fluid, stageOutletP, self.stages[i-1].S)
                nextH           = stageInlet.H + (isenCompression.H - stageInlet.H) / eta_is

            if u2 > self.u2_max:
                print("Tip velocity at stage", i+1, ", exceeds the maximum value: ", \
                      u2, "m/s (>", self.u2_max, "m/s)")
                break

            if R > (1+(phi**2)/(2*Psi)-(Psi/2)):

                print("Degree of reaction at stage", i+1, ", reached the maximum value", \
                     R, " (>", (1+(phi**2)/(2*Psi)-(Psi/2)))
                break

            PR                   = stageOutletP/stageInlet.P

            c1m, c1u, c2m, c2u,\
            w1u, w2u, c1, c2,  \
            c3, w1, w2, u1       = self.velocityKinematics(u2, phi, xi, alpha2, beta1, beta2, deltat)

            rotorOutlet          = self.rotorOutlet(stageInlet.H, stageInlet.S, w1, u1, w2, u2, eta_R)
            b2                   = ((self.massFlow / rotorOutlet.D) / c2m) / (math.pi*D2)

            self.stages          = np.append(self.stages, streamProps.ThermoState("PH", self.fluid, stageOutletP, nextH))
            stageDH              = self.stages[i].H-stageInlet.H
            stagePower           = stageDH*self.massFlow

            Ma1, Ma2R, Ma2S, Ma3 = self.MachNumbers(stageInlet.A, rotorOutlet.A, self.stages[i].A, w1, w2, c2, c3)

            stageInletV          = self.massFlow / self.stages[i].D

            self.thermodynamicDict['Stage '+str(i+1)] = {
                "Stage inlet volumetric [m3/h]"       : round(self.stageInV[i]*3600,2),
                "Rotor outlet pressure [bar]"         : round(rotorOutlet.P/1e5,2),
                "Rotor outlet temperature [\u2103]"   : round(rotorOutlet.T-273.15,2),
                "Diffuser outlet pressure [bar]"      : round(self.stages[i].P/1e5,2),
                "Diffuser outlet temperature [\u2103]": round(self.stages[i].T-273.15,2),
                "Polytropic efficiency [%]"           : round(eta_p*100,2),
                "Isentropic efficiency [%]"           : round(eta_is*100,2)
            }

            self.geometricDict['Stage '+str(i+1)] = {
                "Rotor tip diameter ratio"          : round(deltat,3),
                "Rotor hub diameter ratio"          : round(self.deltah,3),
                "Rotor inlet diameter [cm]"         : round(D2*deltat*100,2),
                "Rotor outlet diameter [cm]"        : round(D2*100,2),
                "Diffuser outlet diameter [cm]"     : round(D3*100,2),
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
                    "Absolute diffuser outlet"           : round(c3,1),
                    "Relative tangential rotor inlet"    : round(w1u,1),
                    "Relative rotor inlet"               : round(w1,1),
                    "Relative tangential rotor outlet"   : round(w2u,1),
                    "Relative rotor outlet"              : round(w2,1),
                    "Blade speed rotor inlet"            : round(u1,1),
                    "Blade speed rotor outlet"           : round(u2,1)
                },
                "Mach Numbers" : {
                    "Ma at rotor inlet"                  : round(Ma1,3),
                    "Ma at rotor outlet"                 : round(Ma2R,3),
                    "Ma at difusser inlet"               : round(Ma2S,3),
                    "Ma at difusser outlet"              : round(Ma3,3)
                },
                "Non-dimensional" : {
                    "Specific speed"                     : round(omega_s,3),
                    "Flow coefficient"                   : round(phi,3),
                    "Flow coefficient'"                  : round(phi_prime,3),
                    "Work coefficient"                   : round(Psi,3),
                    "Isentropic work coefficient"        : round(isPsi,3),
                    "Rotor meridional velocity ratio"    : round(xi,3),
                    "Degree of reaction"                 : round(R,3)
                },
                "Flow angles [\u00b0]" : {
                    "Absolute rotor inlet"               : round(self.alpha1,2),
                    "Absolute rotor outlet"              : round(alpha2,2),
                    "Relative rotor inlet"               : round(beta1,2),
                    "Relative rotor outlet"              : round(beta2,2)
                }
            }

            print(Format.underline + 'Stage '+str(i+1) + Format.end)
            print("Delta enthalpy [kJ/kg]:" + str(round(stageDH/1000,2)))
            print("Compression ratio:" + str(round(PR,2)))
            print("Power [kW]:" + str(round(stagePower/1000,2)))
            print("\033[1mThermodynamic Parameters:\033[0m")
            pprint.pprint(self.thermodynamicDict['Stage '+str(i+1)], sort_dicts=False)
            print("\033[1mGeometric Parameters:\033[0m")
            pprint.pprint(self.geometricDict['Stage '+str(i+1)], sort_dicts=False)
            print("\033[1mKinematic Parameters:\033[0m")
            pprint.pprint(self.kinematicDict['Stage '+str(i+1)], sort_dicts=False)

            if warning == True:
                print("\033[1mWarnings:\033[0m")
                if (omega_s < 0.4) or (omega_s > 1.0):
                    print("Non-optimal specific speed (0.4-1.0): ", round(omega_s,3))
                if (phi < 0.2) or (phi > 0.3):
                    print("Non-optimal flow coefficient (0.2-0.3): ", round(phi,3))
                if (Psi < 0.5) or (Psi > 0.6):
                    print("Non-optimal work coefficient (0.5-0.6): ", round(Psi,3))
                if (R < 0.6) or (R > 0.7):
                    print("Non-optimal degree of reaction (0.6-0.7): ", round(R,3))
                if (deltat < 0.5) or (deltat > 0.75):
                    print("Non-optimal rotor tip diameter ratio (0.5-0.75): ", round(deltat,3))
                if ((self.deltah/deltat) < 0.2) or ((self.deltah/deltat) > 0.7):
                    print("Non-optimal rotor hub/tip ratio (0.2-0.7): ", round(self.deltah/deltat,3))
                if (alpha2 < 60) or (alpha2 > 70):
                    print("Non-optimal rotor outlet absolute flow angle (60\u00b0-70\u00b0): ", round(alpha2,2))
                if (beta2 < 0) or (beta2 > 60):
                    print("Non-optimal rotor outlet relative flow angle (0\u00b0-60\u00b0): ", round(beta2,2))
                if ((w1/w2) < 0.7) or ((w1/w2) > 0.95):
                    print("Non-optimal difussion factor (DF=w1/w2, 0.70-0.95): ", round(w1/w2,2))



    def compDictionary(self):

        self.outlet     = streamProps.ThermoState("PH", self.fluid, self.stages[-1].P, self.stages[-1].H)

        self.DH         = self.outlet.H - self.inlet.H
        self.Power      = self.DH * self.massFlow
        self.compEta_is = self.isen.Power / self.Power
        self.compEta_p  = self.polytropicEfficiency()

        self.compressorDict = {
            "Fluid"                   : self.fluid[9:],
            "Mass flow [kg/s]"        : self.massFlow,
            "Compressor type"         : "Radial compressor",
            "Rotational speed [RPM]"  : self.n,
            "Number of stages"        : self.z,
            "Max. tip velocity [m/s]" : self.u2_max,
            "Inlet pressure [bar]"    : round(self.inlet.P/1e5,2),
            "Inlet temperature [C]"   : round(self.inlet.T-273.15,2),
            "Inlet volumetric [m3/h]" : round(self.inletV*3600,2),
            "Delta enthalpy [kJ/kg]"  : round(self.DH/1000,2),
            "Compression ratio"       : round(self.C_ratio,2),
            "Compressor power [kW]"   : round(self.Power/1000,2),
            "Comp. Isentropic eff. [%]" : round(self.compEta_is*100,2),
            "Comp. Polytropic eff. [%]" : round(self.compEta_p*100,2)
        }

        print(Format.underline + "RADIAL COMPRESSOR FINAL SPECIFICATION" + Format.end)
        pprint.pprint(self.compressorDict, sort_dicts=False)

    def MachNumbers(self, stageInletA, rotorOutletStateA, stageOutletA, w1, w2, c2, c3):

        Ma1       = w1 / stageInletA
        Ma2R      = w2 / rotorOutletStateA
        Ma2S      = c2 / rotorOutletStateA
        Ma3       = c3 / stageOutletA

        return Ma1, Ma2R, Ma2S, Ma3

    def rotorOutlet(self, stageInletH, stageInletS, w1, u1, w2, u2, eta_R):

        rothalpyInlet        = stageInletH + (w1**2)/2 - (u1**2)/2
        rotorOutletH         = rothalpyInlet - (w2**2)/2 + (u2**2)/2
        rotorOutletH_is      = stageInletH + eta_R*(rotorOutletH-stageInletH)
        rotorOutlet_is       = streamProps.ThermoState("HS", self.fluid, rotorOutletH_is, stageInletS)
        rotorOutlet          = streamProps.ThermoState("PH", self.fluid, rotorOutlet_is.P, rotorOutletH)

        return rotorOutlet

    def velocityKinematics(self, u2, phi, xi, alpha2, beta1, beta2, deltat):

        tanAlpha1 = math.tan(self.alpha1*(math.pi/180))
        tanAlpha2 = math.tan(alpha2*(math.pi/180))
        tanBeta1  = math.tan(beta1*(math.pi/180))
        tanBeta2  = math.tan(beta2*(math.pi/180))

        c1m       = u2*phi
        c1u       = c1m*tanAlpha1
        c2m       = u2*xi*phi
        c2u       = c2m*tanAlpha2
        w1u       = u2*phi*tanBeta1
        w2u       = u2*xi*phi*tanBeta2
        c1        = u2*phi*math.sqrt(1+(tanAlpha1**2))
        c3        = c1 #Assumption
        c2        = u2*xi*phi*math.sqrt(1+(tanAlpha2**2))
        w1        = u2*phi*math.sqrt(1+(tanBeta1**2))
        w2        = u2*xi*phi*math.sqrt(1+(tanBeta2**2))
        u1        = u2*deltat

        return c1m, c1u, c2m, c2u, w1u, w2u, c1, c2, c3, w1, w2, u1

    def designParameters(self, n, inletV, DHperStage, inletH, inletP, outletP):

        omega_s   = self.specificSpeed(n, inletV, DHperStage)
        isPsi     = self.isentropicWorkCoeff(omega_s)
        deltat    = self.rotorTipDiameterRatio(omega_s,isPsi)
        u2        = math.sqrt(DHperStage/isPsi)
        D2        = 60*u2/(math.pi*n)
        b1        = self.bladeHeightRotorInlet(D2, deltat)
        c1m       = self.rotorInletMeridionalVelocity(inletV, D2, deltat, b1)
        phi       = c1m/u2
        phi_prime = inletV/(u2*D2)
        D3        = D2*(1.55 + ((deltat**2)-(self.deltah**2))*phi)
        alpha2    = self.absoluteFlowAngles2(omega_s, isPsi)
        eta_p     = self.polytropicEfficiencyEst(omega_s)
        eta_is    = self.polytropic2isentropic(eta_p, inletH, inletP, outletP)
        eta_R     = eta_is #Assumption
        Psi       = isPsi/eta_is
        xi        = self.rotorMeridionalVelocityRatio(Psi, phi, deltat, alpha2)
        R         = self.degreeReaction(Psi, phi, xi, deltat)
        beta1, beta2 = self.relativeFlowAngles(deltat, phi, alpha2, Psi, xi)

        return omega_s, isPsi, deltat, u2, D2, b1, c1m, phi, phi_prime, D3, \
                alpha2, eta_p, eta_R, eta_is, Psi, xi, R, beta1, beta2


    def specificSpeed(self, n, inletV, DHperStage):

        omega_s   = (2*math.pi*n / 60) * (inletV**0.5 / (DHperStage)**0.75)

        return omega_s

    def isentropicWorkCoeff(self, omega_s):

        t1        = 4 * (-0.3 + math.log10(omega_s))
        A         = 1/(1 + math.exp(-t1))
        t2        = 5 * (1.0 + math.log10(omega_s))
        B         = math.exp(-t2)
        isPsi     = 0.55*(1-A) + 0.02*(A) + (0.55-0.45)*B

        return isPsi

    def rotorTipDiameterRatio(self, omega_s, isPsi):

        deltat    = 0.5 + (1.5/math.pi)*(omega_s**2)*(isPsi**1.5)

        return deltat

    def bladeHeightRotorInlet(self, D2, deltat):

        b1        = (D2/2) * (deltat - self.deltah)

        return b1

    def rotorInletMeridionalVelocity(self, inletV, D2, deltat, b1):

        c1m       = (2 * inletV) / (math.pi * D2 * (deltat + self.deltah) * b1)

        return c1m

    def absoluteFlowAngles2(self, omega_s, isPsi):

        A         = ((omega_s**2)*(isPsi**1.5)) / math.pi
        alpha2    = 72 - 0.5*math.log(A) - 585*(A**2)

        return alpha2

    def polytropicEfficiencyEst(self, omega_s):

        w         = math.log10(omega_s / 2.9809)
        eta_p     = 10**(-0.097358 - 0.0800538*w + 0.151771*(w**2) + 0.340467*(w**3))

        return eta_p

    def polytropic2isentropic(self, eta_p, Hin, Pin, Pout):

        inlet      = streamProps.ThermoState("PH", self.fluid, Pin, Hin)
        isOutlet   = streamProps.ThermoState("PS", self.fluid, Pout, inlet.S)

        eta_is     = 0.7 #initial guess
        error      = 10  #initial value

        while error > 0.001:

            outlet_H   = inlet.H + (isOutlet.H - inlet.H) / eta_is
            outlet     = streamProps.ThermoState("PH", self.fluid, Pout, outlet_H)

            s1h        = streamProps.ThermoState("PH", self.fluid, Pin, outlet.H)

            guessEta_p = (s1h.S - outlet.S) / (s1h.S - inlet.S)
            error      = abs(eta_p - guessEta_p)

            eta_is    += 0.001

        return eta_is

    def polytropicEfficiency(self):

        isobarH    = streamProps.ThermoState("PH", self.fluid, self.inlet.P, self.outlet.H)

        compEta_p  = (isobarH.S-self.outlet.S) / (isobarH.S-self.inlet.S)

        return compEta_p

    def rotorMeridionalVelocityRatio(self, Psi, phi, deltat, alpha2):

        xi        = (Psi + phi*deltat*math.tan(self.alpha1*(math.pi/180))) / (phi*math.tan(alpha2*(math.pi/180)))

        return xi

    def degreeReaction(self, Psi, phi, xi, deltat):

        A         = (1-(xi**2)) + (math.tan(self.alpha1*(math.pi/180))**2)*(1-(deltat**2))
        R         = 1 - Psi/2 + (phi**2)*A/(2*Psi) - phi*deltat*math.tan(self.alpha1*(math.pi/180))

        return R

    def relativeFlowAngles(self, deltat, phi, alpha2, Psi, xi):

        beta1     = math.atan((deltat/phi) - math.tan(self.alpha1*(math.pi/180))) * 180/math.pi
        beta2     = math.atan((1-Psi)/(phi*xi) - (deltat/xi)*math.tan(self.alpha1*(math.pi/180))) * 180/math.pi

        return beta1, beta2

    class Isentropic:

        def __init__(self, pair, inlet, outlet_value, massFlow):

            #outlet isentropic thermo props
            if   pair == "PS":
                self.outlet = streamProps.ThermoState("PS", inlet.fluid, outlet_value, inlet.S)

            elif pair == "HS":
                self.outlet = streamProps.ThermoState("HS", inlet.fluid, outlet_value, inlet.S)

            #isentropic power
            self.DH     = self.outlet.H - inlet.H    #J/kg
            self.Power  = self.DH * massFlow         #W

            #isentropic work coefficient
            self.Psi    = 0.45

        def printPower(self):
            print("\u0394h_is [kJ/kg]  = ", round(self.DH/1e3,3))
            print("Power_is [kW]  =", round(self.Power/1e3,3))
