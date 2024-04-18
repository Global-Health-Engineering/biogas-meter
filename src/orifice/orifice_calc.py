import os
import sys
import itertools
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import optimize
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

if sys.platform == 'linux':
    os.environ['RPPREFIX'] = r'/etc/REFPROP'
elif sys.platform == 'darwin':
    os.environ['RPPREFIX'] = r'/Users/jtkaczuk/codes/REFPROP'
else:
    os.environ['RPPREFIX'] = r'c:/Program Files (x86)/REFPROP'


class FluidProps(object):
    fluid = str or tuple

    def __init__(self, fluid):
        self.RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        self.RP.SETPATHdll(os.environ['RPPREFIX'])
        self.MASS_BASE_SI = self.RP.GETENUMdll(0, "MASS BASE SI").iEnum
        self.MOLAR_BASE_SI = self.RP.GETENUMdll(0, "MOLAR BASE SI").iEnum
        self.fluid = self.set_fluid(fluid)
        self.z = {1.0}

    def set_fluid(self, fluid):
        return fluid if type(fluid) == str else '{};{}'.format(*fluid)

    def set_composition(self, z):
        self.z = z

    def set_composition_from_1st_fraction(self, x):
        self.z = [round(x, 8), round(1-x, 8)]

    def get_compressibility_factor(self, p, T):
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'Z',
                                 self.MASS_BASE_SI, 0, 0,
                                 p, T, self.z)
        return res.Output[0]

    def get_p(self, rho, T):
        """
        return pressure [Pa]
        as a function of mass density [kg/m3]
        and temperature [K]
        """
        res = self.RP.REFPROPdll(self.fluid, 'DT', 'P',
                                 self.MASS_BASE_SI, 0, 0,
                                 rho, T, self.z)
        return res.Output[0]

    def get_rhomass(self, p, T):
        """
        return mass density [kg/m3]
        as a function of pressure [Pa]
        and temperature [K]
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'D',
                                 self.MASS_BASE_SI, 0, 0,
                                 p, T, self.z)
        return res.Output[0]

    def get_rhomolar(self, p, T):
        """
        return molar density [mol/m3]
        as a function of pressure [Pa]
        and temperature [K]
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'D',
                                 self.MOLAR_BASE_SI, 0, 0,
                                 p, T, self.z)
        return res.Output[0]

    def get_cp(self, p, T):
        """
        return isobaric specific heat [J/kg/K]
        as a function of pressure [Pa]
        and temperature [K]
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'CP',
                                 self.MASS_BASE_SI, 0, 0,
                                 p, T, self.z)
        return res.Output[0]

    def get_specific_heat_ratio(self, p, T):
        """
        return the specific heat ratio [-]
        as a function of pressure [Pa]
        and temperature [K]
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'CP/CV',
                                 self.MASS_BASE_SI, 0, 0,
                                 p, T, self.z)
        return res.Output[0]

    def get_hmass(self, p, T):
        """return enthalpy [J/kg]"""
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'H',
                                 self.MASS_BASE_SI, 0, 0,
                                 p, T, self.z)
        return res.Output[0]

    def get_T(self, p, h):
        """
        return temperature [K]
        as a function of pressure [Pa]
        and enthalpy in [J/kg]
        """
        res = self.RP.REFPROPdll(self.fluid, 'PH', 'T',
                                 self.MASS_BASE_SI, 0, 0,
                                 p, h, self.z)
        return res.Output[0]

    def get_critical_T(self):
        """ return critical temperature [K] """
        res = self.RP.REFPROPdll(self.fluid, 'PQ', 'TC',
                                 self.MASS_BASE_SI, 0, 0,
                                 1, 1, self.z)
        return res.Output[0]

    def get_critical_p(self):
        """ return critical pressure [Pa] """
        res = self.RP.REFPROPdll(self.fluid, 'TQ', 'PC',
                                 self.MASS_BASE_SI, 0, 0,
                                 1, 1, self.z)
        return res.Output[0]

    def get_speed_of_sound(self, p, T):
        """
        return speed of sound [m/s]
        as a function of pressure [Pa]
        and temperature [K]
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'W',
                                 self.MASS_BASE_SI, 0, 0,
                                 p, T, self.z)
        return res.Output[0]

    def get_JT_coefficient(self, p, T):
        """
        return the Joule-Thomson coefficient [K/Pa]
        as a function of pressure [Pa]
        and temperature [K]
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'JT',
                                 self.MASS_BASE_SI, 0, 0,
                                 p, T, self.z)
        return res.Output[0]

    def get_viscosity(self, p, T):
        """
        return dynamic viscosity [Pa-s]
        as a function of pressure [Pa]
        and temperature [K]
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'VIS',
                                 self.MASS_BASE_SI, 0, 0,
                                 p, T, self.z)
        return res.Output[0]

    def get_thermal_conductivity(self, p, T):
        """
        return thermal conductivity [W/m/K]
        as a function of pressure [Pa]
        and temperature [K]
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'TCX',
                                 self.MASS_BASE_SI, 0, 0,
                                 p, T, self.z)
        return res.Output[0]



class Orifice(object):
    fluid = str or tuple
    D_out = float
    D_in = float
    
    def __init__(self, fluid, D_out, D_in):
        self.fluid = FluidProps(fluid)
        self.D_out = D_out
        self.D_in = D_in
        self.beta = D_in / D_out

    def set_composition(self, x):
        self.fluid.set_composition_from_1st_fraction(x)

    def set_D_in(self, D_in):
        self.D_in = D_in

    def set_D_out(self, D_out):
        self.D_out = D_out

    def get_beta(self):
        return self.beta

    def get_u(self, m_flow, p, T, d):
        """ return fluid velocity in pipe [m/s] """
        return m_flow / self.fluid.get_rhomass(p=p, T=T) / (np.pi * d**2 / 4)

    def get_Re(self, m_flow, p, T, d):
        """ return Reynolds number """
        rho = self.fluid.get_rhomass(p=p, T=T)
        viscosity = self.fluid.get_viscosity(p=p, T=T)
        u = self.get_u(m_flow=m_flow, p=p, T=T, d=d)
        return rho * u * d / viscosity

    def get_discharge_coefficient(self, m_flow, p, T):
        Re_Dout = self.get_Re(m_flow=m_flow, p=p, T=T, d=self.D_out)
        A = (19000 * self.beta / Re_Dout)**0.8
        L_1, L_2prime = 1, 0.47 # D and D/2 tappings
        M_2prime = 2 * L_2prime / (1 - self.beta)
        C = [+ 0.5961,
             + 0.0261 * self.beta**2,
             - 0.216 * self.beta**8,
             + 0.000521 * (1e6 * self.beta / Re_Dout)**0.7,
             + (0.0188 + 0.0063 * A) * self.beta**3.5 * (1e6 / Re_Dout)**0.3,
             + ((0.043 + 0.080 * np.exp(-10 * L_1) - 0.123 * np.exp(-7 * L_1))
                    * (1 - 0.11 * A) * self.beta**4 / (1 - self.beta**4)),
             - 0.031 * (M_2prime - 0.8 * M_2prime**1.1) * self.beta**1.3,
             + 0.011 * (0.75 - self.beta) * (2.8 - self.D_out / 0.0254)]
        if self.D_out < 71.2e-3:
            return sum(C)
        else:
            return sum(C[:-1])

    def get_expansibility_factor(self, p1, p2, T1, T2):
        if p2/p1 > 0.75:
            cp_per_cv = np.mean([self.fluid.get_specific_heat_ratio(p=p1, T=T1),
                                 self.fluid.get_specific_heat_ratio(p=p2, T=T2)])
            return (1 - (0.351 + 0.256 * self.beta**4 + 0.93 * self.beta**8)
                    * (1 - (p2/p1)**(1/cp_per_cv)))
        else:
            return 1

    def get_mass_flow_rate(self, p1, p2, T1, T2):
        """ ISO 5167 calculations of the mass flow rate [kg/s] """
        epsilon = self.get_expansibility_factor(p1=p1, p2=p2, T1=T1, T2=T2)
        rhomass1 = self.fluid.get_rhomass(p=p1, T=T1)
        def q(x):
            return (self.get_discharge_coefficient(m_flow=x, p=(p1+p2)/2, T=(T1+T2)/2)
                    / np.sqrt(1 - self.beta**4) * epsilon * np.pi / 4
                    * self.D_in**2 * np.sqrt(2 * rhomass1 * (p1 - p2)) - x)
        return optimize.fsolve(func=q, x0=1e-3, epsfcn=1e-10)[0]


def set_mpl(fontsize=16):
    mpl.rcParams.update({"font.family": "Times New Roman",
                         "mathtext.fontset": "dejavuserif",
                         "font.size": fontsize,
                         # "axes.linewidth": 1.,
                         "axes.labelsize": fontsize-1,
                         "axes.titlesize" : fontsize-1,
                         "legend.fontsize": fontsize-1,
                         "xtick.top": True,
                         "xtick.bottom": True,
                         "ytick.left": True,
                         "ytick.right": True,
                         "xtick.direction": "in",
                         "ytick.direction": "in",
                         "xtick.major.pad": 5,
                         "ytick.major.pad": 5})
                         # "xtick.major.width": 1,
                         # "xtick.minor.width": 1,
                         # "ytick.major.width": 1,
                         # "ytick.minor.width": 1})


def test_HeNe_EOS_implementation():
    os.environ['RPPREFIX'] = r'/Users/jtkaczuk/codes/REFPROP'
    EOS = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
    EOS.SETPATHdll(os.environ['RPPREFIX'])
    MOLAR_BASE_SI = EOS.GETENUMdll(0, 'MOLAR BASE SI').iEnum
    res = EOS.REFPROPdll('HELIUM;NEON', 'DT', 'P', MOLAR_BASE_SI, 0, 0, 10e3, 200, [0.5, 0.5])
    print(f'implementation results (REFPROP): {res.Output[0]} {res.hUnits}')
    print('correct results (CoolProp)      : 18430775.292600896 Pa')


def orifice_test():
    D_out = 75 # mm
    D_in = 37.5 # mm
    T = 300
    p1 = 1.02e5
    p2 = 1e5
    fluid = 'methane'
    # fluid = 'carbon dioxide'
    orifice = Orifice(fluid, D_out*1e-3, D_in*1e-3)
    m_flow = orifice.get_mass_flow_rate(p1, p2, T, T)
    print("{:.2f} g/s".format(m_flow*1e3))


def main():
    pass


if __name__ == "__main__":
    set_mpl()
    # test_HeNe_EOS_implementation()
    orifice_test()
    # main()

