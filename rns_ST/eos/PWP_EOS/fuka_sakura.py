import numpy as np 
G = 6.673848000000006e-08 
Msun = 1.988399999999998e+33 
C = 29979245800.0
Length = G*Msun/C**2
Time = Length/C
Density = Msun/Length**3

class fuka_sakura:
    def __init__(self, log_p1, Gamma1, Gamma2, Gamma3):
        self.rho1 = 10**14.7 /Density
        self.rho2 = 10**15.0 /Density

        self.p1 = 10**log_p1/Density/C**2
        self.Gamma1 = Gamma1
        self.Gamma2 = Gamma2
        self.Gamma3 = Gamma3

        self.rho13 = 1.0e13  /Density
        self.p13 = 1.5689*10**31/Density/C**2
        self.Gamma0 = 1.35692395
        self.K0 = self.p13 / self.rho13**self.Gamma0

        self.K1 = self.p1 / self.rho1**self.Gamma1
        self.K2 = self.K1 * self.rho1**(self.Gamma1 - self.Gamma2)
        self.K3 = self.K2 * self.rho2**(self.Gamma2 - self.Gamma3)

        self.rho0 = pow(self.K0 / self.K1, 1.0 / (self.Gamma1 - self.Gamma0))

    def fuka(self, eos):
        filename = f'{eos}.polytrope'
        content = f"""# This is a table holding the information for a piecewise polytropic equation of state
# which can be ingested into Margherita to construct the corresponding analytical polytrope
# K0: constant of the zeroth boundary.
# Pressure0: pressure of the zeroth boundary.

num_pieces: 4
rhomin: 1.e-19
rhomax: 1.0
ktab0: {self.K0}
Ptab0: 0

gamma_tab
1.3569199999999999 {self.Gamma1} {self.Gamma2} {self.Gamma3}

rho_tab
0e1 {self.rho0} {self.rho1} {self.rho2}

units: geometrised
"""


        with open(filename, 'w') as file:
            file.write(content)

        print(f"written in {filename}")

    def sakura(self, eos):
        filename = f'{eos}.polytrope'
        content = f"""{self.Gamma0:.16E} {self.Gamma1:.16E} {self.Gamma2:.16E} {self.Gamma3:.16E} {self.rho0:.16E} {self.rho1:.16E} {self.rho2:.16E} {self.K0:.16E}"""
        content = content.replace('E', 'D')
        with open(filename, 'w') as file:
            file.write(content)







    