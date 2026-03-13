import numpy as np 

G = 6.673848000000006e-08
C = 29979245800.0


class EOS4:

    def __init__(self,log_p1,Gamma1,Gamma2,Gamma3):


        self.rho1 = pow(10, 14.7) 
        self.rho2 = pow(10, 15.0) 

        self.p1 = pow(10.0, log_p1) 
        self.Gamma1 = Gamma1
        self.Gamma2 = Gamma2
        self.Gamma3 = Gamma3

        self.rho13 = 1.0e13  
        self.p13 = 1.5689 * 10**31 
        self.Gamma0 = 1.35692395
        self.K0 = self.p13 / self.rho13**self.Gamma0

        self.K1 = self.p1 / pow(self.rho1, self.Gamma1)
        self.K2 = self.K1 * pow(self.rho1, self.Gamma1 - self.Gamma2)
        self.K3 = self.K2 * pow(self.rho2, self.Gamma2 - self.Gamma3)

        self.epsL = 0.0
        self.alphaL = 0.0

        self.rho0 = pow(self.K0 / self.K1, 1.0 / (self.Gamma1 - self.Gamma0))

        self.eps0 = (1.0 + self.alphaL) * self.rho0 + self.K0 / (self.Gamma0 - 1.0) * pow(self.rho0, self.Gamma0)/C**2

        self.alpha1 = self.eps0 / self.rho0 - 1 - self.K1 / (self.Gamma1 - 1) * pow(self.rho0, self.Gamma1 - 1)/C**2
        self.eps1 = (1 + self.alpha1) * self.rho1 + self.K1 / (self.Gamma1 - 1) * pow(self.rho1, self.Gamma1)/C**2
        self.alpha2 = self.eps1 / self.rho1 - 1 - self.K2 / (self.Gamma2 - 1) * pow(self.rho1, self.Gamma2 - 1)/C**2
        self.eps2 = (1 + self.alpha2) * self.rho2 + self.K2 / (self.Gamma2 - 1) * pow(self.rho2, self.Gamma2)/C**2
        self.alpha3 = self.eps2 / self.rho2 - 1 - self.K3 / (self.Gamma3 - 1) * pow(self.rho2, self.Gamma3 - 1)/C**2

        self.p0 = self.K0 * pow(self.rho0, self.Gamma0)
        self.p2 = self.K2 * pow(self.rho2, self.Gamma2)

    def rho2p(self, rho):
        if rho < self.rho0:
            return self.K0 * pow(rho, self.Gamma0)
        elif self.rho0 <= rho < self.rho1:
            return self.K1 * pow(rho, self.Gamma1)
        elif self.rho1 <= rho < self.rho2:
            return self.K2 * pow(rho, self.Gamma2)
        else:
            return self.K3 * pow(rho, self.Gamma3)

    def p2rho(self, p):
        if p < self.p0:
            return pow(p / self.K0, 1.0 / self.Gamma0)
        elif self.p0 <= p < self.p1:
            return pow(p / self.K1, 1.0 / self.Gamma1)
        elif self.p1 <= p < self.p2:
            return pow(p / self.K2, 1.0 / self.Gamma2)
        else:
            return pow(p / self.K3, 1.0 / self.Gamma3)

    def rho2e(self, rho):
        if rho < self.rho0:
            return (1.0 + self.alphaL) * rho + self.K0 / (self.Gamma0 - 1.0) * pow(rho, self.Gamma0)/C**2
        elif self.rho0 <= rho < self.rho1:
            return (1.0 + self.alpha1) * rho + self.K1 / (self.Gamma1 - 1.0) * pow(rho, self.Gamma1)/C**2
        elif self.rho1 <= rho < self.rho2:
            return (1.0 + self.alpha2) * rho + self.K2 / (self.Gamma2 - 1.0) * pow(rho, self.Gamma2)/C**2
        else:
            return (1.0 + self.alpha3) * rho + self.K3 / (self.Gamma3 - 1.0) * pow(rho, self.Gamma3)/C**2

    def p2e(self, p):

        if p < self.p0:
            rho = pow(p / self.K0, 1.0 / self.Gamma0)
            return (1.0 + self.alphaL) * rho + p / (self.Gamma0 - 1)/C**2
        elif self.p0 <= p < self.p1:
            rho = pow(p / self.K1, 1.0 / self.Gamma1)
            return (1.0 + self.alpha1) * rho + p / (self.Gamma1 - 1)/C**2
        elif self.p1 <= p < self.p2:
            rho = pow(p / self.K2, 1.0 / self.Gamma2)
            return (1.0 + self.alpha2) * rho + p / (self.Gamma2 - 1)/C**2
        else:
            rho = pow(p / self.K3, 1.0 / self.Gamma3)
            return (1.0 + self.alpha3) * rho + p / (self.Gamma3 - 1)/C**2
        
    def rho2gamma(self,rho):

        if rho < self.rho0:

            return self.Gamma0
        
        elif self.rho0 <= rho < self.rho1:

            return self.Gamma1

        elif self.rho1 <= rho < self.rho2:

            return self.Gamma2
        
        else:
            return self.Gamma3

class EOS7:

    def __init__(self,log_p1,Gamma1,Gamma2,Gamma3):


        self.rho1 = pow(10, 14.7) 
        self.rho2 = pow(10, 15.0) 

        self.p1 = pow(10.0, log_p1) 
        self.Gamma1 = Gamma1
        self.Gamma2 = Gamma2
        self.Gamma3 = Gamma3

        self.K1 = self.p1 / pow(self.rho1, self.Gamma1)
        self.K2 = self.K1 * pow(self.rho1, self.Gamma1 - self.Gamma2)
        self.K3 = self.K2 * pow(self.rho2, self.Gamma2 - self.Gamma3)

        self.rhoL_1 = 2.62789e12 
        self.rhoL_2 = 3.78358e11 
        self.rhoL_3 = 2.44034e7 
        self.rhoL_4 = 0.0

        self.GammaL_1 = 1.35692
        self.GammaL_2 = 0.62223
        self.GammaL_3 = 1.28733
        self.GammaL_4 = 1.58425

        self.KL_1 = 3.99874e-8 * C**2
        self.KL_2 = 5.32697e+1 * C**2
        self.KL_3 = 1.06186e-6 * C**2
        self.KL_4 = 6.80110e-9 * C**2

        self.epsL_4 = 0.0
        self.alphaL_4 = 0.0
        self.epsL_3 = (1 + self.alphaL_4) * self.rhoL_3 + self.KL_4 / (self.GammaL_4 - 1) * pow(self.rhoL_3, self.GammaL_4) / C**2
        self.alphaL_3 = self.epsL_3 / self.rhoL_3 - 1 - self.KL_3 / (self.GammaL_3 - 1) * pow(self.rhoL_3, self.GammaL_3 - 1) / C**2
        self.epsL_2 = (1 + self.alphaL_3) * self.rhoL_2 + self.KL_3 / (self.GammaL_3 - 1) * pow(self.rhoL_2, self.GammaL_3) / C**2
        self.alphaL_2 = self.epsL_2 / self.rhoL_2 - 1 - self.KL_2 / (self.GammaL_2 - 1) * pow(self.rhoL_2, self.GammaL_2 - 1) / C**2
        self.epsL_1 = (1 + self.alphaL_2) * self.rhoL_1 + self.KL_2 / (self.GammaL_2 - 1) * pow(self.rhoL_1, self.GammaL_2) / C**2
        self.alphaL_1 = self.epsL_1 / self.rhoL_1 - 1 - self.KL_1 / (self.GammaL_1 - 1) * pow(self.rhoL_1, self.GammaL_1 - 1) / C**2

        self.rho0 = pow(self.KL_1 / self.K1, 1.0 / (self.Gamma1 - self.GammaL_1))
        self.eps0 = (1.0 + self.alphaL_1) * self.rho0 + self.KL_1 / (self.GammaL_1 - 1.0) * pow(self.rho0, self.GammaL_1) / C**2
        
        self.alpha1 = self.eps0/self.rho0 - 1 - self.K1/(self.Gamma1 - 1)*pow(self.rho0, self.Gamma1 -1)/C**2
        self.eps1 = (1+self.alpha1)*self.rho1 + self.K1/(self.Gamma1 - 1)*pow(self.rho1, self.Gamma1)/C**2
        self.alpha2 = self.eps1/self.rho1 - 1 - self.K2/(self.Gamma2 - 1)*pow(self.rho1, self.Gamma2 -1)/C**2
        self.eps2 = (1+self.alpha2)*self.rho2 + self.K2/(self.Gamma2 - 1)*pow(self.rho2, self.Gamma2)/C**2
        self.alpha3 = self.eps2/self.rho2 - 1 - self.K3/(self.Gamma3 - 1)*pow(self.rho2, self.Gamma3 -1)/C**2

        self.pL_3 = self.KL_3*pow(self.rhoL_3,self.GammaL_3)
        self.pL_2 = self.KL_2*pow(self.rhoL_2,self.GammaL_2)
        self.pL_1 = self.KL_1*pow(self.rhoL_1,self.GammaL_1)
        self.p0 = self.KL_1*pow(self.rho0,self.GammaL_1)
        self.p2 = self.K2*pow(self.rho2,self.Gamma2)


    def rho2p(self, rho):

        if rho < self.rhoL_3:
            return self.KL_4 * pow(rho, self.GammaL_4)
        
        elif self.rhoL_3 <= rho < self.rhoL_2:
            return self.KL_3 * pow(rho, self.GammaL_3)
        
        elif self.rhoL_2 <= rho < self.rhoL_1:
            return self.KL_2 * pow(rho, self.GammaL_2)
        
        elif self.rhoL_1 <= rho < self.rho0:
            return self.KL_1 * pow(rho, self.GammaL_1)
        
        elif self.rho0 <= rho < self.rho1:
            return self.K1 * pow(rho, self.Gamma1)
        
        elif self.rho1 <= rho < self.rho2:
            return self.K2 * pow(rho, self.Gamma2)
        
        else:
            return self.K3 * pow(rho, self.Gamma3)


    def rho2e(self, rho):

        if rho < self.rhoL_3:
            return (1.0 + self.alphaL_4) * rho + self.KL_4 / (self.GammaL_4 - 1.0) * pow(rho, self.GammaL_4)/C**2
        
        elif self.rhoL_3 <= rho < self.rhoL_2:
            return (1.0 + self.alphaL_3) * rho + self.KL_3 / (self.GammaL_3 - 1.0) * pow(rho, self.GammaL_3)/C**2
        
        elif self.rhoL_2 <= rho < self.rhoL_1:
            return (1.0 + self.alphaL_2) * rho + self.KL_2 / (self.GammaL_2 - 1.0) * pow(rho, self.GammaL_2)/C**2
        
        elif self.rhoL_1 <= rho < self.rho0:
            return (1.0 + self.alphaL_1) * rho + self.KL_1 / (self.GammaL_1 - 1.0) * pow(rho, self.GammaL_1)/C**2
        
        elif self.rho0 <= rho < self.rho1:
            return (1.0 + self.alpha1) * rho + self.K1 / (self.Gamma1 - 1.0) * pow(rho, self.Gamma1)/C**2
        
        elif self.rho1 <= rho < self.rho2:
            return (1.0 + self.alpha2) * rho + self.K2 / (self.Gamma2 - 1.0) * pow(rho, self.Gamma2)/C**2
        
        else:
            return (1.0 + self.alpha3) * rho + self.K3 / (self.Gamma3 - 1.0) * pow(rho, self.Gamma3)/C**2

        
    def rho2gamma(self,rho):

        if rho < self.rhoL_3:
            return self.GammaL_4
        
        elif self.rhoL_3 <= rho < self.rhoL_2:
            return self.GammaL_3
        
        elif self.rhoL_2 <= rho < self.rhoL_1:
            return self.GammaL_2
        
        elif self.rhoL_1 <= rho < self.rho0:
            return self.GammaL_1
        
        elif self.rho0 <= rho < self.rho1:
            return self.Gamma1
        
        elif self.rho1 <= rho < self.rho2:
            return self.Gamma2
        
        else:
            return self.Gamma3




    