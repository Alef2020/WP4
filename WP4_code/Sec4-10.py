#SECTION 4.10


import math
import numpy as np
from scipy.optimize import minimize

#Equations for the compliances
class AttachedPart:
    def __init__(self,E, t, Dfo, Dfi):
        self.E = E
        self.t = t
        self.Dfo = Dfo
        self.Dfi = Dfi

    
    def compliance(self):
        return 4 * self.t / (self.E * math.pi * (self.Dfo**2 - self.Dfi**2))
    
class Fastener:
    def __init__(self, E, d_nom, d_thread_minor, k_h=0.5, k_eng=0.4, k_n=0.4): #add justification in WP4 for why we choose these default values, with regards to the nut shape in reference B
        self.E = E
        self.d_nom = d_nom
        self.d_thread_minor = d_thread_minor
        self.k_h = k_h
        self.k_eng = k_eng
        self.k_n = k_n

    def compliance(self):
        import math

        d = self.d_nom
        A_shank  = math.pi * (d / 2)**2
        A_thread = math.pi * (self.d_thread_minor / 2)**2

        Lh_sub        = self.k_h   * d
        L_engaged     = self.k_eng * d
        Ln_sub        = self.k_n   * d

        # simple split; all shank, no thread. "A fastener is selected such that the unthreaded shank spans the entire grip length, so threads are not engaged in the joint region. Thus, the engaged region is modelled entirely as shank."
        L_shank   = L_engaged      
        L_thread  = 0.0            

        delta = 0.0
        delta += Lh_sub   / (self.E * A_shank)
        delta += L_shank  / (self.E * A_shank)
        delta += L_thread / (self.E * A_thread)
        delta += Ln_sub   / (self.E * A_shank)

        return delta

    

#force ratio equation
def force_ratio(delta_a, delta_b):
    return delta_a / (delta_a + delta_b)








#Fixed Variables:
E_backplate = 70e9      # Pa
E_wall      = 70e9      # Pa
t_backplate = 0.004     # m  (t2)
t_wall      = 0.003     # m  (t3)
Dfo         = 0.012     # m  (outer diameter under head)
Dfi         = 0.006     # m  (inner diameter under head / hole diameter)


#Fastener Material
E_fastener = 210e9      # Pa
thread_minor_factor = 0.85


#Optimization Section


#Function for D_NOM
def phi_diamater(d_nom):


    delta_a = AttachedPart(E_backplate, t_backplate, Dfo, Dfi).compliance() + AttachedPart(E_wall, t_wall, Dfo, Dfi).compliance()
    delta_b = Fastener(E_fastener, d_nom, d_nom * thread_minor_factor).compliance()

    phi = force_ratio(delta_a, delta_b)
    return phi, delta_a, delta_b


#Objective for minimization

    #choose phi target to be 0.5

def objective(d_nom):
    phi_target = 0.5
    phi, _, _ = phi_diamater(d_nom)
    return (phi - phi_target)**2 # or abs(phi - phi_target)


#Optimization run
def optimize_d_nom():
    bounds = [(0.001, 0.03)]  # reasonable bounds for d_nom
    x0 = [0.01]

    result = minimize(objective, x0, method = 'SLSQP', bounds=bounds, options = {'maxiter':10000})

    d_opt = result.x[0]
    print ("Optimal d_nom (m):", d_opt)

    phi_opt, delta_a_opt, delta_b_opt = phi_diamater(d_opt)
    print ("Optimal phi:", phi_opt)
    print ("Delta_a (m/N):", delta_a_opt)
    print ("Delta_b (m/N):", delta_b_opt)

if __name__ == "__main__":
    optimize_d_nom()

