#SECTION 4.10


import math

class AttachedPart:
    def __init__(self,E, t, Dfo, Dfi):
        self.E = E
        self.t = t
        self.Dfo = Dfo
        self.Dfi = Dfi

    
    def compliance(self):
        return 4 * self.t / (self.E * math.pi * (self.Dfo**2 - self.Dfi**2))
    
class Fastener:
    def __init__(self,E, d_shank, d_thread_minor, Lh_sub, L_engaged_shank, L_engaged_thread, Ln_sub):
        self.E = E
        self.d_shank = d_shank
        self.d_thread_minor = d_thread_minor
        self.Lh_sub = Lh_sub
        self.L_engaged_shank = L_engaged_shank
        self.L_engaged_thread = L_engaged_thread
        self.Ln_sub = Ln_sub


    def compliance(self):
        
        # area
        A_shank = math.pi * (self.d_shank / 2)**2
        A_thread = math.pi * (self.d_thread_minor / 2)**2

        #summation of L/EA

        delta = 0
        delta += self.Lh_sub / (self.E * A_shank)
        delta += self.L_engaged_shank / (self.E * A_shank)
        delta += self.L_engaged_thread / (self.E * A_thread)
        delta += self.Ln_sub / (self.E * A_shank)
        return delta
    
def force_ratio(delta_a, delta_b):
    return delta_b / (delta_a + delta_b)

# def main():
#     # Example values
#     E_attached = 200e9  # Pa
#     t_attached = 0.01   # m
#     Dfo_attached = 0.1  # m
#     Dfi_attached = 0.08 # m

#     E_fastener = 210e9  # Pa
#     d_shank_fastener = 0.02  # m
#     d_thread_minor_fastener = 0.018  # m
#     Lh_sub_fastener = 0.03  # m
#     L_engaged_shank_fastener = 0.04  # m
#     L_engaged_thread_fastener = 0.05  # m
#     Ln_sub_fastener = 0.02  # m

#     attached_part = AttachedPart(E_attached, t_attached, Dfo_attached, Dfi_attached)
#     fastener = Fastener(E_fastener, d_shank_fastener, d_thread_minor_fastener,
#                         Lh_sub_fastener, L_engaged_shank_fastener,
#                         L_engaged_thread_fastener, Ln_sub_fastener)

#     delta_a = attached_part.compliance()
#     delta_b = fastener.compliance()

#     ratio = force_ratio(delta_a, delta_b)

#     print(f"Attached Part Compliance: {delta_a:.6e} m/N")
#     print(f"Fastener Compliance: {delta_b:.6e} m/N")
#     print(f"Force Ratio (Fastener to Total): {ratio:.6f}")

if __name__ == "__main__":
    main()