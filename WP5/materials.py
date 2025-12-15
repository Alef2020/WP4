#material properties
class MaterialProperties(): #class for material properties
    def __init__(self, E,  rho = 0, strength = 0, poisson_ratio = 1/3): #more properties could be added is needed
        self.Youngs_modulus = E
        self.density = rho
        self.stress_yield = strength
        self.poisson_ratio = poisson_ratio



#actual Materials
m2014_T6    = MaterialProperties(72.4*(10**9),2800,414*(10**6))
m7075_T6    = MaterialProperties(71.7*(10**9),2810,503*(10**6))
m2024_T2    = MaterialProperties(73.1*(10**9),2780,324*(10**6))
m2024_T3    = MaterialProperties(73.1*(10**9),2780,345*(10**6))
m2024_T4    = MaterialProperties(71.0*(10**9),2780,324*(10**6))
mAZ1916_T6  = MaterialProperties(44.8*(10**9),1810,145*(10**6))
m356_T6     = MaterialProperties(72.4*(10**9),2680,138*(10**6))
m4130_steel = MaterialProperties(205*(10**9) ,7850,435*(10**6))
m8630_steel = MaterialProperties(200*(10**9) ,7850,550*(10**6))

Material_lst = [
    m2014_T6,
    m2024_T2,
    m2024_T3,
    m2024_T4,
    m7075_T6,
    m4130_steel,
    mAZ1916_T6,
    m356_T6,
    m8630_steel
]
