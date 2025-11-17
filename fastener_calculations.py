# A different program should convert load cases of a particular fastener configuration into a list of load cases for each fastner (in the coordinate system of the fastener)

class MaterialProperties(): #class for material properties
    def __init__(self, E = [0,0,0],  rho = 0):
        self.Youngs_modulus = E
        self.density = rho
        self.stress_yield =  0 

#for example some materials
Alumunium = MaterialProperties(E = [70e9,70e9,70e9], rho = 2700) #example material
Steel = MaterialProperties(E = [200e9,200e9,200e9], rho = 7850) #example material

class Load_Case(): #Define one load case
    def __init__(self, Fx = 0, Fy = 0, Fz = 0, Mx = 0, My = 0, Mz = 0):
        self.Force_X = Fx
        self.Force_Y = Fy
        self.Force_Z = Fz
        self.Moment_X = Mx
        self.Moment_Y = My
        self.Moment_Z = Mz

class Fastener_Configuration(): #Define one design configuration
    def __init__(self, d_h = 1, d_f = 1, material = Steel):
        self.fastener_diameter = d_f
        self.head_diameter = d_h
        self.material = material 

class Fastener(): #Define one fastener with its properties and load cases
    def __init__(self, configuration = Fastener_Configuration(), load_cases = [Load_Case()], position = [0,0,0]):
        self.position = position
        self.configuration = configuration
        self.load_cases = load_cases

Total_load_cases = [Load_Case(0,0,0,0,0,0)] #should eventually be a list of all possible load cases

Fastener_positions = [ #defines the positions of the fasteners relative to teh c.g. origin
    [0,0,0],
    [1,0,0],
    [0,0,1],
    [1,0,1]
]

n_f = len(Fastener_positions)
fastener_area = A_i = 0.002 #m^2
area_moment = [0,0,0]
for Fastener_pos in Fastener_positions:
    area_moment[0] += A_i*Fastener_pos[0]**2
    area_moment[1] += A_i*Fastener_pos[1]**2 
    area_moment[2] += A_i*Fastener_pos[2]**2


FastnerList = [] #list of fasteners objects containing material properties, dimensions and Load cases

for pos in Fastener_positions: #This loop goes through each loadcase for the entire structue, calculates the forces on each fastener and stores them in a fastener object
    Local_load_cases = []
    for LC in Total_load_cases:
        F_p = Load_Case(0,0,0,0,0,0)
        F_p.Force_X = 0
        F_p.Force_Y = LC.Force_Y / n_f + LC.Moment_Z*pos[3]*A_i/area_moment[2]
        F_p.Force_Z = 0
        Local_load_cases.append(F_p)
    fastener = Fastener(configuration = Fastener_Configuration(d_h = 0.01, d_f = 0.005, material = Steel), load_cases = Local_load_cases, position = pos)
    FastnerList.append(fastener)

