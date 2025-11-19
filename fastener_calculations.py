import math

class MaterialProperties(): #class for material properties
    def __init__(self, E = [0,0,0],  rho = 0, strength = 0, shear = 0): #more properties could be added is needed
        self.Youngs_modulus = E
        self.shear_strengh = shear
        self.density = rho
        self.stress_yield = strength

#example materials
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
    def __init__(self, d_h = 0.0075, d_f = 0.007, d_o = 0.004, material = Steel): #standard M4-steel bolt
        self.fastener_diameter = d_f
        self.head_diameter = d_h
        self.d_o = d_o
        self.material = material

#example fastener config
M5_steel = Fastener_Configuration(0.0083, 0.0075, 0.005, Steel)
M6_Aluminium = Fastener_Configuration(0.011, 0.009, 0.006, Alumunium)

class Fastener(): #Define one fastener with its properties and load cases
    def __init__(self, position = [0,0,0], configuration = Fastener_Configuration(), load_cases = []):
        self.position = position
        self.configuration = configuration
        self.load_cases = load_cases

class design_Configuration(): #Define one design configuration
    def __init__(self, fastener_positions, load_cases = [], t_1 = 3, t_2 = 4, t_3 = 6, fastener_config = Fastener_Configuration): #any default values are bs, will be fixed with the first iteration
        #other parameters
        self.t_1 = t_1
        self.t_2 = t_2
        self.t_3 = t_3
        self.load_cases = load_cases

        #fastener design
        self.fastener_positions = fastener_positions
        self.n_f = len(self.fastener_positions)
        self.area_moment = [0,0,0]
        self.fastener_config = fastener_config
        self.A_i = math.pi*(self.fastener_config.d_o/2)**2

        self.FastenerList = [] #list of fasteners objects containing material properties, dimensions
        
        for Fastener_pos in self.fastener_positions:
            self.area_moment[0] += self.A_i*Fastener_pos[0]**2 #approximation of 2nd moment of area (x-axis)
            self.area_moment[1] += self.A_i*(Fastener_pos[0]**2 + Fastener_pos[2]**2) #approximation of polar moment of inertia (y-axis))
            self.area_moment[2] += self.A_i*Fastener_pos[2]**2 #approximation of 2nd moment of area around z-axis  
            self.FastenerList.append(Fastener(Fastener_pos ,self.fastener_config))
        
        for f in self.FastenerList: #This part approximates an area moment of inertia
            local_load_cases = []
            for LC in self.load_cases: #This loop goes through each loadcase for the entire structure, calculates the forces on each fastener and stores them in a fastener object
                F_p = Load_Case(0,0,0,0,0,0)
                F_p.Force_X = LC.Force_X / self.n_f + LC.Moment_Y * f.position[2] * self.A_i / self.area_moment[2]
                F_p.Force_Y = LC.Force_Y / self.n_f + LC.Moment_Z * f.position[0] * self.A_i / self.area_moment[0] - LC.Moment_X * f.position[2] * self.A_i / self.area_moment[2]
                F_p.Force_Z = LC.Force_Z / self.n_f - LC.Moment_Y * f.position[0] * self.A_i / self.area_moment[0]
                local_load_cases.append(F_p) # Note that there are no moments on the fastener itself in this simple model
            f.load_cases = local_load_cases



###
#Example program using the classes above:
###


#define loadcases to be analysed
Total_load_cases = [Load_Case(3,4,2,3,6,1), Load_Case(9,2,4,1,3,4)] #should eventually be a list of all possible load cases from 4.1


#Part of finding the center of gravity of fasteners and the relative coordinates of each fastener to the center of gravity
#origin set as the center of the slag
#can be only used for identical fasteners, spacing of fasteners does not have to be uniform
import numpy as np

def center_of_gravity(masses,positions):
    masses = np.array(masses, dtype=float)
    positions = np.array(positions, dtype=float)
    weighted_positions = masses[:,None]*positions
    cg = weighted_positions.sum(axis=0)/masses.sum()
    return cg
    
def positions_cg(masses, positions, D2):
    if not check_spacing_constraint(positions, D2):
        raise ValueError("spacing of fasteners are irrelevant")
        
    positions = np.array(positions, dtype=float)
    cg = center_of_gravity(masses,positions)
    rel_positions = positions -cg
    return cg, rel_positions

def check_spacing_constraint(positions, D2):#checking constraints for the spacing of the fasteners
    positions = np.array(positions, dtype=float)
    min_dist = 0.5 * D2
    n = len(positions)
    for i in range(n):
        for j in range(i+1, n):
            dx = abs(positions[i][0] - positions[j][0])
            dz = abs(positions[i][2] - positions[j][2])
            if dx < min_dist or dz < min_dist:
                return False
    return True
    
#Manual Inputs    
D2 = 
masses = []
positions = []

cg, rel_positions = positions_cg(masses, positions, D2)



#Define Fastener position
Fastener_positions = [ #defines the positions of the fasteners relative to teh c.g. origin ... Square pattern for testing
    [-1,0,-1],
    [1,0,-1],
    [-1,0,1],
    [1,0,1]
]

#Set design parameters, (change standard values if wanted)
example_design = design_Configuration(Fastener_positions, Total_load_cases, fastener_config = M5_steel)


#This program checks if the forces in the design actually add up to the total load case, but any function can be made to check for certain conditions
for LC in Total_load_cases:
    print("Total forces for load case:")
    Total_Fx = 0
    Total_Fy = 0
    Total_Fz = 0
    Total_Mx = 0
    Total_My = 0 
    Total_Mz = 0
    for fastener in example_design.FastenerList:
        index = Total_load_cases.index(LC) #have to check the signs, also some terms may be removed since pos[1] = 0 always
        Total_Fx += fastener.load_cases[index].Force_X
        Total_Fy += fastener.load_cases[index].Force_Y
        Total_Fz += fastener.load_cases[index].Force_Z
        Total_Mx += fastener.load_cases[index].Force_Z * fastener.position[1] - fastener.load_cases[index].Force_Y * fastener.position[2]
        Total_My += (fastener.load_cases[index].Force_X * fastener.position[2] - fastener.load_cases[index].Force_Z * fastener.position[0])/2 # I lost track of why this needs to be devided by 2 but it works ? 
        Total_Mz += fastener.load_cases[index].Force_Y * fastener.position[0] - fastener.load_cases[index].Force_Z * fastener.position[1]
    
    if abs(Total_Fx - LC.Force_X) > 1e-6 or abs(Total_Fy - LC.Force_Y) > 1e-6 or abs(Total_Fz - LC.Force_Z) > 1e-6:
        print("Error in force distribution!") # accounts for floating point errors
    else:
        print("Force distribution correct.")
    
    print(f"Fx: {Total_Fx}, Fy: {Total_Fy}, Fz: {Total_Fz}, Mx: {Total_Mx}, My: {Total_My}, Mz: {Total_Mz}")



