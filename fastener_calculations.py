import math

class Load_Case(): #Define one load, this is more for clarity, it could also just be an array
    def __init__(self, Fx = 0, Fy = 0, Fz = 0, Mx = 0, My = 0, Mz = 0):
        self.Force_X = Fx
        self.Force_Y = Fy
        self.Force_Z = Fz
        self.Moment_X = Mx
        self.Moment_Y = My
        self.Moment_Z = Mz    

class MaterialProperties(): #class for material properties
    def __init__(self, E,  rho = 0, strength = 0, shear = 0): #more properties could be added is needed
        self.Youngs_modulus = E
        self.density = rho
        self.stress_yield = strength
        self.shear_strenght = shear if shear != 0 else 1.5*self.stress_yield**2 # weird relation as given in WP4, I think its not right yet

#example materials
Alumunium = MaterialProperties(E = 70e9, rho = 2700) #example material
Steel = MaterialProperties(E = 200e9, rho = 7850) #example material  

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
    def __init__(self, fastener_positions, load_cases = [], t_1 = 3, t_2 = 4, t_3 = 6, fastener_config = Fastener_Configuration):
        #other parameters these can be changed for other purposes
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

        #This section of the code calculates the c.g. of the fastener pattern and ajust the positions to make that the origin
        CG = [0,0.0]
        CG[0] = sum([pos[0] for pos in self.fastener_positions])/self.n_f #can prob be done better with np
        CG[2] = sum([pos[2] for pos in self.fastener_positions])/self.n_f
        for Fastner_pos in self.fastener_positions:
            Fastner_pos[0] -= CG[0]
            Fastner_pos[2] -= CG[2]
        
        for LC in self.load_cases:
            LC.Moment_X = -CG[2]*LC.Force_Y 
            LC.Moment_Z = CG[0]*LC.Force_Y #adjust moments based on new origin

        #This section calculates the area moment of inertia (sort of; the equivelant calculation to get the distribution of forces on each fastener) and creates fastener objects for each position
        self.FastenerList = [] #list of fasteners objects containing material properties, dimensions
        for Fastener_pos in self.fastener_positions:
            self.area_moment[0] += self.A_i*Fastener_pos[0]**2 #approximation of 2nd moment of area (x-axis)
            self.area_moment[1] += self.A_i*(Fastener_pos[0]**2 + Fastener_pos[2]**2) #approximation of polar moment of inertia (y-axis))
            self.area_moment[2] += self.A_i*Fastener_pos[2]**2 #approximation of 2nd moment of area around z-axis  
            self.FastenerList.append(Fastener(Fastener_pos ,self.fastener_config))

        #This section goes through each loadcase for the entire structure, calculates the forces on each fastener and stores them in a fastener object
        for f in self.FastenerList: 
            local_load_cases = []
            for LC in self.load_cases: 
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



import numpy as np

#function of checking the spacing constraint of fasteners
def check_spacing_constraint(fastener_positions, D2):#checking constraints for the spacing of the fasteners
    positions = np.array(fastener_positions, dtype=float)
    min_dist = 2 * D2
    max_dist = 3 * D2
    n = len(positions)
    for i in range(n):
        for j in range(i+1, n):
            
            dz = abs(positions[i][2] - positions[j][2])
            if  dz < min_dist  or dz > max_dist:
                return False
    return True
       




#Define Fastener position
Fastener_positions = [ #defines the positions of the fasteners relative to teh c.g. origin ... Square pattern for testing
    [-1,0,-1],
    [1,0,-1],
    [-1,0,1],
    [1,0,1]
]

#Set design parameters, (change standard values if wanted)
example_design = design_Configuration(Fastener_positions, Total_load_cases, fastener_config = M5_steel)


def pull_through_check(design_option):
            for fastener in design_option.FastenerList:
                for LC in fastener.load_cases:
                    A_bearing = math.pi * (design_option.fastener_config.head_diameter**2 - design_option.fastener_config.d_o**2) /4
                    bearing_stress = LC.Force_Y / A_bearing
                    if bearing_stress > design_option.fastener_config.material.stress_yield:
                        return False
            return True


for width in range(10, 20, 0.1):
    for height in range(10, 20, 0.1):
        e =0.3
        fastener_positions = [
            [ width/2-e, 0,  height/2-e],
            [-width/2+e, 0,  height/2-e],
            [ width/2-e, 0, -height/2+e],
            [-width/2+e, 0, -height/2+e],
            [0, 0, -height/2+e],
            [0, 0, -height/2+e]
        ]

        design_option = design_Configuration(fastener_positions, Total_load_cases, fastener_config = M5_steel)
        if pull_through_check(design_option):
            print(f"Design passes for width: {width} cm and height: {height} cm")










