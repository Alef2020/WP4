import math
import numpy as np

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
        self.shear = shear if shear != 0 else 1.5*self.stress_yield**2 # weird relation as given in WP4, I think its not right yet
        self.specific_cost = 0.01 #assumed cost per kg, could be changed based on material

#example materials
Alumunium = MaterialProperties(E = 70e9, rho = 2700) #example material
Steel = MaterialProperties(E = 200e9, rho = 7850) #example material  

class Fastener_Configuration(): #Define one design configuration
    def __init__(self, Head_diameter = 0.0075, shank_diameter = 0.007, Nut_diameter = 0.004, Length = 1, material = Steel, mass = 0): #standard M4-steel bolt
        self.d_s = shank_diameter
        self.d_h = Head_diameter
        self.d_n = Nut_diameter
        self.L = Length
        self.material = material
        self.mass = mass if mass != 0 else 1.3*(math.pi*(self.d_s/2)**2)*material.density #approximate mass based on volume of cylinder of 3cm length or given mass
        self.cost = 0.5 #assumed cost per fastener, could be changed based on material and size

#example fastener config
M5_steel = Fastener_Configuration(0.0083, 0.0075, 0.005, Steel)
M6_Aluminium = Fastener_Configuration(0.011, 0.009, 0.006, Alumunium)

class Fastener(): #Define one fastener with its properties and load cases
    def __init__(self, position = [0,0,0], configuration = Fastener_Configuration(), load_cases = []):
        self.position = position
        self.configuration = configuration
        self.load_cases = load_cases

class design_Configuration(): #Define one design configuration
    def __init__(self, fastener_positions, load_cases = [], width = 3, height= 2, t_2 = 4, t_3 = 6, Material = Steel, SC_material = Alumunium, fastener_config = Fastener_Configuration):
        #global design parameters
        self.load_cases = load_cases
        
        #Lug plate parameters
        self.t_2 = t_2 #thickness of lug plate
        self.t_3 = t_3 #thickness of mounting plate
        self.lug_width = width #width and height of lug plate
        self.lug_height = height #width and height of lug plate (NOT THICKNESS)
        self.material = Material #material of lug plate
        self.material_wall = SC_material #material of spacecraft wall, could be changed for different

        #fastener design
        self.fastener_positions = fastener_positions
        self.n_f = len(self.fastener_positions)
        self.area_moment = [0,0,0]
        self.fastener_config = fastener_config
        self.A_i = math.pi*(self.fastener_config.d_s/2)**2

        
        #This section of the code calculates the c.g. of the fastener pattern and ajust the positions to make that the origin
        CG = [0,0,0]
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
            self.FastenerList.append(Fastener(Fastener_pos, self.fastener_config))

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

    def mass(self):#function to determine the mass of the fasteners in the design
        fastener_total = sum([f.configuration.mass for f in self.FastenerList])
        Volume_plate = self.t_2*(self.lug_width*self.lug_height - self.n_f*self.A_i)
        plate_mass = Volume_plate * self.material.density
        return fastener_total + plate_mass
    
    def cost(self):
        fastener_cost = self.n_f * self.fastener_config.cost #assumed cost per fastener
        plate_cost = self.mass()*self.material.specific_cost 
        return fastener_cost + plate_cost
    
design = design_Configuration(fastener_positions = [[1,0,1], [-1,0,1], [1,0,-1], [-1,0,-1]], load_cases = [Load_Case(Fx=1000, Fy=500, Fz=2000)], width = 3, height= 2, t_2 = 4, t_3 = 6, Material = Steel, SC_material = Alumunium, fastener_config = M5_steel)