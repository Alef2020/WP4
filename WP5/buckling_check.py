# Buckling check
# By: 6126324
#------------------------


#importing packages
import math
import scipy.optimize as opt
#importing other python files
from materials import *
material = Material_lst[1]


#Global parameters
#sets fuel tank volumes:
v_1 = 0.662 #N2O4
v_2 = 1.106 #N2H4

#Other parameters
applied_force = 1
pressure = 100
youngs_modulus = material.Youngs_modulus
poisson_ratio=material.poisson_ratio

#defining functions as specified in WP5
def get_euler_column_buckling_stress(params):
    radius,thickness,length = params
    area_moment_of_inertia = math.pi/64 * ((radius+thickness)**2-radius**2)
    area = 2*math.pi*radius*thickness
    return math.pi**2 * youngs_modulus * area_moment_of_inertia/(area*length**2)

def get_shell_buckling_stress(params):
    radius, thickness, length = params
    def get_k(_lambda):
        return _lambda + 12/(math.pi**4) * length**4/(radius**2 * thickness**2) * (1-poisson_ratio**2)/_lambda
    k = opt.minimize_scalar(get_k).fun
    q = pressure/youngs_modulus * (radius/thickness)**2
    return (1.983-0.983*math.e**(-23.14*q))*k*math.pi**2*youngs_modulus/(12*(1-poisson_ratio**2))*(thickness/length)**2

# finding safety margin for both types of buckling
def get_euler_column_buckling_safety_margin(params):
    radius, thickness, length = params
    area_moment_of_inertia = math.pi/64 * ((radius+thickness)**2-radius**2)
    area = 2*math.pi*radius*thickness
    applied_stress = applied_force/area
    return applied_stress/get_euler_column_buckling_stress(params)-1

def get_shell_buckling_safety_margin(params):
    radius, thickness, length = params
    applied_stress = applied_force/(2*math.pi*radius*thickness) # area found using thin-wall assumption
    return applied_stress/get_shell_buckling_stress(params)-1

# gives all buckling information at once in a dictionary
def get_minimum_buckling_safety_margin(params):
    radius, thickness, length = params
    area_moment_of_inertia = math.pi / 64 * ((radius + thickness) ** 2 - radius ** 2)
    area = 2 * math.pi * radius * thickness
    euler_column_buckling_stress = get_euler_column_buckling_stress(params)
    euler_column_buckling_safety_margin = get_euler_column_buckling_safety_margin(params)

    shell_buckling_stress = get_shell_buckling_stress(params)
    shell_buckling_safety_margin = get_shell_buckling_safety_margin(params)

    minimum_buckling_stress = min(euler_column_buckling_stress,shell_buckling_stress)
    minimum_buckling_safety_margin = min(euler_column_buckling_safety_margin,shell_buckling_safety_margin)
    buckling_dictionary = {
        "euler column buckling stress": euler_column_buckling_stress,
        "shell buckling stress": shell_buckling_stress,
        "minimum buckling stress": minimum_buckling_stress,
        "euler column buckling safety margin": euler_column_buckling_safety_margin,
        "shell buckling safety margin": shell_buckling_safety_margin,
        "minimum buckling safety margin": minimum_buckling_safety_margin
    }

    #prints minimum safety margin in color
    if minimum_buckling_safety_margin < 0:
        print("\033[91mMinimum safety margin: " + str(minimum_buckling_safety_margin))
    else:
        print("\033[92mMinimum safety margin " + str(minimum_buckling_safety_margin))
    return buckling_dictionary


#Volume constraints
def tank_integration_constraint(params):
    R, t, L = params
    return v_2-(4*math.pi*(R**3)/3) #checks if largest fuel tank can be integrated in to a shell structure with radius R

def volume_constraint(params):
    R, t, L = params
    return (L*math.pi*(R**2) - 2*math.pi*(R**3)/3)-v_1-v_2 #checks if the total volume in the shell is larger than the the fuel tanks (also accounting for spherical caps that reduce volume efficiency)
