# Buckling check
# By: 6126324
#------------------------


#importing packages
import math
import scipy.optimize as opt
#importing other python files
from materials import *
material = MaterialProperties(0,0,0,0)


#defining functions as specified in WP5
def get_euler_column_buckling_stress(params, youngs_modulus = material.Youngs_modulus):
    radius,thickness,length,pressure,applied_force = params
    area_moment_of_inertia = math.pi/64 * ((radius+thickness)**2-radius**2)
    area = 2*math.pi*radius*thickness
    return math.pi**2 * youngs_modulus * area_moment_of_inertia/(area*length**2)

def get_shell_buckling_stress(params, youngs_modulus = material.Youngs_modulus,poisson_ratio = material.poisson_ratio):
    radius, thickness, length, pressure, applied_force = params
    def get_k(_lambda):
        return _lambda + 12/(math.pi**4) * length**4/(radius**2 * thickness**2) * (1-poisson_ratio**2)/_lambda
    k = opt.minimize_scalar(get_k).fun
    q = pressure/youngs_modulus * (radius/thickness)**2
    return (1.983-0.983*math.e**(-23.14*q))*k*math.pi**2*youngs_modulus/(12*(1-poisson_ratio**2))*(thickness/length)**2

# finding safety margin for both types of buckling
def get_euler_column_buckling_safety_margin(params, youngs_modulus = material.Youngs_modulus):
    radius, thickness, length, pressure, applied_force = params
    area_moment_of_inertia = math.pi/64 * ((radius+thickness)**2-radius**2)
    area = 2*math.pi*radius*thickness
    applied_stress = applied_force/area
    return applied_stress/get_euler_column_buckling_stress(params, youngs_modulus)-1

def get_shell_buckling_safety_margin(params, youngs_modulus = material.Youngs_modulus,poisson_ratio = material.poisson_ratio):
    radius, thickness, length, pressure, applied_force = params
    applied_stress = applied_force/(2*math.pi*radius*thickness) # area found using thin-wall assumption
    return applied_stress/get_shell_buckling_stress(params,youngs_modulus,poisson_ratio)-1

# gives all buckling information at once in a dictionary
def get_minimum_buckling_safety_margin(params, youngs_modulus = material.Youngs_modulus,poisson_ratio=material.poisson_ratio):
    radius, thickness, length, pressure, applied_force = params
    area_moment_of_inertia = math.pi / 64 * ((radius + thickness) ** 2 - radius ** 2)
    area = 2 * math.pi * radius * thickness
    euler_column_buckling_stress = get_euler_column_buckling_stress(params, youngs_modulus)
    euler_column_buckling_safety_margin = get_euler_column_buckling_safety_margin(params, youngs_modulus)

    shell_buckling_stress = get_shell_buckling_stress(params, youngs_modulus,poisson_ratio)
    shell_buckling_safety_margin = get_shell_buckling_safety_margin(params, youngs_modulus,poisson_ratio)

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