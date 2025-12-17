from math import pi
import numpy as np
import pprint

##############################################
# Vars
##############################################

material_properties = {
    "7075-T6": {
        "density": 2810,
        "E": 71.7 * 10**9,
        "sigma_y": 503 * 10**6,
        "sigma_u": 572 * 10**6,
        "cost/kg": 10,
    },
    "AZ31B H24": {
        "density": 1770,
        "E": 165 * 10**9,
        "sigma_y": 200 * 10**6,
        "sigma_u": None,  # not provided
        "cost/kg": 42.52,  # 50 USD
    },
    "Ti-6Al-4V": {
        "density": 4430,
        "E": 3 * 10**9,
        "sigma_y": 855 * 10**6,
        "sigma_u": None,
        "cost/kg": 38.27,  # 45 USD
    },
    "AMS 7906": {
        "density": 1850,
        "E": 5 * 10**9,
        "sigma_y": 241 * 10**6,
        "sigma_u": None,
    },
    "T800H/epoxy [0,±45,90]": {
        "density": 1600,
        "E": 4 * 10**9,
        "sigma_y": 660 * 10**6,
        "sigma_u": None,
    },
    "Kevlar49/epoxy [0,±45,90]": {
        "density": 1400,
        "E": 5 * 10**9,
        "sigma_y": 150 * 10**6,
        "sigma_u": None,
    },
}
R = 0.597
L = 5.4
internal_pressure = 500000
# hoop_stress = pr/t
# sigma_long = pr/2t
##############################################
# Constraints
##############################################


def calc_thickness(
    material_name,
    material_yield_stress_long_direction,
    material_yield_stress_hoop_direction,
):
    t1_long = (internal_pressure * R) / (material_yield_stress_long_direction * 2)

    t1_hoop = (internal_pressure * R) / material_yield_stress_hoop_direction
    t1 = t1_long if t1_long > t1_hoop else t1_hoop

    return t1


def calc_mass_cylinder(
    thickness,
    material_density,
    material_cost_per_kilogram,
):
    cylinder_mass = material_density * pi * (R**2 - (R - thickness) ** 2) * L
    end_caps_mass = material_density * 4 / 3 * pi * (R**3 - (R - thickness) ** 3)
    total_mass = cylinder_mass + end_caps_mass
    return total_mass


def calc_costs(mass, material_cost_per_kilogram):
    return mass * material_cost_per_kilogram


##############################################
# Results
##############################################

results = [["material", "thickness", "mass", "costs"]]
for material_name, material_data in material_properties.items():
    material_yield_stress_hoop_direction = material_data["sigma_y"]
    material_yield_stress_long_direction = material_data["sigma_y"]
    material_density = material_data["density"]
    material_cost_per_kilogram = (
        material_data["cost/kg"] if "cost/kg" in material_data else -1000
    )

    thickness_internal_pressure = calc_thickness(
        material_name,
        material_yield_stress_long_direction,
        material_yield_stress_hoop_direction,
    )
    thickness = min(thickness_internal_pressure, 1000)
    mass = calc_mass_cylinder(thickness, material_density, material_cost_per_kilogram)
    costs = calc_costs(mass, material_cost_per_kilogram)
    results.append([material_name, thickness, mass, costs])
pprint.pprint(results)
