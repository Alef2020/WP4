# heat stress calculator
# assumes linear heat gradient
# 6126324

import math
# thermal stresses
# requires compliance of lug and of bolt from 4.10
# requires stiffness area of bolt from 4.3
#placeholder variables used here:
compliance_lug = 1
compliance_bolt = 1


# finding thermal changes


# defining variables

# known variables
solar_panel_max_temperature = 453
solar_panel_min_temperature = 408

inside_wall_max_temperature = 298
inside_wall_min_temperature = 263

# assumed variables

inner_diameter = 0.01

thermal_expansion_coefficient_fastener = 1e-6
thermal_expansion_coefficient_clamped_part = 1.5e-6

younghs_modulus_bolt = 50e9
younghs_modulus_backplate = 70e9

outer_diameter = 0.3
inner_diameter = 0.2

t_3 = 0.01
t_2 = 0.01

total_thickness_of_attachment = 0.03

solar_panel_beam_length = 0.1


neutral_temperature = 288





stiffness_area = math.pi * (inner_diameter/2)**2






# defining functions
def find_difference(a,b): # works
    return abs(a-b)
def find_total_length(): # works
    return t_3 +t_2 + total_thickness_of_attachment + solar_panel_beam_length
def find_temperature(min_temperature, max_temperature, position):
    return max_temperature - find_difference(min_temperature, max_temperature)*position/find_total_length()

def show_only_outlying_values(input_list):
    output_dict = {}
    if min(input_list) < 0:
        output_dict['Minimum value: '] = min(input_list)
    if max(input_list) > 0:
        output_dict['Maximum values: '] = max(input_list)
    return output_dict

# finding the 4 linearly interpolated temperatures. Any possible temperature will be between these values.
temperature_1 = find_temperature(inside_wall_max_temperature, solar_panel_max_temperature, solar_panel_beam_length)
temperature_2 = find_temperature(inside_wall_min_temperature, solar_panel_min_temperature, solar_panel_beam_length)
temperature_3 = find_temperature(inside_wall_min_temperature, solar_panel_max_temperature, solar_panel_beam_length)
temperature_4 = find_temperature(inside_wall_max_temperature, solar_panel_min_temperature, solar_panel_beam_length)

# finding the maximum and minimum possible temperature
max_temperature = max(temperature_1, temperature_2, temperature_3, temperature_4)
min_temperature = min(temperature_1, temperature_2, temperature_3, temperature_4)

# finding the highest and lowest(most negative) possible temperature differences
highest_deltaT = max_temperature-neutral_temperature
lowest_deltaT = min_temperature-neutral_temperature
deltaT = [highest_deltaT, lowest_deltaT]

# finding the outlying values
deltaT = show_only_outlying_values(deltaT)

#finding max thermal stresses

#compliance_lug = 4*t_2/(younghs_modulus_backplate * math.pi * (outer_diameter**2-inner_diameter**2))

force_ratio = compliance_lug/(compliance_bolt+compliance_lug) # CHECK THIS IS NOT INVERTED

thermal_force_max = (thermal_expansion_coefficient_clamped_part - thermal_expansion_coefficient_fastener) * deltaT.get(
    'Maximum values: ') * younghs_modulus_bolt * stiffness_area * (1 - force_ratio)



#print the max and min values
try:
    thermal_force_max = (thermal_expansion_coefficient_clamped_part-thermal_expansion_coefficient_fastener) * deltaT.get('Maximum values: ') * younghs_modulus_bolt * stiffness_area * (1-force_ratio)
    print("Max force: ", thermal_force_max)
except TypeError:
    pass
try:
    thermal_force_min = (thermal_expansion_coefficient_clamped_part-thermal_expansion_coefficient_fastener) * deltaT.get('Minimum values: ') * younghs_modulus_bolt * stiffness_area * (1-force_ratio)
    print("Min force: ", thermal_force_min)
except TypeError:
    pass