import numpy as np
from scipy.interpolate import interp1d


material_properties = {"2014-T6": { "density": 2800, "E": 72.4*(10**9), "sigma_y": 414*(10**6), "sigma_u": 483*(10**6)},
                       "7075-T6": { "density": 2810, "E": 71.7*(10**9), "sigma_y": 503*(10**6), "sigma_u": 572*(10**6)},
                       "2014-T2": { "density": 2780, "E": 73.1*(10**9), "sigma_y": 324*(10**6), "sigma_u": 469*(10**6)},
                       "2024-T3": { "density": 2780, "E": 73.1*(10**9), "sigma_y": 345*(10**6), "sigma_u": 483*(10**6)},
                       "2024-T4": { "density": 2780, "E": 71.0*(10**9), "sigma_y": 324*(10**6), "sigma_u": 469*(10**6)},
                       "AZ1916-T6": { "density": 1810, "E": 44.8*(10**9), "sigma_y": 145*(10**6), "sigma_u": 275*(10**6)},
                       "356-T6": { "density": 2680, "E": 72.4*(10**9), "sigma_y": 138*(10**6), "sigma_u": 207*(10**6)},
                       "4130-steel": { "density": 7850, "E": 205*(10**9), "sigma_y": 435*(10**6), "sigma_u": 670*(10**6)},
                       "8630-steel": { "density": 7850, "E": 200*(10**9), "sigma_y": 550*(10**6), "sigma_u": 620*(10**6)}
                 
                        }
                        

def interpolation(filename, i=0, j=1):
    x_vals = []
    y_vals = []

    with open(filename, "r") as f:
        # read all lines except header
        lines = f.readlines()[1:]

    for line in lines:
        # replace decimal comma and semicolon
        line = line.replace(",", ".").replace(";", " ")
        parts = line.strip().split()

        # skip lines that do not contain enough columns
        if len(parts) <= max(i, j):
            continue

        try:
            x = float(parts[i])
            y = float(parts[j])
            x_vals.append(x)
            y_vals.append(y)
        except ValueError:
            # skip lines with invalid numeric entries
            continue

    # convert to numpy arrays
    x = np.array(x_vals)
    y = np.array(y_vals)

    # create interpolation function
    f = interp1d(x, y, kind='cubic', fill_value='extrapolate')

    return f
def Kt( curve, W_D ):

    if curve == "Curve 1":
        f_interp = interpolation("Kt Curves.csv", 0, 1 )
    if curve == "Curve 2":
        f_interp = interpolation("Kt Curves.csv", 2, 3 )
    if curve == "Curve 3":
        f_interp = interpolation("Kt Curves.csv", 4, 5 )
    if curve == "Curve 4":
        f_interp = interpolation("Kt Curves.csv", 6,7 )
    if curve == "Curve 5":
        f_interp = interpolation("Kt Curves.csv", 8, 9 )
    if curve == "Curve 6":
        f_interp = interpolation("Kt Curves.csv", 10, 11 )
    if curve == "Curve 7":
        f_interp = interpolation("Kt Curves.csv", 12, 13 )
    

    K_t = f_interp(W_D)

    return K_t

def K_Bry( t_D, e_D):

    if t_D <= 0.07:
        f_interp = interpolation("K_bry Curves.csv", 0, 1)
        K_bry = f_interp(e_D)
    elif t_D <= 0.09:
        f_interp = interpolation("K_bry Curves.csv", 2, 3)
        K_bry = f_interp(e_D)
    elif t_D <= 0.11:
        f_interp = interpolation("K_bry Curves.csv", 4, 5)
        K_bry = f_interp(e_D)
    elif t_D <= 0.135:
        f_interp = interpolation("K_bry Curves.csv", 6, 7)
        K_bry = f_interp(e_D)
    elif t_D <= 0.175:
        f_interp = interpolation("K_bry Curves.csv", 8, 9)
        K_bry = f_interp(e_D)
    elif t_D <= 0.25:
        f_interp = interpolation("K_bry Curves.csv", 10, 11)
        K_bry = f_interp(e_D)
    elif t_D <= 0.35:
        f_interp = interpolation("K_bry Curves.csv", 12, 13)
        K_bry = f_interp(e_D)
    elif t_D <= 0.5:
        f_interp = interpolation("K_bry Curves.csv", 14, 15)
        K_bry = f_interp(e_D)
    else:
        f_interp = interpolation("K_bry Curves.csv", 16, 17)
        K_bry = f_interp(e_D)

    return K_bry
def K_ty(curve, Av_Abr):

    if curve == "Curve 1":
        f_interp = interpolation("K_ty Curves.csv", 0, 1)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 2":
        f_interp = interpolation("K_ty Curves.csv", 2, 3)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 3":
        f_interp = interpolation("K_ty Curves.csv", 4, 5)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 4":
        f_interp = interpolation("K_ty Curves.csv", 6, 7)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 5":
        f_interp = interpolation("K_ty Curves.csv", 8, 9)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 6":
        f_interp = interpolation("K_ty Curves.csv", 10, 11)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 7":
        f_interp = interpolation("K_ty Curves.csv", 12, 13)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 8":
        f_interp = interpolation("K_ty Curves.csv", 14, 15)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 9":
        f_interp = interpolation("K_ty Curves.csv", 16, 17)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 10":
        f_interp = interpolation("K_ty Curves.csv", 18, 19)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 11":
        f_interp = interpolation("K_ty Curves.csv", 20, 21)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 12":
        f_interp = interpolation("K_ty Curves.csv", 22, 23)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 13":
        f_interp = interpolation("K_ty Curves.csv", 24, 25)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 14":
        f_interp = interpolation("K_ty Curves.csv", 26, 27)
        Kty = f_interp(Av_Abr)
    return Kty




def load_check(material, curve_Kt,flanges,  c, F_xx, F_yy, F_zz, M_xx, M_yy, M_zz, W, D_1, D_p, t_1, e, l, curve_Kty = "Curve 3"):

    Av = 6/(3/(0.5*t_1*(W-D_1*np.sin(np.pi/4)))+1/(0.5*t_1*(W-D_1*np.sin(np.pi/4)))+1/(0.5*t_1*(W-D_1))+1*((e-D_1/2)*t_1))
    ratio = Av/(D_p*t_1) # area ratio

    W_D = W/D_1
    t_D = t_1/D_1
    e_D = e/D_1

    F_tu = material_properties[material]["sigma_u"] #ultimate
    F_ty = material_properties[material]["sigma_y"] #yield

    F_a = F_xx/flanges #forces per flange (x=a, y=b, z=c)
    F_b = F_yy/flanges
    F_c = F_zz/flanges

    M_a = F_c*l + M_xx/flanges #moment per flange
    M_c = F_a*l + M_zz/flanges

    I_xx = (1/12)*t_1*(W**3) #moment of Inertia
    I_zz = (1/12)*W*(t_1**3)
    
    
    sigma_max_root = F_b/(W*t_1) + (M_a*(W/2))/I_xx + (M_c*(t_1/2))/I_zz #max stress at the root
    P_tu = Kt( curve_Kt, W_D)*F_tu*t_1*(W-D_1)*c #force it can take per failure condition
    P_Bry = K_Bry(t_D, e_D)*F_tu*D_p*t_1
    P_ty = K_ty(curve_Kty ,ratio)*F_ty*D_p*t_1


    return F_ty-sigma_max_root, P_tu-F_b, P_bry-F_b, P_ty-F_c

def mass_flange(material, W, D_1, t_1, l, e):

    V = W*l*t_1 + (1/2)*np.pi*(e**2)*t_1 - np.pi*(D_1**2)*t_1
    m = V*material_properties[material]["density"]

    return m
    

#Forces per lug:

F_x = 15000 #N
F_y = 15000 #N
F_z = 15000 #N
M_x = 15 #Nm
M_y = 0 #Nm
M_z = 15 #Nm



#Geometry m

W = 0.03
D_1 = 0.02
D_p = 0.0198 #Diamter of the bolt itself
t_1 = 0.2
e = W/2 # distance center of D_1 to the end of the flange (currently a circular end)
l = 0.1 # flange length

flanges = 2 #number of flanges



#Material Properties

material = "2014-T6" #choose a material from the dictionary material_properties            
curve_Kt = "Curve 1" #Change only the number (based on which material and kinda on geometery)
curve_Kty = "Curve 3"#Change only the number


c = 0.6 #reduction due to ultimate instead or yield (tension)





















