#When you add new variables:
#1. Check they don't already exist with a different name
#2. Give it a resonable name
#3. Add units in comment in square brackets
#4. Explain what the variable describes



#-----------------------------
#Set up of variables - beginning
#-----------------------------

#the values below are arbitrary for now
t = 10 #[mm] thickness of the flange, reffer to fig. 4.1
D1 = 10 #[mm] diameter of the hole in flange, reffer to fig. 4.1
E = 7.31e9#[Pa] Young's modulus

#-----------------------------
#Set up of variables - ending
#-----------------------------


#----------------------------
#Example section - beginning
#---------------------------

#calculates the product of thickness of flange and Young's modulus
print(t/E,"In arbitrary units")


#----------------------------
#Example section - ending
#---------------------------


#----------------------------
#4.1 section - beginning
#---------------------------
array_width = 1.303840481
array_length = 1.303840481

sc_moment_arm = [0.50722, 1, 0.0187]
assembly_loads = []
torques = [0.0019215, 0.0018200, 0.000275913]
array_forces = map(zip(torques, sc_moment_arm)torques, lambda x: x[1]/x[2])


# f_z = 6.5*g


#----------------------------
#4.1 section - ending
#---------------------------

