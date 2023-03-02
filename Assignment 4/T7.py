from matplotlib import pyplot as plt

##Constants#####################################################################
T = np.linspace(100,2000,1000)  #Is this the same as np.linspace(100,2000, 1000)?
kB =1
atm = 101325
beta = (kB * T)**-1

entropies_for_task_7 = np.loadtxt("/content/entropies_for_task_7.txt")

##Parameters for CO and O2
#mass_of_atoms       = np.array([12.011, 15.999])
#mass_of_molecules   = np.array([12.011 + 15.999, 15.999 + 15.999])  #CO & O2 respectively in U
surface_areas       = np.array([67.994, 61.422, 57.465])

adsorbtion_energyCO =  np.array([-0.124, -1.434, -1.586])
adsorbtion_energyO2 =  np.array([ 0.325, -1.041, -1.751])


#adsorbtion_energyCO  = np.array([14.394,15.947,16.012])  #These are placeholder values! 
#adsorbtion_energyO2  = np.array([2.949,4.298, 4.851]) 
activation_energy    =  0.22 -0.3 * (adsorbtion_energyO2 + adsorbtion_energyCO) 


entropy_CO = entropies_for_task_7[0,:]
entropy_O2 = entropies_for_task_7[1,:]

materials = ["Au", "Pt", "Rh"]
#Reaction rate:
v= 10**12 #The so called "Complicated term"



#Pre-allocate lists:
theta_CO      =  np.zeros((3,len(T)))
theta_O       =  np.zeros((3,len(T)))
reaction_rate =  np.zeros((3,len(T)))


for i in range(len(materials)):
  material = materials[i]
  K_CO = np.exp(-entropy_CO[0]/kB)*np.exp(-beta * adsorbtion_energyCO[i])
  K_O2 = np.exp(-entropy_O2[1]/kB)*np.exp(-beta *2* adsorbtion_energyO2[i])

  #Fractional coverage: 

  theta_CO[i,:] =  (K_CO * ((K_O2)**0.5 - (1+K_CO)))/ (K_O2 - (1+ K_CO)**2) 


  theta_O[i,:]  = (K_O2 - np.sqrt(K_O2) * (1 + K_CO))/(K_O2 - (1 + K_CO)**2)
  #Reaction rate: 
  reaction_cross_section = np.multiply(theta_O[i,:],  theta_CO[i,:]) * np.exp(-beta*activation_energy[i])*v

  material_cross_section = surface_areas[i]**-1

  reaction_rate[i,:] = np.dot(reaction_cross_section, material_cross_section)

#Plot results
fig = plt.figure()
fig1 = fig.add_subplot(2,1,1)
for i in range(len(materials)):
    material = materials[i]
    fig1.semilogy(T, reaction_rate[i,:], label="Reaction rate for: "  + material)
fig1.legend()


fig,ax = plt.subplots()

for i in range(len(materials)):
    material = materials[i]
    ax.semilogy(T, theta_CO[i,:], label="Fraction coverage of CO for: "  + material)
    ax.semilogy(T, theta_O[i,:], label="Fraction coverage of O for: "  + material)

ax.legend()