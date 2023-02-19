##### Initialization #####
import numpy as np
from ase.io.trajectory import Trajectory
from matplotlib import pyplot as plt

plt.rcParams['text.usetex'] = True # Allows for TeX in figures

# Importing the trajectories
traj_AIMD = Trajectory('Data/Dynamics.traj')       # Our AIMD
traj_Na   = Trajectory('Na-aimd/NaCluster24.traj') # Git trajectory w Na+
traj_H2O  = Trajectory('Na-aimd/cluster24.traj')   # Git trajectory w/o Na+


##### Finding the RDF and the first solvation shell #####
def RDF(traj,a_id=72,start=0,stop=None,r1=3,r2=8):
    '''
    Arguments:
        traj  - Trajectory object
        a_id  - ID of atom to be analyzed
        start - Initial step of trajectory
        stop  - Final step of trajectory
        r1,r2 - Inner and outer approx. radii of the first solvation shell

    Returns:
        r      - Radius-array
        gPrime - Array of RDF(r)-values
        SH     - Coordination number
    '''
    if stop == None:
        stop = len(traj) # Default stop is the final snapshot
    
    d = np.zeros([])         # Initial distance array
    V = traj[0].get_volume() # Cell volume
    rho = V/len(traj[0])     # Particle density

    w = traj[0].get_cell_lengths_and_angles()[0] # Cell width

    # Distance from Na+ to O for each timestep
    for i in range(start,stop):
        atoms = traj[i]
        d_i = atoms.get_distances(a_id,
                atoms.get_atomic_numbers()==8,mic=True) # Distances to oxygen
                
        d_i = np.where(d_i<w,d_i,d_i-w) # Applying the boundary condition
        d_i = d_i[d_i != 0] # Disregard distances to the atom itself

        d = np.append(d,d_i) # Add to list
    
    # Finding the RDF g(r)
    hist = np.histogram(d,bins=100) # Data binning
    dr   = hist[1][1] - hist[1][0]
    r    = hist[1][1:] - dr/2
    dn_r = hist[0]/(stop-start)

    gPrime = dn_r*(4*np.pi*r**2*dr*rho)**-1
    
    # Finding the first solvation shell
    minapprox = int((r1-r[0])/dr) # Approximation of first minimum
    min = np.argmin(gPrime[minapprox-r2:minapprox+r2]) + minapprox - r2

    SH = sum(hist[0][:min])/(stop-start) # Solvation shell
     
    return [r,gPrime], SH

# Evaluating the RDFs and the coordination numbers
data,SH         = RDF(traj_AIMD)
data_Na,SH_Na   = RDF(traj_Na,start=3500)
data_H2O,SH_H2O = RDF(traj_H2O,a_id=4,start=2000,r2=5)

print(f'Our Na:    SH = {SH}')
print(f'Their Na:  SH = {SH_Na}')
print(f'Their H20: SH = {SH_H2O}') 


##### Plotting #####
save_figs = False

# Individual figures
fig1,ax1 = plt.subplots(figsize=(7,7)) # Our AIMD
fig2,ax2 = plt.subplots(figsize=(7,7)) # Git AIMD w Na
fig3,ax3 = plt.subplots(figsize=(7,7)) # Git AIMD w/o Na

ax    = [ax1,ax2,ax3]
radii = [data[0],data_Na[0],data_H2O[0]]
rdf   = [data[1],data_Na[1],data_H2O[1]]
label = [r'Our data (Na$^+$)',r'Git data (Na$^+$)',r'Git data (H$_2$O)']

for i in range(3):
    ax[i].plot(radii[i][20:],rdf[i][20:],label=label[i])
    ax[i].tick_params(axis = 'both', labelsize = 15)
    ax[i].set_xlabel(r"Radius [Å]",fontsize=20)
    ax[i].set_ylabel(r"g(r)",fontsize=20)
    ax[i].set_xlim([radii[i][20],radii[i][-1]])
    ax[i].set_ylim([0,1.05*np.amax(rdf[i][20:])])
    ax[i].grid()
    ax[i].legend(fontsize=19)

# Superimposed plot
fig_all,ax_all = plt.subplots(figsize=(10,7))

ax_all.plot(radii[0][20:],rdf[0][20:],label=r'Our data (Na$^+$)')
ax_all.plot(radii[1][20:],rdf[1][20:],label=r'Git data (Na$^+$)')
ax_all.plot(radii[2][20:],rdf[2][20:],label=r'Git data (H$_2$O)')

ax_all.tick_params(axis = 'both', labelsize = 15)
ax_all.set_xlabel(r"Radius [Å]",fontsize=20)
ax_all.set_ylabel(r"g(r)",fontsize=20)
ax_all.set_xlim([radii[0][20],radii[0][-1]])
ax_all.set_ylim([0,0.016])
ax_all.grid()
ax_all.legend(fontsize=19)

# Saving the figures
if save_figs == True:
    fig1.savefig('Figures/RDF our NA.pdf',bbox_inches='tight')
    fig2.savefig('Figures/RDF their Na.pdf',bbox_inches='tight')
    fig3.savefig('Figures/RDF Only H2O.pdf',bbox_inches='tight')
    fig_all.savefig('Figures/RDF superimposed.pdf',bbox_inches='tight')