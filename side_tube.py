from pyfluids import Fluid, FluidsList, Input
import numpy as np
import matplotlib.pyplot as plt

def htc_tube_side(mdot_tube,d,T_in_C,T_out_C,p_in_Pa,n):
    #Assume that temperature gradients accross the fully developed flow are small but present in direction of flow
    sections = n
    #Divide L into 100 steps, assume linear temperature distribution
    A = np.pi * (d/2)**2
    mean_htc = 0
    h_list=[]
    for i in range(sections):
        T = T_in_C + i/sections * (T_out_C-T_in_C)
        h2 = Fluid(FluidsList.Hydrogen).with_state(
            Input.pressure(p_in_Pa), Input.temperature(T))
        Re = mdot_tube * d / A / h2.dynamic_viscosity
        
        # Constant temperature boundary conditions are chosen as the capacity ratio is low (C_ratio  = 0.1)
        if Re < 2.3e3:#Laminar
            Nu = 3.657 # From Shah 2003 Table 7.3 on Page 476 for constant wall temp
            #Nu = 4.364 # From Shah 2003 Table 7.3 on Page 476 for constant heat flux
        else:
            Pr = h2.prandtl
            Nu = 0.024 * Re ** 0.8 * Pr ** 0.4 # From Shah 2003 Table 7.6

       
        h=h2.conductivity * Nu / d
        #h = conduct[i] * Nu/d
        h_list.append(h)
        mean_htc += h/sections
        #print(f"Re : {Re:.2e}, local htc: {h:3.0f} W/m2/K, conductivity: {h2.conductivity:.3f} W/m/K, Pr: {h2.prandtl:.2f}")
    
    return mean_htc, h_list

#Test
n = 30
mdot_H2 = 0.063/4 #kg/s
N_tubes = 354
p_in = 100e5 #Pa
tube_ID = (19.05-2*2.11)/1e3 #m
pass_nb = 2
mdot_tube = mdot_H2/N_tubes*pass_nb
A_ratio = 19.05e-3/tube_ID
T_in = 40-273.15 #K
h_mean, h_list = htc_tube_side(mdot_tube,tube_ID,T_in,568.47-273.15,p_in,n)
#print(f"Area averaged htc (outter area based): {h_mean/A_ratio:.2f} W/m2/K")

# From Aspen
htc_aspen = np.array([22.2,28.8,36.9,46.6,57.9,69.4,80.6,92.2,92.7,107,120.6,133.1,133.5,141.9,149.8,154.1,148.6,143.8,139.6,135.8])
bulk_temp_aspen = np.array([40.27,57.31,78.31,103.57,133.17,166.59,202.85,240.79,242.2,292.5,342.1,389.48,390.75,424.34,455.68,484.45,509.96,532.07,551.33,568.19])
k = np.array([])
Pr = np.array([])
Re = np.array([])

G = mdot_tube/ (np.pi/4*tube_ID**2)
for T in bulk_temp_aspen:
    h2 = Fluid(FluidsList.Hydrogen).with_state(
            Input.pressure(p_in), Input.temperature(T-273.15))
    k = np.append(k,[h2.conductivity])
    Pr = np.append(Pr,[h2.prandtl])
    Re = np.append(Re,[G*tube_ID/h2.dynamic_viscosity])

Nu = htc_aspen * tube_ID / k

print(f"Nusselt adjustment factor from Aspen: {np.mean(Nu):.1f} Aspen, books give 3.657; correct by a factor of {np.mean(Nu)/3.657:.3f}")


'''
plt.figure(1)
plt.plot(np.arange(n)/(n-1), np.array(h_list)/A_ratio,"blue",label="Tube side (H2)")
plt.plot([0,1], [274/0.85,274/0.85],"red",label="Shell side (flue)")
plt.ylim(0,550)
#plt.savefig('plots/HTC_tubes_internal_low_dp.svg', format='svg')
plt.show()
'''

'''
plt.figure()

#plt.plot([0,1], [3.67,3.67],label="Model")


St=Nu/Re/Pr
plt.scatter(Re,St*Pr**(2/3),marker='+', label = r"Aspen $j$")
plt.scatter(Re, 3.657/Re*Pr**(-1/3),marker='+',label=r'Aspen $j$ if Nu = cst')
plt.plot(Re, 8.34*Re**(-0.9),linestyle='--',label=r'K&L $f$',color='black')
plt.plot(Re, 8.34/2*Re**(-0.9),linestyle=':',label=r'K&L $f/2$',color='black')
plt.plot(Re, 1.11*Re**(-0.7923),label=r'K&L $j$',color='black')


combined_data = np.column_stack((Re, St*Pr**(2/3)))

# Save to CSV
np.savetxt('arrays.csv', combined_data, delimiter=',', header='Reynolds,j', comments='')

Temps_K = [40,	87.96,	136.09,	184.45,	233.01,	281.3,	329.38,	377.32,	425.09,	472.87,	520.66,	568.47]
conduct = [0.0319,0.0628,	0.0938,	0.1207,	0.1459,	0.1697,	0.1919,	0.2138,	0.2342,	0.2546,	0.274,	0.2932] 
plt.ylabel(r"$j= St \cdot Pr^{2/3}$ (-)")
plt.xlabel(r"Reynolds number")
plt.title("Tube side (H2) Heat Transfer Coefficient Correl")
plt.legend(loc="lower right")
plt.xlim([600,15000])
plt.ylim([0.002,0.03])
plt.yscale('log')
plt.xscale('log')
plt.tight_layout()
plt.grid(True)
#plt.scatter(Re[0]/1e3,St[0]*(3.28)**(2/3))
plt.savefig('plots/St_tubes_internal.svg', format='svg')
plt.show()

'''