import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import PercentFormatter
from side_shell_geom_param import heat_transfer_coefficient, pressure_drop, calculate_area_fraction
from shell_n_tube_mass import SnT_mass
from side_tube import htc_tube_side
from pyfluids import Fluid, FluidsList, Input
from epsilon_NTU import epsilon_ntu
from side_shell import calculate_Nu, calculate_Hg




#General HX inputs


k_shell = 0.051 #Flue/air side thermal conductivity in W/m-K Aspen claims varies from 0.0523 at 700K to 0.0496 at 657K
mu_shell = 33e-6 # Pa-s dynamic (mu) viscosity varies from 34 at 700K to 32 at exit from Aspen
Pr_shell = 0.74 # From Aspen
rho_in = np.mean([0.72]) #  inlet 0.72, approx constant goes up to 0.73 at end
p_in_flue = 1.5e5

#Values for my HX - all from Aspen
Ds =  0.6 #Shell inner diameter
d_cut2centre = 58.7/1e3 #Baffle cut to  centre, related to l_c in Shah
do = 19.05/1e3 #Tube outter diameter
D_otl = 0.5873 #Outter tube limit
pt = 1.25 * do # tube pitch (centre to centre)

d_top = 141.2/1e3 #D_s/2 - l_c in Shah: Distance from on top tube row to centre of shell, from tubesheet layout is Aspen
d_bot = 182.44/1e3 #Same for lower tube row


#New inpunts, from side_shell also (tube distances)
Xl = 3**0.5 /2 #Tube longitudinal (parrallel to flow) distance - non dimensional Xl* in Shah
Xt = 1.25 #Tube transverse (normal to flow) distance - non dimensional Xt* in Shah

N_b = 2 #Number of baffles
w_p = 15.88/1e3# Width of gaps between tubes (horizontal pass lane width in Aspen HX Geom then Bundle)
N_p =  1# Number of gaps between tubes
delta_tb = 0.4/1e3 #diameter Clearance between baffle hole and tube outter diameter
delta_sb = 4.76/1e3 ##diameter Clearance between baffle hole and shell inner diameter
N_sealpairs = 2 # Number of sealing strip pairs

######## CALCULATIONS for heat transfer  ########
d_baffle =  d_cut2centre # Distance from  baffle cut to baffle cut, from Tubesheet layout in Aspen, effectively D_s/2 - l_c
Dctl = D_otl - do #Definition in Shah p591
#Shah  8.113 to find fraction of area without tubes on top and bottom
theta_top_tubes = 2*np.arccos(2*d_top /Dctl)  #Shah (8.114)
c1 = calculate_area_fraction(theta_top_tubes) #Area fraction on top where there are no tubes
theta_bot_tubes = 2*np.arccos(2*d_bot /Dctl)  #Shah (8.114)
c2 = calculate_area_fraction(theta_bot_tubes) #Area fraction on bottom where there are no tubes

Nt, Nt_2b, Nt_w_top, Nt_w_bot, N_r_cw, N_r_cc, N_t_b = 354, 116,101,137,6, (4+6)/2,235
#Number of tubes total, passing through all/2 baffles, in top window, in bottom window, 
#rows in cross flow in window, rows in pure crossflow and tubes passing through one baffle

#Independent length variables
L_shell_noflange = 3144/1e3 # in m, distance between to flanges, effective length of tubes for htc
tubesheet_thickness = 149.52/1e3
tube_protrusion_past_each_tubesheet = 3/1e3
L_tube = L_shell_noflange + 2*tubesheet_thickness + 2*tube_protrusion_past_each_tubesheet
baffle_centre_spacing_frac = 0.22

#Dependent length variables
L_b_c = baffle_centre_spacing_frac*L_shell_noflange#685/1e3 #axial Distance between baffles mm
L_b_i = L_shell_noflange*(1-baffle_centre_spacing_frac)/2#1229.97/1e3 #axial Distance from flange to first baffle


c_p_flue = 1130 #J/kg/K from Aspen, shouldn't change

Re =[]
j_corr = []
f_corr = []

#mdot_shell = 10/N_HX #kg/s
massflows = np.geomspace(0.1,10,50)# [0.3,0.4,0.5,0.6,0.7,1,1.5,2,3,3.144,4,5]
for mdot_shell in massflows:

    Re_shell, Nu, h_ideal, J_c, J_l, J_b, J_s, J_r, h_corrected =heat_transfer_coefficient(N_t_b, N_r_cw, N_r_cc,D_otl, Dctl, c1, c2,N_p, w_p, Ds,N_b, d_baffle,L_b_i,L_b_c,d_bot,d_top,delta_sb,delta_tb, do, Xl,Xt, mdot_shell, k_shell, mu_shell, Pr_shell,N_sealpairs)
    dp_shell, dp_io_frac, dp_window_frac, dp_crossflow_frac, Hg, Re_shell, zeta_b, zeta_l, zeta_s = pressure_drop(Re_shell,rho_in,N_r_cw, N_r_cc,Nt_w_top,Nt_w_bot,N_t_b,D_otl, Dctl, pt,Xt,Xl,N_p, w_p, Ds, do, L_b_c,L_b_i,N_b,d_baffle,delta_tb,delta_sb, mdot_shell, mu_shell,N_sealpairs)
    Re.append(Re_shell)
    G = Re_shell * mu_shell / do
    #j_corr.append(h_corrected * do * Pr_shell**(2/3) /G/c_p_flue )
    j_corr.append(h_corrected * Pr_shell**(2/3) /G/c_p_flue )
    A_htc = Nt * L_shell_noflange * np.pi * do
    fudge_factor = 28.885/24.196
    A_o_cr = fudge_factor * L_b_c * (Ds - D_otl + Dctl *( 1-1/Xt) )
    A_htc_2_A_cross_ratio =A_htc/A_o_cr
    #A_htc_2_A_cross_ratio = 226 / 0.1
    f_corr.append(2*rho_in * dp_shell/G**2/A_htc_2_A_cross_ratio)


Re = np.array(Re)

#My Shell and Tube spacing
#Re = np.geomspace(600, 15e3, 50)

Xl = 3**(0.5)/2
Xt = 1.25
Pr=Pr_shell
N_r=10
j=[]
f=[]
for Red in Re:
        Nu = calculate_Nu(Red,Pr,Xl,Xt,inline=False,Nr=N_r)
        Hg = calculate_Hg(Red, Xl, Xt, inline=False,Nr=N_r)
        j.append( Nu / Red / Pr * Pr**(2/3) )
        f.append(2*(Xt-1)/np.pi * Hg/Red**2 )
#plt.plot(Re, j, "k-",label="Heat transfer (ideal)",color='C1')
#plt.plot(Re, f, "k--",label="Friction (ideal)",color='C1')
#plt.plot(Re, j_corr, "k-",label="Heat transfer (corr.)",color='C2')
#plt.plot(Re, f_corr, "k--",label="Friction (corr.)",color='C2')
m_corr =  np.array(f_corr)/2/np.array(j_corr)
m = np.array(f)/2/np.array(j)

print(f"Tube bank: {np.mean(m[Re > 1e4]):.2f}, Shell and Tube: {np.mean(m_corr[Re > 1e4]):.2f}")
plt.plot(Re, m_corr,label="Friction per heat transfer corr.")
plt.plot(Re, m,label="Friction per heat transfer")

# Plot settings
plt.xlabel(r'$Re_d$')
plt.ylabel(r'$j$ and $f$')
plt.title('Tube Bank in Crossflow (ideal) vs Shell and Tube (B-D corr.)')
plt.grid(True)
plt.legend()

for target in [2.5,10/3,5,10]:
    index = np.argmin(np.abs(massflows - target))
    plt.plot([Re[index],Re[index]],[2e-3,2.2e-3],"k--")

plt.xscale('log')
plt.yscale('log')
#plt.ylim([2e-3,0.1])
#plt.savefig('plots/Corr_uncorr_.svg', format='svg')
plt.show()



'''
# [0.3,0.4,0.5,0.6,0.7,1,1.5,2,3,3.144,4,5]
aspen_lengths = [4, 3.145,3,2.5,2,1.5,1]
aspen_mass = [26.3,24.7,24.4,23.5,22.5,21.59,20.5]
aspen_dp = np.array([4.7,5.91,6.12,7.54,8.21,11.73,19.47])/100
aspen_NTU = [24.2,20.1,19.3,16.7,13.6,10.6,7.2]


# Create figure and axis
fig, ax1 = plt.subplots()

# Plot data on the left y-axis
#l1= ax1.scatter(lengths, epsilons, marker='+', label=r'$\varepsilon=q/q_{max}$')
#l2 =ax1.scatter(lengths, dps, marker='x',color='C1', label=r'$\Delta p/p_{in}$ ')
l2, =ax1.plot(lengths, dps, color='C1', label=r'$\Delta p/p_{in}$ ')
ax1.set_ylabel('Pressure drop', color='C1')
ax1.tick_params(axis='y', labelcolor='black')

ax1.set_ylim(0,1)
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
# Create a secondary y-axis for masses
ax2 = ax1.twinx()
#l1=ax2.scatter(lengths,ntus, marker='+',color='C0', label=r'NTU (-) ')
l1,=ax2.plot(lengths,ntus, color='C0', label=r'$N_{tu}\propto UA$ (-) ')
#l3 = ax2.scatter(lengths, masses, marker='o', color='grey', facecolors='grey', label='Mass (t)')
l3, = ax2.plot(lengths, masses, color='grey', label='Mass (t)')
#ax2.set_ylabel('Mass (t) and NTU(-)')
# Manually set the two parts of the label with different colors
ax2.text(1.1, 0.3, 'Mass (t)', transform=ax2.transAxes, verticalalignment='center', color='grey', rotation=90)
ax2.text(1.1, 0.45, ' and ', transform=ax2.transAxes, verticalalignment='center', rotation=90)
ax2.text(1.1, 0.6, 'NTU (-)', transform=ax2.transAxes, verticalalignment='center', color='C0', rotation=90)
#ax2.tick_params(axis='y', labelcolor='grey')
ax2.set_ylim(0,27)
ax2.set_xlim(0,5)
l4=ax2.scatter(aspen_lengths,aspen_mass,marker='+',color='grey',label='Aspen validation')

lines= [l1,l2,l3,l4]
labels =[l.get_label() for l in lines]
ax1.legend(lines, labels, loc='center right')
ax1.scatter(aspen_lengths,aspen_dp,marker='+',color='C1')

ax2.scatter(aspen_lengths,aspen_NTU,marker='+',color='C0')
# Set x-axis label
ax1.set_xlabel('Tube length (m)')
fig.suptitle("Shell and Tube shortening: Impact on performance")
fig.tight_layout()
plt.savefig('plots/Lengthening.svg', format='svg')
plt.show()



#print(f"Max heat transfer {N_HX*q_max_shell_W/1e3:.2f} kW, U={U:.0f}W/m2/K, A={A_htc:.0f}m2, NTU = {NTU:.1f} eps ={eps:.2%}")

'''