import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import PercentFormatter
from side_shell_geom_param import heat_transfer_coefficient, pressure_drop, calculate_area_fraction
from shell_n_tube_mass import SnT_mass
from side_tube import htc_tube_side
from pyfluids import Fluid, FluidsList, Input
from epsilon_NTU import epsilon_ntu




#General HX inputs
N_HX = 4
mdot_shell = 10/N_HX #kg/s
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

#Tube side inputs
n_tube_steps = 30

mdot_H2 = 0.063/N_HX #kg/s
pass_nb = 2
p_inH2 = 100e5 #Pa
mdot_1tube = mdot_H2/Nt*pass_nb # If 2 pass, only half the tubes get used
t_tube_single = 2.11/1e3 #Tube single thickness in m
tube_ID = do - 2*t_tube_single # m (19.05-2*2.11)/1e3 
Aspen_Nusselt_fudge_factor = 2.033
A_ratio = 19.05e-3/tube_ID

T_in = 40-273.15 #C
c_p_flue = 1130 #J/kg/K from Aspen, shouldn't change
h2_in = Fluid(FluidsList.Hydrogen).with_state(
    Input.pressure(p_inH2), Input.temperature(T_in)
)

q_max_shell_W = mdot_H2 * (h2_in.heating_to_temperature(700-273.15).enthalpy - h2_in.enthalpy)

masses =[]
dps =[]
epsilons =[]
ntus=[]
lengths = np.linspace(0.3,5,50)# [0.3,0.4,0.5,0.6,0.7,1,1.5,2,3,3.144,4,5]
for L_shell_noflange in lengths:
    #Independent length variables
    #L_shell_noflange = 3144/1e3 # in m, distance between to flanges, effective length of tubes for htc
    tubesheet_thickness = 149.52/1e3
    tube_protrusion_past_each_tubesheet = 3/1e3
    L_tube = L_shell_noflange + 2*tubesheet_thickness + 2*tube_protrusion_past_each_tubesheet
    baffle_centre_spacing_frac = 0.22

    #Dependent length variables
    L_b_c = baffle_centre_spacing_frac*L_shell_noflange#685/1e3 #axial Distance between baffles mm
    L_b_i = L_shell_noflange*(1-baffle_centre_spacing_frac)/2#1229.97/1e3 #axial Distance from flange to first baffle
    A_htc =N_HX* Nt * L_shell_noflange * np.pi * do



    Re_shell, Nu, h_ideal, J_c, J_l, J_b, J_s, J_r, h_corrected =heat_transfer_coefficient(N_t_b, N_r_cw, N_r_cc,D_otl, Dctl, c1, c2,N_p, w_p, Ds,N_b, d_baffle,L_b_i,L_b_c,d_bot,d_top,delta_sb,delta_tb, do, Xl,Xt, mdot_shell, k_shell, mu_shell, Pr_shell,N_sealpairs)
    dp_shell, dp_io_frac, dp_window_frac, dp_crossflow_frac, Hg, Re_shell, zeta_b, zeta_l, zeta_s = pressure_drop(Re_shell,rho_in,N_r_cw, N_r_cc,Nt_w_top,Nt_w_bot,N_t_b,D_otl, Dctl, pt,Xt,Xl,N_p, w_p, Ds, do, L_b_c,L_b_i,N_b,d_baffle,delta_tb,delta_sb, mdot_shell, mu_shell,N_sealpairs)

    m_total, m_bundle, m_1head, m_shell = SnT_mass(L_shell_noflange=L_shell_noflange)
    m_HX = N_HX * m_total
    masses.append(m_HX/1e3)

    T_out_guess = 568.47-273.15 #C
    T_out_result = 700-273.15 # C Initialise a random value

    it_max = 100
    it=0
    while abs(T_out_guess - T_out_result)>1 and it<it_max:

        h_mean_i, h_list = htc_tube_side(mdot_1tube,tube_ID,T_in,T_out_guess,p_inH2,n_tube_steps)

        h_i = Aspen_Nusselt_fudge_factor*h_mean_i/A_ratio # Make htc based on outter area

        U = 1/ ( 1/h_i + 1/h_corrected) # Neglect resistance due to tube walls
        c_p_H2 = (h2_in.heating_to_temperature(T_out_guess).enthalpy - h2_in.enthalpy)/(T_out_guess-T_in)
        C_min = mdot_H2  * c_p_H2
        C_ratio = C_min / (mdot_shell*c_p_flue)

        NTU = U * A_htc / C_min /N_HX
        eps = epsilon_ntu(NTU,C_ratio,exchanger_type='shell_and_tube', shell_passes=1)

        enthalpy_out = h2_in.enthalpy + eps*q_max_shell_W / mdot_H2
        T_out_result = h2_in.heating_to_enthalpy(enthalpy_out).temperature
        T_err = T_out_guess - T_out_result
        #print(f"Delta T: {T_err:.1f} K Guessed outlet temperature: {T_out_guess+273.15:.1f} K, actual {T_out_result+273.15:.1f} K")
        it+=1
        T_out_guess = T_out_result + 0.05 * T_err
    dps.append(dp_shell/p_in_flue)
    epsilons.append(eps)
    ntus.append(NTU)
    print(f"mass {m_HX/1e3:.1f}t, A_htc {A_htc:.2f} m2, U={U:.0f}W/m2/K,NTU={NTU:.1f},"\
          f" dp/pin ={dp_shell/p_in_flue:.2%} ({dp_io_frac*dp_shell/p_in_flue:.0%} i/o, {dp_crossflow_frac*dp_shell/p_in_flue:.0%} xflow) ")
    #(c_p_H2 = {c_p_H2/1e3:.1f} kJ/kg/K) eps ={eps:.0%},

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

