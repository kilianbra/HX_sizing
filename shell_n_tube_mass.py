import numpy as np
from side_shell_geom_param import calculate_area_fraction

def SnT_mass(do=19.05/1e3,L_shell_noflange=3144/1e3,t_tube_single= 2.11/1e3):

    rho_material = 7830 #kg/m3 for carbon steel
    #Basic geometry

    N_t = 354 # Number of tubes total
    #L_tube = 3450/1e3 # Tube length total in m
    #do = 19.05/1e3 #Tube outer diameter in m

    d_i_shell = 600/1e3 #shell inner diameter in m
    t_d_shell = 24/1e3 # Shell diameter/double thickness (i.e. two times the physical thickness)

    d_o_head = 764/1e3 # Head outter diameter
    d_i_head = 600/1e3 # Head inner diameter

    d_o_inlet_ss_nozzle = 406.4/1e3
    h_inlet_ss_nozzle = 149.28/1e3
    d_o_exit_ss_nozzle = 406.4/1e3
    h_exit_ss_nozzle = 149.28/1e3

    #Tubes
    #t_tube_single = 2.11/1e3 #Tube single thickness in m
    t_tubesheet = 149.52/1e3 # Tubesheet thickness in m
    l_tube_extra = 3/1e3 # Protrudes on each end by 3mm
            
    N_b = 2	#Number of baffles (-)
    t_b = 9.52/1e3	#Baffle thickness in m
    delta_tb = 	0.4/1e3	#Baffle-tube clearance (double) m
    delta_sb = 	4.76/1e3 #Baffle-shell clearance (double) m
            
    N_tierod = 6 #Tie rod number (-)
    d_tierod = 	9.55/1e3 #Tie rod diameter m
    N_sealpairs = 2	#Sealing strip pairs (-)
                 
    t_ss_nozzle = 	9.5/1e3	#SS nozzles thickness (single) m
            
    #L_shell_noflange =3144/1e3# Shell flangeless length	m
    t_flange =150.5/1e3 # Flange thickness m
    L_cyl_head = 375.5/1e3 #Cyl (head part) length	m
    w_support =	150/1e3#Shell support width	m
    l_support =537/1e3 # Shell support length m
    h_sup_intube = 497/1e3 #Tube centre to bottom support m
            
    d_baffle = 58.7069/1e3 #Baffle centre to outer baffle cut m
            
    t_seal = 6.35/1e3 #Sealing Strips thickness	m
    w_seal = 16.28/1e3 #Sealing Strips width m
            
#   Nt_2b =	115	#Number of Tubes passing through both baffles 
    N_t_b = 235 # Number of tubes passing through one baffle (avg)

    L_tube = L_shell_noflange + 2* t_tubesheet + 2*l_tube_extra






    #Calculate mass of shell
    pi = np.pi
    m_shell_only = rho_material* pi /4 * ( (d_i_shell+t_d_shell)**2 - d_i_shell**2 ) * L_shell_noflange
    #Ignore the bit of the supports which curves around the shell, assume it is rectangular (small mass)
    m_1shell_support= rho_material * (h_sup_intube-(d_i_shell+t_d_shell)/2) * l_support * w_support

    m_shell_nozzles = rho_material * pi/4 *( ((d_o_inlet_ss_nozzle)**2 - (d_o_inlet_ss_nozzle-t_ss_nozzle)**2) * h_inlet_ss_nozzle + 
                                             ((d_o_exit_ss_nozzle)**2 - (d_o_exit_ss_nozzle-t_ss_nozzle)**2) * h_exit_ss_nozzle ) 
    #The shell flange is hollow, and contains the tubesheet, thicknesses are assumed from others
    perim_flange = pi * d_o_head
    area_flange_from_tubes = pi/4 * (d_o_head**2 - d_i_shell**2)
    m_shell_1flange_hollow = rho_material * (perim_flange * t_flange * t_d_shell + area_flange_from_tubes * t_d_shell/2) 

    m_shell = m_shell_only + 2 * m_1shell_support + m_shell_nozzles + 2 * m_shell_1flange_hollow

    #Calculate mass of the tank heads, which contain a flange, a cylindrical and an elliptical (2:1 ratio) section
    m_head_ellip = rho_material * 0.1309 * (d_o_head**3 - d_i_head**3) # Using a value true for all 2:1 elliptical heads
    m_head_cyl = rho_material * pi/4*  (d_o_head**2 - d_i_head**2) * L_cyl_head
    m_head_flange_solid = rho_material * pi/4*  (d_o_head**2 - d_i_shell**2) * t_flange

    m_1head = m_head_cyl + m_head_ellip + m_head_flange_solid

    #Calculate tube bundle masses
    m_1tube = rho_material * pi/4 * (do**2 - (do-2*t_tube_single)**2) * L_tube
    m_1tubesheet = rho_material * pi/4 * ( (d_o_head-t_d_shell)**2 - N_t * do**2 ) * t_tubesheet

    theta_ctl = 2*np.arccos(2*d_baffle / d_i_shell)
    F_w  = calculate_area_fraction(theta_ctl)  #Area fraction of window
    m_1baffle = rho_material * pi/4 * ( (1-F_w)*(d_i_shell-delta_sb)**2 - N_t_b * (do+delta_tb)**2) * t_b
    m_1sealinstrip = rho_material * t_seal * w_seal * L_tube
    m_1tierod = rho_material * pi/4 * d_tierod**2 * L_tube

    m_bundle = N_t * m_1tube + 2 * m_1tubesheet + N_sealpairs * 2 * m_1sealinstrip + N_tierod * m_1tierod + N_b * m_1baffle

    m_total = m_bundle + 2* m_1head + m_shell

    return m_total, m_bundle, m_1head, m_shell


L_shell_noflange =3144/1e3
m_total, m_bundle, m_1head, m_shell = SnT_mass(L_shell_noflange=L_shell_noflange)
print(f"Mass for one physical HX: \t\t {m_total/1e3:.2f} t ({m_shell/1e3:.2f} t shell, 2x {m_1head/1e3:.2f} t heads,{m_bundle/1e3:.2f} t bundle)")