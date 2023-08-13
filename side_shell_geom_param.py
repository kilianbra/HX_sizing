import numpy as np
from side_shell import calculate_Nu, calculate_Hg

def calculate_tube_count(Dctl, Ct, pt, c=0,w_p=0):
    """
    Calculate the total number of tubes in an single tube pass shell and tube exchanger.
    Shah 2003 8.5.1
    Parameters:
    Dctl (float): Diameter of the circle through the centers of the outermost tubes.
    Ct (float): Correction factor for tube layout. 0.866 for 30° and 60° tube layouts, 1.00 for 45° and 90° tube layout.
    pt (float): Tube pitch.
    c (float): empty cross sectional area fraction
    c is to be determined by (8.110) if there is an impingement plate on one side or if tube field is removed on both sides

    Returns:
    Nt (float): Total number of tubes in an exchanger.
    """
    Nt = ((np.pi / 4) * Dctl ** 2 * (1 - c) - Dctl*w_p)  / (Ct * pt ** 2) 
    return Nt

def calculate_window_area(Ds, b, lc):
    """
    Calculate the gross window area (i.e., without tubes in the window) or the area of a segment corresponding to the window section.
    
    Parameters:
    Ds (float): Shell-side inside diameter.
    b (float): Angle in radians between two radii intersected at the inside shell wall with the baffle cut.
    lc (float): Clearance.

    Returns:
    Afr_w (float): Gross window area.
    """
    Afr_w =  Ds ** 2 / 4 * (b / 2 - (1 - 2 * lc / Ds) * np.sin(b / 2))
    return Afr_w

def calculate_area_fraction(theta_ctl):
    """
    Calculate the fraction of the number of tubes in one window section encircled by the centerline of the outer tube row.
    
    Parameters:
    theta_ctl (float, radians): angle of chord determining the area fraction to be excluded

    Returns:
    Fw (float): Fraction of the number of tubes in one window section.
    """
    Fw = theta_ctl / (2 * np.pi) - np.sin(theta_ctl) / (2*np.pi) #Shah  8.113
    return Fw


def tube_and_row_numbers(Dctl, c1, c2, pt, w_p, Ds, d_baffle,d_top,d_bot, do, Xl,inline=False):
    """
    Calculate tube and tube row numbers for a shell and tube exchanger.
    
    Parameters:
    - Dctl (float): Diameter of the circle through the centers of the outermost tubes.
    - pt (float): Tube pitch.
    - c1 (float): Area fraction on top where there are no tubes.
    - c2 (float): Area fraction on bottom where there are no tubes.
    - w_p (float): Width of gaps between tubes.
    - Ds (float): Shell-side inside diameter.
    - d_baffle (float): Distance from baffle cut to baffle cut.
    - do (float): Tube outer diameter.
    - Xl (float): Tube longitudinal distance.
    - L_b_c (float): Axial Distance between baffles.

    - Ct (float): Correction factor for tube layout. 0.866 for 30° and 60° tube layouts, 1.00 for 45° and 90° tube layout.
    
    
    Returns:
    tuple: Nt, Nt_2b, Nt_w_top, Nt_w_bot, N_r_cw, N_r_cc, N_t_b
    """
    Ct = 1 if inline else 0.866


    #Calculate the Number of tubes from the top/bottom limits of the tube bundle in Aspen
    Nt = calculate_tube_count(Dctl, Ct=Ct, pt=pt, c=c1+c2)

    theta_ctl = 2*np.arccos(2*d_baffle / Dctl)
    F_w  = calculate_area_fraction(theta_ctl)  #Area fraction of window
    Nt_2b = calculate_tube_count(Dctl, Ct=Ct, pt=pt, c=2*F_w,w_p=w_p)
    F_gap = Dctl*w_p/((np.pi / 4) * Dctl ** 2) #When gap in between middle
    F_bundle =  1 - c1 - c2 - F_gap
    Nt_w_top = (F_w-c1) * Nt / F_bundle #Need to adapt Shah (8.115) as part of window has no tubes, need to account for the empty area fraction c1
    Nt_w_bot =  (F_w-c2) * Nt /F_bundle

    #Find effective number of tubes in cross flow in the window section
    lc_tube_edge_avg = (d_top - d_baffle + d_bot - d_baffle)/2 # Average distance between baffle cut and outter tube row
    lc = (Ds/2-d_baffle) #See Figure 8.9 in Shah, distance from inner shell to baffle cut
    l_c_used = min(0.4*(lc- (Ds-Dctl)/2), lc_tube_edge_avg) # Need to account for the fact that there are only tubes up to a certain point
    N_r_cw = 2*l_c_used/ (Xl*do) #Number of tube rows where cross flow dominates: either effective of real, whatever is lowest. Shah (8.119)

    #Find number of tube rows in pure crossflow (between both baffles)
    N_r_cc = (Ds-2*lc - w_p)/(Xl*do) #Shah (8.121)

    #Number of tubes in one baffle is
    N_t_b = Nt - (Nt_w_bot+Nt_w_top)/2 #Use average as different in top and bottom

    
    return Nt, Nt_2b, Nt_w_top, Nt_w_bot, N_r_cw, N_r_cc, N_t_b


def heat_transfer_coefficient(N_t_b, N_r_cw, N_r_cc,D_otl, Dctl, c1, c2,N_p, w_p, Ds,N_b, d_baffle,L_b_i,L_b_c,d_bot,d_top,delta_sb,delta_tb, do, Xl,Xt, mdot_shell, k_shell, mu_shell, Pr_shell,N_sealpairs):
    """
    Calculate the corrected heat transfer coefficient and related parameters.
    
    Parameters:
    ... [All the input parameters]
    
    Returns:
    tuple: Re_shell, Nu, h_ideal, J_c, J_l, J_b, J_s, J_r, h_corrected
    """
    #Find the cross flow area in between two baffles
    #If 45 deg or 60 deg bundles need to check if pt/do >= 1.707 or 3.7032, if not then need to change (1-1/Xt) to 2*(pt/d0-1)/Xt
    #If tubes have FINS need to modify also!!
    A_o_cr = L_b_c * (Ds - D_otl + Dctl *( 1-1/Xt) ) # Shah (8.122), these Xt are the non dimensional ones, Xt* in Shah
    #Mass velocity/Reynolds numbers
    G_shell = mdot_shell / A_o_cr
    Re_shell = G_shell*(do) / mu_shell

    #Nusselt for crossflow over a bank of tubes (only works for more than 10 rows!)
    Nu = calculate_Nu(Re_shell,Pr_shell,Xl,Xt,False,Nr=N_r_cc)
    h_ideal = Nu * k_shell/ (do)

    theta_top_tubes = 2*np.arccos(2*d_top /Dctl)  #Shah (8.114)
    c1 = calculate_area_fraction(theta_top_tubes) #Area fraction on top where there are no tubes
    theta_bot_tubes = 2*np.arccos(2*d_bot /Dctl)  #Shah (8.114)
    c2 = calculate_area_fraction(theta_bot_tubes) #Area fraction on bottom where there are no tubes

    #Correction for average heat transfer in the window: how many tubes are between baffle tips: 1.0 for no tubes in windowz, 1.15 for small baffle cuts, 0.65 large baffle cuts
    #Near 1.0 for well designed HX
    theta_ctl = 2*np.arccos(2*d_baffle / Dctl)
    F_w  = calculate_area_fraction(theta_ctl)  #Area fraction of window
    #Find area fraction  where just crossflow
    F_gap = Dctl*w_p/((np.pi / 4) * Dctl ** 2) #When gap in between middle
    F_c = 1- 2*F_w  - F_gap #Shah (8.120)
    #Fraction of tubes where just cross flow 
    F_tubes_crossflow = F_c/ (1-c1-c2-F_gap)
    J_c = 0.55 + 0.72 * F_tubes_crossflow #Don't use F_c as this is area fraction, we care about the fraction of tubes in crossflow This differs due to gaps on the oustide with no tubes

    #Axial Leakage through and around baffles (streams A & E), 0.7-0.8 is typical
    A_o_tb = np.pi/4 * ((do+delta_tb)**2-do**2)*N_t_b # Leakage through baffle modified from Shah (8.129)
    #Calculate the area of window (where baffle allows flow to pass)
    lc = (Ds/2-d_baffle) #See Figure 8.9 in Shah, distance from inner shell to baffle cut
    theta_b = 2*np.arccos(1- 2*lc / Ds) # Shah (8.112)
    A_o_sb = np.pi * Ds * delta_sb/2 * (1-theta_b/2/np.pi) # Shah (8.130)
    r_s = A_o_sb / (A_o_sb+A_o_tb)
    r_lm = (A_o_sb + A_o_tb)/A_o_cr
    J_l = 0.44 * (1-r_s) + (1-0.44 * (1-r_s))*np.exp(-2.2*r_lm)

    # bundle partition bypass (streams F & C), 0.7 for large clearances, 0.9 for good clearances or with good clearance strips
    N_ss_plus = N_sealpairs / N_r_cc # number of seal pairs N_ss in Shah, if large enough blocks J_b from being <1
    C = 1.35 if Re_shell <= 100 else 1.25
    r_b = (Ds-D_otl+0.5*N_p*w_p)*L_b_c/A_o_cr # Shah (8.127)
    J_b = 1 if N_ss_plus>=0.5 else np.exp(-C*r_b*(1- (2*N_ss_plus)**(1./3)  ))

    #Correction factor for baffle spacing at inlet and outlet, due to nozzles being further away typically between 0.85 and 1
    L_plus = L_b_i/L_b_c
    #Need to determing if laminar or turbulent? base on Re_s? assume same as 100 limit
    n = 0.6 if Re_shell > 100 else 1./3
    J_s = (N_b - 1 + 2* L_plus**(1-n)) / (N_b - 1 + 2* L_plus)

    # Correction factor for adverse temperature buildup in laminar flows, only for Re_shell <100
    N_r_c = N_r_cc + N_r_cw
    if Re_shell >= 100:
        J_r = 1
    elif Re_shell <= 20:
        J_r = (10/N_r_c)**0.18
    else:
        J_r = (10/N_r_c)**0.18 + (Re_shell - 20)/80 * (1- (10/N_r_c)**0.18 )


    h_corrected=h_ideal*J_b*J_c*J_l*J_r*J_s

    return Re_shell, Nu, h_ideal, J_c, J_l, J_b, J_s, J_r, h_corrected


def pressure_drop(Re_shell,rho_in,N_r_cw, N_r_cc,Nt_w_top,Nt_w_bot,N_t_b,D_otl, Dctl, pt,Xt,Xl,N_p, w_p, Ds, do, L_b_c,L_b_i,N_b,d_baffle,delta_tb,delta_sb, mdot_shell, mu_shell,N_sealpairs):
    """
    Calculate the shell-side pressure drop and related parameters.
    
    Parameters:
    ... [All the input parameters]
    
    Returns:
    tuple: dp_shell, dp_io_percent, dp_window_percent, dp_crossflow_percent, Hg, Re_shell
    """
    ###Pressure drop##
    #Hagen number is a non dimensional pressure drop, avoids velocity dependence but density and viscosity sensistive
    Hg_bundle = calculate_Hg(Re_shell,Xl,Xt,inline=False,Nr=N_r_cc)
    Hg_window = calculate_Hg(Re_shell,Xl,Xt,inline=False,Nr=N_r_cw)

    #Ideal bundle pressure drop
    dp_b_id = mu_shell**2 / rho_in * N_r_cc * Hg_bundle / do**2 #Shah  (6.37)

    #Bypass correction factor (Flows C and F) that avoid tubes in crossflow
    D = 4.5 if Re_shell<=100 else 3.7
    N_r_c = N_r_cc + N_r_cw    
    N_ss_plus =  N_sealpairs/ N_r_c
    A_o_cr = L_b_c * (Ds - D_otl + Dctl *( 1-1/Xt) ) # Shah (8.122), these Xt are the non dimensional ones, Xt* in Shah
    #Bypass area for streams C&F per crossflow section: C is around outside of tubes, F is in between gap at centerline of tubes
    r_b = (Ds-D_otl+0.5*N_p*w_p)*L_b_c/A_o_cr # Shah (8.127)
    zeta_b = 1 if N_ss_plus >= 0.5 else np.exp(-D*r_b*(1-(2*N_ss_plus)**(1./3)))

    #Leakage correction factor (streams A & E) that go around or through the baffle
    A_o_tb = np.pi/4 * ((do+delta_tb)**2-do**2)*N_t_b # Leakage through baffle modified from Shah (8.129)
    #Calculate the area of window (where baffle allows flow to pass)
    lc = (Ds/2-d_baffle) #See Figure 8.9 in Shah, distance from inner shell to baffle cut
    theta_b = 2*np.arccos(1- 2*lc / Ds) # Shah (8.112)
    A_o_sb = np.pi * Ds * delta_sb/2 * (1-theta_b/2/np.pi) # Shah (8.130)
    r_s = A_o_sb / (A_o_sb+A_o_tb)
    r_lm = (A_o_sb + A_o_tb)/A_o_cr
    p=-0.15 * (1+r_s) + 0.8
    zeta_l = np.exp(-1.33*(1+r_s)*r_lm**p)

    #Inlet/outlet correction factor
    #Need to determing if laminar or turbulent? base on Re_s? assume same as 100 limit
    n_prime = 1.0 if Re_shell <= 100 else 0.2
    zeta_s = 2 * (L_b_c/L_b_i) **(2-n_prime)

    A_fr_w = calculate_window_area(Ds,theta_b,lc) # Shah (8.111)


    #Calculate the number of tubes in the window
    A_frontal_tubes_top = np.pi*do**2/4 * Nt_w_top

    A_frontal_tubes_bot = np.pi*do**2/4 * Nt_w_bot
    #Area occupied by tubes in a window on average
    A_fr_t = (A_frontal_tubes_top + A_frontal_tubes_bot)/2 # Take average

    #Find the flow area in the window for shell side flow
    A_o_w = A_fr_w - A_fr_t #Shah (8.117)
    #Find the hydraulic diameter of window section
    theta_ctl = 2*np.arccos(2*d_baffle / Dctl)
    F_w  = calculate_area_fraction(theta_ctl)
    D_h_window = 4 * A_o_w / (np.pi * do * (Nt_w_bot + Nt_w_top)/2  + np.pi*Ds*F_w/2/np.pi) #Modified equation (8.118) to account for top and bottom


    #In the window, the mass flow goes through a different area than in cross flow, so needs a corrected mass velocity
    G_window =mdot_shell*np.sqrt(1/ (A_o_cr*A_o_w) ) #Shah  (6.41) Average mass velocity through window

    if Re_shell>100:
        dp_w_id = (2+0.6*N_r_cw) * G_window**2/2/rho_in #Shah  (6.39a)
    else:
        dp_w_id = 26 * G_window * mu_shell/rho_in * (N_r_cw/(pt-do) + L_b_c/D_h_window**2) + G_window**2 / rho_in #Shah  (6.39b)

    #Shah uses the same cross flow pressure drop and scales it
    #Maybe one should use a cross flow pressure drop accounting for the number of rows in cross flow?
    #This could be included by determining another dp_b_id with Hg_window
    #Hg_window = calculate_Hg(Re_shell,Xl,Xt,inline=False,Nr=N_r_cw)
    dp_i_o = 2 * dp_b_id * (1+ N_r_cw/N_r_cc)*zeta_b*zeta_s #Shah  (6.42)

    dp_shell = zeta_l * ((N_b-1)*dp_b_id*zeta_b + N_b*dp_w_id) + dp_i_o

    dp_io_frac = dp_i_o / dp_shell
    dp_window_frac= zeta_l * ((N_b-1)*dp_b_id*zeta_b) / dp_shell
    dp_crossflow_frac = zeta_l * N_b*dp_w_id / dp_shell
    
    return dp_shell, dp_io_frac, dp_window_frac, dp_crossflow_frac, Hg_bundle, Re_shell, zeta_b, zeta_l, zeta_s



'''

######## INPUTS########
#General HX inputs
mdot_shell = 10/4 #kg/s
k_shell = 0.051 #Flue/air side thermal conductivity in W/m-K Aspen claims varies from 0.0523 at 700K to 0.0496 at 657K
mu_shell = 33e-6 # Pa-s dynamic (mu) viscosity varies from 34 at 700K to 32 at exit from Aspen
Pr_shell = 0.74 # From Aspen
rho_in = np.mean([0.72]) #  inlet 0.72, approx constant goes up to 0.73 at end

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
frac_spacing = 0.22
L_b_i = 3144*(1-frac_spacing)/2/1e3#1229.97/1e3 #axial Distance from flange to first baffle
L_b_c = frac_spacing*3144/1e3#685/1e3 #axial Distance between baffles mm
N_b = 2 #Number of baffles
w_p = 15.88/1e3# Width of gaps between tubes (horizontal pass lane width in Aspen HX Geom then Bundle)
N_p =  1# Number of gaps between tubes
delta_tb = 0.4/1e3 #diameter Clearance between baffle hole and tube outter diameter
delta_sb = 4.76/1e3 ##diameter Clearance between baffle hole and shell inner diameter
N_sealpairs = 2 # Number of sealing strip pairs N_ss in Shah

######## CALCULATIONS for heat transfer  ########
d_baffle =  d_cut2centre # Distance from  baffle cut to baffle cut, from Tubesheet layout in Aspen, effectively D_s/2 - l_c
Dctl = D_otl - do #Definition in Shah p591
#Shah  8.113 to find fraction of area without tubes on top and bottom
theta_top_tubes = 2*np.arccos(2*d_top /Dctl)  #Shah (8.114)
c1 = calculate_area_fraction(theta_top_tubes) #Area fraction on top where there are no tubes
theta_bot_tubes = 2*np.arccos(2*d_bot /Dctl)  #Shah (8.114)
c2 = calculate_area_fraction(theta_bot_tubes) #Area fraction on bottom where there are no tubes


aspen_htc=230.5
aspen_Re_shell = 1.4e4
aspen_dp = 0.087

#Nt, Nt_2b, Nt_w_top, Nt_w_bot, N_r_cw, N_r_cc, N_t_b =tube_and_row_numbers(Dctl, c1, c2, pt, w_p, Ds, d_baffle,d_top,d_bot, do, Xl,inline=False):
Nt, Nt_2b, Nt_w_top, Nt_w_bot, N_r_cw, N_r_cc, N_t_b = 354, 116,101,137,6, (4+6)/2,235
#Number of tubes total, passing through all/2 baffles, in top window, in bottom window, 
#rows in cross flow in window, rows in pure crossflow and tubes passing through one baffle
Re_shell, Nu, h_ideal, J_c, J_l, J_b, J_s, J_r, h_corrected = heat_transfer_coefficient(N_t_b, N_r_cw, N_r_cc,D_otl, Dctl, c1, c2,N_p, w_p, Ds,N_b, d_baffle,L_b_i,L_b_c,d_bot,d_top,delta_sb,delta_tb, do, Xl,Xt, mdot_shell, k_shell, mu_shell, Pr_shell,N_sealpairs)
dp_shell, dp_io_frac, dp_window_frac, dp_crossflow_frac, Hg, Re_shell, zeta_b, zeta_l, zeta_s = pressure_drop(Re_shell,rho_in,N_r_cw, N_r_cc,Nt_w_top,Nt_w_bot,N_t_b,D_otl, Dctl, pt,Xt,Xl,N_p, w_p, Ds, do, L_b_c,L_b_i,N_b,d_baffle,delta_tb,delta_sb, mdot_shell, mu_shell,N_sealpairs)
'''
'''
print(f"Total number of tubes in the exchanger:\t {Nt:.1f}, Aspen: 354")
#Estimate number of tubes passing through baffle overlap area
print(f"tubes passing through both baffles: \t {Nt_2b:.1f}, Aspen: 116")

print(f"Reynolds and Nusselt numbers for SS: \t Nu = {Nu:.2e}; Re = {Re_shell:.2e} (Aspen: {aspen_Re_shell:.2e})")
print(f"Ideal and corrected htc: \t \t {h_ideal:.0f} and {h_corrected:.0f} W/m2/K (Aspen: {aspen_htc:.1f})")
print(f"J coef (htc correction):\t\t J_c = {J_c:.2f} (0.65 - 1.15), J_l= {J_l:.2f} (0.7-0.8), J_b= {J_b:.2f} (0.7-0.9), J_s = {J_s:.2f} (0.85-1), J_r = {J_r:.2f} (<1)")


print(f"zeta coefficients (dp correction): \t z_l {zeta_l:.2f} (0.4-0.5), z_b {zeta_b:.2f} (0.5-0.8), z_s {zeta_s:.2f} (0.5-2)")
print(f"pressure drop shell side: \t\t Dp = {dp_shell/1e5:.3f} bar ({dp_io_frac:.0%} I/O, {dp_crossflow_frac:.0%} xflow," \
      f"{dp_window_frac:.0%} window );(Aspen: {aspen_dp:.3f})")




#Unused dump



'''

'''
# Calculate the gross window area
Afr_w = calculate_window_area(Ds, b, lc)
print(f"Gross window area: {Afr_w} m^2")

# Calculate the fraction of the number of tubes in one window section
Fw = calculate_area_fraction(Dctl)
print(f"Fraction of the number of tubes in one window section: {Fw}")

# Calculate the number of tubes in the window section
Nt_w = Fw * Nt
print(f"Number of tubes in the window section: {Nt_w}")

# Calculate the area occupied by tubes in the window section
Afr_t = np.pi / 4 * do ** 2 * Fw * Nt
print(f"Area occupied by tubes in the window section: {Afr_t} m^2")

# Calculate the net flow area in one window section
Ao_w = Afr_w - Afr_t
print(f"Net flow area in one window section: {Ao_w} m^2")

# Calculate the hydraulic diameter for the window section
Dh_w = 4 * Ao_w / (np.pi * do * Nt_w + np.pi * Ds * (b / 2))
print(f"Hydraulic diameter for the window section: {Dh_w} m")

# Calculate the number of effective tube rows in crossflow in each window
Nr_cw = 0.8 * (17.7 / 1000) / (86.7 / 1000) * (1 / 2 * (Ds - Dctl))
print(f"Number of effective tube rows in crossflow in each window: {Nr_cw}")
'''