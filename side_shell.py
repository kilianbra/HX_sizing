import numpy as np
import matplotlib.pyplot as plt

#Assumes that there are at least 10 rows of tubes

def calculate_Hglam(Red, Xl, Xt, inline=False):
    #Shah 2003 Fundamentals of HX Equation 7.110
    if inline or (Xl >= 0.5 * (2 * Xt + 1) ** 0.5):
        Hglam = 140 * Red * ((Xl ** 0.5 - 0.6) ** 2 + 0.75) / (Xt ** 1.6 * (4 * Xt * Xl / np.pi - 1))
    else:
        Xd = (Xt ** 2  + Xl ** 2) **0.5
        Hglam = 140 * Red * ((Xl ** 0.5 - 0.6) ** 2 + 0.75) / (Xd ** 1.6 * (4 * Xt * Xl / np.pi - 1))
    test=Hglam
    assert isinstance(test, (int, float)) and not isinstance(test, complex), f"The variable Hglam must be a real number but is {test}"
    return Hglam

def calculate_Hgturb_i(Red, Xt, Xl):
    #Shah 2003 Fundamentals of HX Equation 7.111
    Hgturb_i = ((0.11 + 0.6 * (1 - 0.94 / Xl) ** 0.6 / (Xt - 0.85) ** 1.3) * 10**(0.47 * (Xl / Xt - 1.5)) + 0.015 * (Xt - 1) * (Xl - 1) ) *  Red  **(2-0.1*Xl/Xt)
    test=Hgturb_i
    assert isinstance(test, (int, float)) and not isinstance(test, complex), f"The variable Hgturb_i must be a real number but is {test}"
    return Hgturb_i

def calculate_Hgturb_s(Red, Xt, Xl,Nr):
    #Shah 2003 Fundamentals of HX Equation 7.112
    #Need to correct for number of tube rows
    if Nr>10:
        phi_t_n=0
    else: # Shah Equation (7.114)
        if Xl >= 0.5 * np.sqrt(2*Xt + 1):
            phi_t_n = (1/Nr - 1/10) / (2*Xt**2)
        else:
            Xd = (Xt ** 2  + Xl ** 2) **0.5
            phi_t_n = 2*  (1/Nr - 1/10) *((Xd-1)/(Xt*(Xt-1)))**2

    Hgturb_s = ((1.25 + 0.6 / (Xt - 0.85) ** 1.08) + 0.2 * (Xl / Xt - 1) ** 3 - 0.005 * (Xt / Xl - 1) ** 3) * Red ** 1.75  + phi_t_n * Red ** 2
    if Red>250000:
        Hgturb_s = Hgturb_s * (1 + (Red - 250000) / 325000) #Shah 2003 Fundamentals of HX Equation 7.113
    test=Hgturb_s
    assert isinstance(test, (int, float)) and not isinstance(test, complex), f"The variable Hgturb_s must be a real number but is {test}"
    return Hgturb_s

def calculate_Hg(Red, Xl, Xt, inline=False,Nr=11):
    #Shah 2003 Fundamentals of HX Equation 7.109
    if inline:
        if not (1.25 <= Xt <= 3):
            raise ValueError(f'Xt value {Xt:.2} is not within the range (1.25, 3)')
        if not (1.2 <= Xl <= 3.0):
            raise ValueError(f'Xl value {Xl:.2} is not within the range (1.2, 3.0)')
    else:
        X_d = (Xt ** 2  + Xl ** 2) **0.5
        if not (1.25 <= Xt <= 3):
            raise ValueError(f'Xt value {Xt:.2} is not within the range (1.25, 3)')
        if not (0.6 <= Xl <= 3):
            raise ValueError(f'Xl value {Xl:.2} is not within the range (0.6, 3)')
        if not (5 <= Nr):
            raise ValueError(f'Nr value {Nr:.2} is not greater than 5')
        if X_d is None or X_d <= 1.25:
            raise ValueError(f'X_d value {X_d:.2} is not greater than 1.25')

    Hglam = calculate_Hglam(Red, Xl, Xt, inline)
    if inline:
        Hgturb = calculate_Hgturb_i(Red, Xt, Xl)
    else:
        Hgturb = calculate_Hgturb_s(Red, Xt, Xl,Nr=Nr)
    if inline:
        Hg = Hglam + Hgturb * (1 - np.exp(1 - (Red + 1000) / 2000))
    else:
        Hg = Hglam + Hgturb * (1 - np.exp(1 - (Red + 200) / 1000))

    test=Hg
    assert isinstance(test, (int, float)) and not isinstance(test, complex), f"The variable Hg must be a real number but is {test}"
    
    return Hg

def calculate_Nu(Red, Pr, Xl, Xt, inline=False,Nr=11):
    '''Shah 2003 Fundamentals of HX Equation 7.118
    Flow normal to a tube bundle by Zukauskas 1987 (shell side of shell and tube)
    1< Re_d < 3e5     0.7 < Pr < 700  
    valid (Stag)     1.25 < Xt < 3        0.6 < Xl < 3 and X_d > 1.25 (staggered)
    valid (inline)   1.25 < Xt < 3        1.2 < Xl < 3.0
    Correlations based on 7.9 < d < 73 mm
    Red : Reynolds number based on tube diameter
    Pr: Prandlt number of the fluid
    Xt: tangential (normal to flow) spacing of tubes, normalised by tube outter diameter Xt* in Shah
    Xl: longitudinal (parrallel to flow) spacing of tubes, normalised by tube outter diameter Xl* in Shah
    '''
    Hg = calculate_Hg(Red, Xl, Xt, inline,Nr=Nr)
    if inline:
        Lq = 1.18 * Hg * Pr * (4*Xt/np.pi - 1)/Xl
    elif Xl>= 1:
        Xd = (Xt ** 2  + Xl ** 2) **0.5
        Lq = 0.92 * Hg * Pr * (4*Xt/np.pi- 1)/Xd
    else:
        Xd = (Xt ** 2  + Xl ** 2) **0.5
        Lq = 0.92 * Hg * Pr * (4*Xt*Xl/np.pi - 1)/Xd/Xl
    test=Lq
    assert isinstance(test, (int, float)) and not isinstance(test, complex), f"The variable Lq must be a real number but is {test}"
    
    #Shah 2003 Fundamentals of HX Equation 7.117
    #by Martin 2002
    if inline:
        Nu = 0.404 * Lq **(1./3) * ((Red+1)/(Red+1000))**0.1
    else:
        Nu = 0.404 * Lq **(1./3)
    test = Nu
    assert isinstance(test, (int, float)) and not isinstance(test, complex), f"The variable Nu must be a real number but is {test}"
    return Nu

# For my HX
Xt = 1.25
Xl = 3**0.5 /2
d_o_tube = 19.05/1e3
mdot_flue = 10 /4
A_ocr = 0.1
mu_flue = 0.033e-3
Red = mdot_flue/A_ocr * d_o_tube / mu_flue
#Red = 326
Pr = 0.74
Nu = calculate_Nu(Red,Pr,Xl,Xt,False)
k = 0.051 #Flue/air side thermal conductivity in W/m-K Aspen claims varies from 0.0523 at 700K to 0.0496 at 657K

#For external flow, or in thermal entrance region Nu = hL/k or hD_h/k -> use tube outter diameter for shell side nusselt
h = Nu * k/d_o_tube


Pr_opt = [0.74]#np.linspace(0.6, 1.1, 0.1)  


def calc_Nu2(Red,Pr,inline=True):
    # Heat transfer from Tubes in Crossflow A. Zukauskas 1987
    #Ignoring the wall to bulk temperature gradient
    if inline:
        if 1.6 <= Red <= 40:
            Nu = 0.9 * Red**0.4 * Pr ** 0.36 # Eqn (40)
        elif 40< Red <= 1e3:
            Nu = 0.52 * Red**0.5 * Pr ** 0.36 # Eqn (41)
        else:
            Nu = np.nan
    else:
        if 1.6 <= Red <= 40:
            Nu = 1.04 * Red**0.4 * Pr ** 0.36 # Eqn (38)
        elif 40< Red <= 1e3:
            Nu = 0.71 * Red**0.5 * Pr ** 0.36 # Eqn (39)
        else:
            Nu = np.nan
    return Nu

'''
Xl_i = 1.25
# Calculate the Prandtl number for each pressure
#Compare with Kays and London Inputs
Pr = 0.74
N_r = 15
#Kays and London spacing
Xl = 1.25
Xt = 1.5
Re = np.geomspace(600, 15e3, 50)
j=[]
f=[]
for Red in Re:
        Nu = calculate_Nu(Red,Pr,Xl,Xt,inline=False,Nr=N_r)
        Hg = calculate_Hg(Red, Xl, Xt, inline=False,Nr=N_r)
        j.append( Nu / Red / Pr * Pr**(2/3) )
        f.append(2*(Xt-1)/np.pi * Hg/Red**2 )
plt.plot(Re, j, "k-",label="Heat transfer (K&L)",color='C0')
plt.plot(Re, f, "k--",label="Friction (K&L)",color='C0')
#My Shell and Tube spacing
Re = np.geomspace(600, 15e3, 50)
Xl = 3**(0.5)/2
Xt = 1.25
j=[]
f=[]
for Red in Re:
        Nu = calculate_Nu(Red,Pr,Xl,Xt,inline=False,Nr=N_r)
        Hg = calculate_Hg(Red, Xl, Xt, inline=False,Nr=N_r)
        j.append( Nu / Red / Pr * Pr**(2/3) )
        f.append(2*(Xt-1)/np.pi * Hg/Red**2 )
plt.plot(Re, j, "k-",label="Heat transfer (Aspen geom)",color='C1')
plt.plot(Re, f, "k--",label="Friction (Aspen geom)",color='C1')


# Plot settings
plt.xlabel(r'$Re_d$')
plt.ylabel(r'$j$ and $f$')
plt.title('Tube bank in crossflow Nu correlations')
plt.grid(True)
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.savefig('plots/KL_n_aspen.svg', format='svg')
plt.show()
'''
""" 
plt.figure(1)
for Pr in Pr_opt:
    Nu = []
    Nu_i =[]
    Nu2 = []
    for Red in Re:
        Nu.append(calculate_Nu(Red,Pr,Xl,Xt,inline=False))
        Nu_i.append(calculate_Nu(Red,Pr,Xl_i,Xt,inline=True))
        Nu2.append(calc_Nu2(Red,Pr,inline=False))
    plt.plot(Re, Nu, label=f'Stag. Shah Corr Pr={Pr:.2f} ')
    plt.plot(Re, Nu_i, label=f'Inline ')
    plt.plot(Re, Nu2, label=f'Stag. Zukau Corr Pr={Pr:.2f} ')

mask = (np.array(Re) > 1e4) & (np.array(Re) < 1e6)
Re_filtered = np.array(Re)[mask]
Nu_filtered = np.array(Nu)[mask]

# Find the best fit coefficients
slope, intercept = np.polyfit(np.log10(Re_filtered), np.log10(Nu_filtered), 1)
mask = (np.array(Re) > 1e3) & (np.array(Re) < 1e6)
Re_filtered = np.array(Re)[mask]
plt.plot(Re_filtered, 10**(intercept + slope*np.log10(Re_filtered)), 'k--', label = rf'Best fit: $Nu= {10**intercept:.2f} Re^{{ {slope:.2f} }}$')



#print("done")

#Nu = calculate_Nu(Red,Pr,Xl,Xt,inline=False)
#Nu2 = 1.04 * Red**0.4 * Pr ** 0.36
#print(f"Nu: {Nu:.3e}, Simple Nu:  {Nu2:.3e} at Re {Red:.1e}")
plt.figure(2)

Hg = []
Hg_inline = []
for Red in Re:
    Hg.append(calculate_Hg(Red, Xl, Xt, inline=False))
    Hg_inline.append(calculate_Hg(Red, Xl_i, Xt, inline=True))
# Filter data for the range of interest
mask = (np.array(Re) > 1e4) & (np.array(Re) < 1e6)
Re_filtered = np.array(Re)[mask]
Hg_filtered = np.array(Hg)[mask]

# Find the best fit coefficients
slope, intercept = np.polyfit(np.log10(Re_filtered), np.log10(Hg_filtered), 1)

# Plot the data and best fit line
plt.plot(Re, Hg, label=f'Stag Xl={Xl:.2f}, Xt={Xt:.2f} ')
plt.plot(Re, Hg_inline, label=f'Inline Xl={Xl_i:.2f}, Xt={Xt:.2f}')
mask = (np.array(Re) > 1e3) & (np.array(Re) < 1e6)
Re_filtered = np.array(Re)[mask]
plt.plot(Re_filtered, 10**(intercept + slope*np.log10(Re_filtered)), 'k--', label = rf'Best fit: $Hg= {10**intercept:.2f} Re^{{ {slope:.2f} }}$')


plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$Re_d$')
plt.ylabel(r'Hagen number $Hg$')
plt.legend()
plt.grid(True)
plt.title(r'Tube bank in crossflow $\Delta p$ correlations (Shah)')

#plt.text(5e4,   5e4**slope * 10**intercept, f'Slope: {slope:.3f}\nIntercept: {intercept:.3f}', bbox=dict(facecolor='white', alpha=0.5))

plt.savefig('plots/Hg_crossflow_over_tubes.svg', format='svg')
#plt.show()
 """