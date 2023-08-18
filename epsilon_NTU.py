import numpy as np
import matplotlib.pyplot as plt

def epsilon_ntu(NTU, C_ratio, exchanger_type='double_pipe', flow_type='parallel', shell_passes=1):
    """
    Calculate the effectiveness (epsilon) of a heat exchanger using the Number of Transfer Units (NTU) method.
    
    Parameters:
    NTU (float/np array): Number of Transfer Units
    C_ratio (float): Capacity rate ratio 0 <= C_min/C_max <= 1
    exchanger_type (str, optional): Type of the heat exchanger. Can be 'double_pipe', 'cross_flow', or 'shell_and_tube'. Default is 'double_pipe'.
    flow_type (str, optional): Arrangement of the flow. Can be 'parallel', 'counter', 'unmixed', 'Cmax_mixed', 'Cmin_mixed', or 'both_mixed'. Default is 'parallel'.
    shell_passes (int, optional): Number of shell passes in case of shell-and-tube exchanger. Default is 1.
    epsilon_1p (float, optional): Effectiveness of each case in case of shell-and-tube exchanger. Default is 1.
    
    Returns:
    epsilon (float): Effectiveness of the heat exchanger
    """
    assert 0 <= C_ratio <= 1, f"C_ratio should be between 0 and 1 (inclusive), but is {C_ratio:.1f}"
    assert isinstance(shell_passes, int) and shell_passes >= 1, f"shell_passes should be an integer greater than or equal to 1 but is {shell_passes:.0f}"
    valid_exchanger_types = ['double_pipe', 'cross_flow', 'shell_and_tube']
    valid_flow_types = ['parallel', 'counter', 'unmixed', 'Cmax_mixed', 'Cmin_mixed', 'both_mixed']
    assert exchanger_type in valid_exchanger_types, f"Invalid exchanger_type. Must be one of {valid_exchanger_types}"
    assert flow_type in valid_flow_types, f"Invalid flow_type. Must be one of {valid_flow_types}"



    tol = 1e-9
    if C_ratio < 0 + tol:  #Close enough to zero, doesn't matter the type as C_min fluid doesn't change temperature
        epsilon = 1 - np.exp(-NTU)  # Kays & London (2-13a)
    
    elif exchanger_type == 'double_pipe': #Two concentric tubes aligned with eachother
        if flow_type == 'parallel':
            epsilon = (1 - np.exp(-NTU * (1 + C_ratio))) / (1 + C_ratio) # Kays & London (2-14)
        elif flow_type == 'counter':
            if C_ratio < 1 - tol: #Far enough from 1
                epsilon = (1 - np.exp(-NTU * (1 - C_ratio))) / (1 - C_ratio * np.exp(-NTU * (1 - C_ratio))) # Kays & London (2-13)
            else: #Close enough to 1
                epsilon = NTU / (NTU + 1) # Kays & London (2-13b)

    elif exchanger_type == 'cross_flow':
        if flow_type == 'unmixed':
            epsilon = 1 - np.exp(NTU**0.22 / C_ratio * (np.exp(-C * NTU**0.78) - 1)) #CHECK REFERENCE Eqn (3.24) in Lopez
        elif flow_type == 'Cmax_mixed':
            #epsilon = 1 / C_ratio * (1 - np.exp(1 - C_ratio) * (1 - np.exp(-NTU))) # Lopez (3.25) WRONG
            epsilon = 1 / C_ratio * (1 - np.exp(-C_ratio * (1-np.exp(-NTU))))  # Kays & London (2-16)
        elif flow_type == 'Cmin_mixed':
            epsilon = 1 - np.exp(-1 / C_ratio * (1 - np.exp(-C * NTU))) # Kays & London (2-15)
        elif flow_type == 'both_mixed':
            epsilon = 1 / (1 / (1 - np.exp(-NTU)) + C_ratio / (1 - np.exp(-C * NTU)) - 1 / NTU) # Kays & London (2-17)

    elif exchanger_type == 'shell_and_tube':

        ntu_p = NTU / shell_passes # Number of heat transfer units in case of shell-and-tube exchanger
        #1 pass epsilon, from Kays & London (2-19) Assumes  mixing  between passes, but Kays & London says this doesn't affect much
        C_ratio_sqrt = np.sqrt(1 + C_ratio**2)
        epsilon_1p = np.where(
            NTU < 10,
            2 / (1 + C_ratio + C_ratio_sqrt * (1 + np.exp(-NTU * C_ratio_sqrt)) / (1 -  np.exp(-NTU * C_ratio_sqrt))),
            2 / (1 + C_ratio + C_ratio_sqrt) )
        if C_ratio > 1 - tol: #Close enough to 1 that will be considered as such
            epsilon = shell_passes * epsilon_1p / (1 + (shell_passes - 1) * epsilon_1p) # Kays & London (2-18a) works for all  shell_pass values
        elif shell_passes == 1:
            epsilon = epsilon_1p
        else: # Kays & London (2-18)
            epsilon = (( (1 - epsilon_1p * C_ratio) / (1 - epsilon_1p) )**shell_passes - 1) / ((( (1 - epsilon_1p * C_ratio) / (1 - epsilon_1p) )**shell_passes) - C_ratio) 
            
    return epsilon


#Test values
C_min = 0.063 *  14500
C_ratio = 0.063 *  14500  /  10 / 1130
exchanger_type = 'shell_and_tube'
shell_passes = 1
ls = ["k--","k-", "k-."]
for i,C  in  enumerate([0,C_ratio,1]):   
    NTU  = np.linspace(0.1,20,100)
    epsilon = epsilon_ntu(NTU, C, exchanger_type, shell_passes=shell_passes)
    #plt.plot(NTU, epsilon, ls[i], label=rf'$C_r$ = {C:.2f}')
    if C==C_ratio:
        eps_max = 0.99*epsilon[-1]
        indices = np.where(epsilon < eps_max)
        if indices[0].size > 0:
            largest_NTU = NTU[indices[-1][-1]]
        else:
            largest_NTU = None
        #print(f"The highest NTU which is 1% lower than eps_max = {epsilon[-1]:.2%} is NTU={largest_NTU:.2f}")


#Rating on one HX (low dp)
epsilon = [0.7,0.75,0.8, 0.85, 0.9, 0.92, 0.94,0.95,0.952]
U_A_req = np.array([18.3*64.9,20.7*67.2,23.9*69,28.6*70.3, 37.1*70.6, 43.3*70.1, 54.9*69.5, 74*69.3, 93.5*68.4])
oversize_fact = np.array([14.59,12.9,11.16,9.3, 7.18, 6.16, 4.86, 3.6, 2.85])
dp = np.array([0.08807,0.08779,0.08747,0.08718,0.0873,0.08767,0.08834,0.0884])/1.5
U_A_max = U_A_req * oversize_fact
'''
print(f"Average pressure drop: {np.mean(dp):.2%} +/- {np.std(dp)*100:.2f} of p_in")

#plt.scatter(U_A_req/C_min, epsilon,label = "Aspen req. area",marker="+")
plt.scatter(U_A_max/C_min, epsilon,label = "Aspen act. area",marker="x")
plt.grid(True)
plt.xlabel(r'Number of Transfer Units  $ N_{TU}  = UA/C_{min} (-)$')
plt.ylabel(r'Effectiveness $\varepsilon =  q/q_{max}  (-)$')
plt.legend()
#plt.scatter([385.7],[0.986])
plt.title('Single pass Shell & Tube HX effectiveness curve')
plt.savefig('plots/eps_NTU_with_Aspen.svg', format='svg')
plt.show() 

'''
