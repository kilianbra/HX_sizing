import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import PercentFormatter

# Data
labels = ['Model'] #'Aspen',

# Aspen data
mass_aspen = 1
shell_htc_aspen = 1
tube_htc_aspen = 1
pressure_drop_aspen = 1

# Data values
shell_numb = 4
baseline_mass = 6.175
baseline_htc_o = 230.5
baseline_htc_i_outter_area =  99
baseline_dp_p_in_air_shell_side = 0.0869/1.5

mass_data = 4.9/baseline_mass

shell_htc_data = 192/baseline_htc_o #149 before modification of A_o_cr
shell_htc_mean = shell_htc_data/0.85
shell_htc_min = shell_htc_data/3
shell_htc_max = shell_htc_data/0.5

tube_htc = 54.3/baseline_htc_i_outter_area
#For pressure drop, use inlet and outlet density for I/O, and then add xflow and window at average
shell_dp_data = 0.201/(1.2)**3/(baseline_dp_p_in_air_shell_side * 1.5)
shell_dp_error_adjusted = shell_dp_data/0.95
shell_dp_min = shell_dp_data/2
shell_dp_max = shell_dp_data/0.95


# Bar positions
x = np.arange(len(labels))
width = 0.1

fig, ax = plt.subplots()

# Plotting bars $h_o (W/m^2/K)$ \Delta p/p_{in} 
rects1 = ax.bar(x - 3*width/2, [mass_data-1], width, label='Mass (t)', color='black')
rects2 = ax.bar(x - width/2, [shell_htc_data-1], width, label='Shell side htc $(W/m^2/K)$', color='red')
rects3 = ax.bar(x + width/2, [ tube_htc-1], width, label='Tube side htc $(W/m^2/K)$', color='blue')
rects4 = ax.bar(x + 3*width/2, [shell_dp_data-1], width, label="Shell side $\Delta p/p_{in} (\%)$ Aspen $V_{avg}$", color='green')
rects4prime = ax.bar(x + 5/2 * width, [ shell_dp_data*1.2**3-1], width,label="Shell side $\Delta p/p_{in} (\%)$ Model $A_c$",color='green', alpha=0.5)



# Adding text labels
#ax.text(x[0] - 3*width/2, mass_aspen/2, rf"{baseline_mass*shell_numb:.1f}", ha='center', va='center', color='white')
#ax.text(x[0] - width/2, shell_htc_aspen/2, f"{baseline_htc_o:.0f}", ha='center', va='center', color='white')
#ax.text(x[0] + width/2, tube_htc_aspen/2, f"{baseline_htc_i_outter_area:.0f}", ha='center', va='center', color='white')
#ax.text(x[0] + 3*width/2, pressure_drop_aspen/2, f"{baseline_dp_p_in_air_shell_side:.1%}", ha='center', va='center', color='white')

# Setting labels, title, and custom x-axis tick labels
#ax.set_xlabel('Agreement with Aspen')
ax.set_ylabel(r'Relative difference to Aspen $\frac{x-x_{Aspen}}{x_{Aspen}}$')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend(loc="upper left")

#ax.errorbar([1-width/2],[shell_htc_mean],ms=20,color="k",marker="+",yerr= [[shell_htc_data - shell_htc_min], [shell_htc_max - shell_htc_data]], capsize=10)
#ax.errorbar([1+3*width/2],[shell_dp_data],ms=20,color="k",marker="|",yerr= [[shell_dp_data - shell_dp_min], [shell_dp_max - shell_dp_data]], capsize=10)
#ax.scatter([1+3*width/2], [shell_dp_error_adjusted],s=250,color="k",marker="+",lw=1)
#ax.errorbar([1+3*width/2], [shell_dp_error_adjusted],ms=20,color="k",marker="+",yerr= [[shell_dp_error_adjusted - shell_dp_ro_in], [-shell_dp_error_adjusted + shell_dp_ro_out_aspen]], capsize=10,ecolor="red")
#ax.errorbar([1+3*width/2], [shell_dp_error_adjusted],ms=20,color="k",marker="+",yerr= [[shell_dp_error_adjusted - shell_dp_ro_in/0.95], [-shell_dp_error_adjusted + shell_dp_ro_out_aspen/0.95]], capsize=10)
plt.ylim(-1,1.5)
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
fig.tight_layout()
plt.savefig('plots/aspen_v_model_lowdp.svg', format='svg')
plt.show()