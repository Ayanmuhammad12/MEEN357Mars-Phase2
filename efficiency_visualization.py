# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 18:28:44 2026

@author: ricodesi88
efficiency_visualization.py
MEEN 357 Design Project,Phase 2, Task 5
Plots DC motor efficiency as a function of motor torque, using a cubic
spline fit through the experimentally measured torque-efficiency data.
"""
#doing our imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


#Define constants as used in Phase 1

motor = {'torque_stall': 170,
         'torque_noload': 0,
         'speed_noload': 3.80,
         'mass': 5.0,
         'effcy_tau': np.array([0,10,20,40,70,165]),
         'effcy': np.array([0,0.55,0.75,0.71,0.50,0.05])}
wheel = {'radius':0.30, 'mass':1.0}
speed_reducer = {'type': 'reverted',
                 'diam_pinion': 0.04,
                 'diam_gear': 0.07,
                 'mass': 1.5}
wheel_assembly= {'wheel': wheel,
                 'speed_reducer':speed_reducer,
                 'motor':motor}
rover={'wheel_assembly': wheel_assembly,
       'chassis':{'mass': 659},
       'science_payload':{'mass':75},
       'power_subsys':{'mass':90}}
planet = {'g': 3.72}

# Pull the torque efficiency data out of the rover dict

effcy_tau = rover['wheel_assembly']['motor']['effcy_tau']
effcy = rover['wheel_assembly']['motor']['effcy']

# Bulding the interpolation function using cubic spine
effcy_fun = interp1d(effcy_tau,effcy,kind = 'cubic')


#evaluating efficiency at 100 evenly spaced torque values.
#Making sure to space them between the minimum and maximum torque given in the rover dictionary.
tau_samples = np.linspace(effcy_tau.min(),effcy_tau.max(),100)

effcy_samples = effcy_fun(tau_samples)

#Plotting
fig, ax = plt.subplots(figsize = (8,5))
ax.plot(effcy_tau, effcy, '*', markersize=12, label = 'Measured data')
ax.plot(tau_samples,effcy_samples,'-',label='Cubic spline fit')
ax.set_xlabel('Motor Shaft Torque [Nm]')
ax.set_ylabel('Motor Efficiency [-]')
ax.set_title('Motor Efficiency vs Torque')
ax.legend()
ax.grid(True)
fig.tight_layout()
fig.savefig('efficiency_visualization.png')
plt.close('all')
