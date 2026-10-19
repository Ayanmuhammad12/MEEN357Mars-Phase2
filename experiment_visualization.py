#visualizing the terrain by running interp1d 

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from define_experiment import experiment1

experiment, end_event = experiment1()
alpha_dist = experiment['alpha_dist']
alpha_deg  = experiment['alpha_deg']
alpha_fun = interp1d(experiment['alpha_dist'], experiment['alpha_deg'], kind = 'cubic', fill_value = 'extrapolate')

positions = np.linspace(min(alpha_dist), max(alpha_dist), 100)

angles = alpha_fun(positions)

plt.plot(alpha_dist, alpha_deg, 'r*', label = 'Original Data')

plt.plot(positions, angles, 'b-', label = 'Cubic Spline')

plt.xlabel('Position (m)')
plt.ylabel('Terrain Angles (degree)')
plt.title('Mars Rover Experimental Terrain')
plt.legend()
plt.grid(True)
plt.savefig('experiment_visualization.png', dpi=150)
plt.close('all')