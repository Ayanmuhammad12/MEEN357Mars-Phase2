from define_experiment import experiment1
import numpy as np
from subfunctions import simulate_rover
import matplotlib.pyplot as plt

# ---- Rover + planet dicts (standalone-runnable) ----
motor = {
    'torque_stall':  170,
    'torque_noload':   0,
    'speed_noload':  3.80,
    'mass':          5.0,
    'effcy_tau':  np.array([0, 10, 20, 40, 70, 165]),
    'effcy':      np.array([0, 0.55, 0.75, 0.71, 0.50, 0.05])
}
wheel          = {'radius': 0.30, 'mass': 1.0}
speed_reducer  = {'type': 'reverted', 'diam_pinion': 0.04,
                  'diam_gear': 0.07, 'mass': 1.5}
wheel_assembly = {'wheel': wheel, 'speed_reducer': speed_reducer, 'motor': motor}
rover = {
    'wheel_assembly':  wheel_assembly,
    'chassis':         {'mass': 659},
    'science_payload': {'mass': 75},
    'power_subsys':    {'mass': 90}
}
planet = {'g': 3.72}

experiment, end_event = experiment1()

end_event['max_distance']=1000
end_event['max_time'] = 10000
end_event['min_velocity'] = 0.01

rover = simulate_rover(rover,planet, experiment, end_event)
telem = rover['telemetry']

#plots
fig, axes = plt.subplots(3,1,figsize=(8,10))

axes[0].plot(telem['Time'], telem['position'])
axes[0].set_xlabel('Time [s]')
axes[0].set_ylabel('Position [m]')
axes[0].set_title('Rover Position vs Time')
axes[1].plot(telem['Time'], telem['velocity'])
axes[1].set_xlabel('Time [s]')
axes[1].set_ylabel('Velocity [m/s]')
axes[1].set_title('Rover Velocity vs Time')
axes[2].plot(telem['Time'], telem['power'])
axes[2].set_xlabel('Time [s]')
axes[2].set_ylabel('Power [W]')
axes[2].set_title('Motor Power vs Time')

fig.tight_layout()
fig.savefig('rover_experiment1.png', dpi=150)
plt.close('all')