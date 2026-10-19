#subfunctions 
###  SIGN CONVENTION: POSITIVE (FORWARD OR DOWNHILL) ; NEGATIVE (BACKWARD OR UPHILL)

#These subfunctions have been imported from phase 1 to include the new subfunctions
###
import numpy as np
from scipy.special import erf
from scipy.interpolate import interp1d

def get_mass(rover):
    #Computes total mass of rover
    if type(rover) != dict:
        raise Exception('get_mass - Input rover must be a dict')
    
    #finding mass
    chassis_mass = rover['chassis']['mass']
    wheel_mass = rover['wheel_assembly']['wheel']['mass']
    motor_mass = rover['wheel_assembly']['motor']['mass']
    payload_mass = rover['science_payload']['mass']
    reducer_mass = rover['wheel_assembly']['speed_reducer']['mass']
    power_mass = rover['power_subsys']['mass']
    m = chassis_mass + payload_mass + power_mass + (6 * (wheel_mass + motor_mass + reducer_mass))

    return m



def get_gear_ratio(speed_reducer):
    #Returns the speed reduction ratio
    if type(speed_reducer) != dict:
            raise Exception('get_gear_reducer - Input speed_reducer must be a dict')

    #checing if 'type' = 'reverted'
    if speed_reducer['type'].lower() != 'reverted':
        raise Exception('get_gear_reducer - type of reducer is not what is expected')

    #reducer
    pinion_diam = speed_reducer['diam_pinion']
    gear_diam = speed_reducer['diam_gear']
    Ng = (gear_diam/pinion_diam)**2
    return Ng


def tau_dcmotor(omega,motor):
    #Returns the motor shaft torque given shaft speed and motor specs
    if type(motor) != dict:
        raise Exception('tau_dcmotor - Input motor must be a dict')
    if not isinstance(omega, (int, float, np.ndarray)):
        raise Exception('tau_dcmotor - Input omega must be a scalar or a numpy array')
    if isinstance(omega, np.ndarray) and omega.ndim != 1:
        raise Exception('tau_dcmotor - omega array must be a 1D numpy array or a scalar')

    tau_stall = motor['torque_stall']
    tau_noload = motor['torque_noload']
    omega_noload = motor['speed_noload']
    slope = ((tau_stall-tau_noload)/omega_noload)

    #torque if scalar
    if isinstance(omega, (int, float)):
         if omega < 0:
              return float(tau_stall)
         elif omega > omega_noload:
              return 0.0
         else:
              return float(tau_stall - (slope * omega))

    #torque for array
    tau = tau_stall - (slope * omega)
    tau[omega < 0] = tau_stall
    tau[omega > omega_noload] = 0

    return tau


def F_drive(omega,rover):
    #Returns the force applied to the rover by the drive system given drive system and shaft speed
    if type(rover) != dict:
        raise Exception('F_drive - Input rover must be a dict')
    if not isinstance(omega, (int, float, np.ndarray)):
        raise Exception('F_drive - Input omega must be a scalar or a numpy array')

    motor = rover['wheel_assembly']['motor']
    speed_reducer = rover['wheel_assembly']['speed_reducer']
    r = rover['wheel_assembly']['wheel']['radius']
    Ng = get_gear_ratio(speed_reducer)
    tau = tau_dcmotor(omega,motor)

    Fd = 6 * (Ng * tau)/r

    return Fd
 

def F_gravity(terrain_angle, rover, planet):
    '''magnitude of the force component acting on the rover in the direction of its
translational motion due to gravity as a function of terrain inclination angle and rover
properties'''
    if not isinstance(terrain_angle, (int, float, np.ndarray)):
        raise Exception('F_gravity - Input terrain angle must be a scalar or a numpy array')
    if np.any(np.array(terrain_angle) > 75) or  np.any(np.array(terrain_angle) < -75):
        raise Exception('F_gravity - Terrain angle must be between -75 and 75 degrees')
    if type(rover) != dict or type(planet) != dict:
            raise Exception('F_gravity - Input rover and planet must both be dict')
    m = get_mass(rover)
    g = planet['g']

    radang = np.deg2rad(terrain_angle)
    Fgt = -m * g * np.sin(radang)
    
    return Fgt

def F_rolling(omega, terrain_angle, rover, planet, Crr):
    #magnitude of force component due to rolling resistances given the terrain inclination angle, rover properties, and a
#rolling resistance coefficient
    if not isinstance(omega, (float, int, np.ndarray)):
        raise Exception('F_rolling - Omega must be a scalar or numpy array')
    if not isinstance(terrain_angle,(float, int, np.ndarray)):
         raise Exception('F_rolling - terrain angle must be a scalar or numpy array')
    if np.shape(np.array(omega)) != np.shape(np.array(terrain_angle)):
         raise Exception('F_rolling - omega and terrain angle must have the same shape')
    if np.any(np.array(terrain_angle) > 75) or  np.any(np.array(terrain_angle) < -75):
        raise Exception('F_rolling - terrain angle must be between -75 and 75 degrees')
    if type(rover) != dict:
        raise Exception('F_rolling - Input rover must be a dict')
    if type(planet) != dict:
            raise Exception('F-rolling - Input planet must be a dict')
    if not (isinstance(Crr, (int, float)) and Crr > 0):
        raise Exception('F_rolling - Crr must be a positive scalar')
    
    
    m = get_mass(rover)
    g = planet['g']
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    r = rover['wheel_assembly']['wheel']['radius']

    radang = np.deg2rad(terrain_angle)

    v_rover = (omega/Ng) * r

    Frr_simple = m * g * np.cos(radang) * Crr
    Frr = -erf(40*v_rover) * Frr_simple
    
    return Frr

def F_net(omega, terrain_angle, rover, planet, Crr):
    # magnitude of net force acting on the rover
    if not isinstance(omega, (float, int, np.ndarray)):
        raise Exception('F_net - Omega must be a scalar or numpy array')
    if not isinstance(terrain_angle,(float, int, np.ndarray)):
        raise Exception('F_net - terrain angle must be a scalar or numpy array')
    if np.shape(np.array(omega)) != np.shape(np.array(terrain_angle)):
        raise Exception('F_net - omega and terrain angle must have the same shape')
    if np.any(np.array(terrain_angle) > 75) or  np.any(np.array(terrain_angle) < -75):
        raise Exception('F_net - terrain angle must be between -75 and 75 degrees')
    if type(rover) != dict:
        raise Exception('F_net - Input rover must be a dict')
    if type(planet) != dict:
        raise Exception('F_net - Input planet must be a dict')
    if not (isinstance(Crr, (int, float)) and Crr > 0):
        raise Exception('F_net - Crr must be a positive scalar')

    Fd = F_drive(omega,rover)
    Fgt = F_gravity(terrain_angle, rover, planet)
    Frr = F_rolling(omega, terrain_angle, rover, planet, Crr)

    Fnet = Fd + Fgt + Frr

    return Fnet


"""
Everything below this line is new for phase 2 of the project
"""


def motorW(v, rover):
    # finds rotational speed of the motor shaft in rad/s
    if type(rover) != dict:
        raise Exception('motorW input "rover" must be a dictionary')
        
    if not isinstance(v, (int, float, np.ndarray)):
        raise Exception('motorW velocity must be a scalar or numpy array')
        
    if isinstance(v, np.ndarray) and v.ndim != 1:
        raise Exception('motorW velocity must be a 1D array')
        
    r = rover['wheel_assembly']['wheel']['radius']
    
    speed_reducer = rover['wheel_assembly']['speed_reducer']
    
    Ng = get_gear_ratio(speed_reducer)
    
    w = (v /r) * Ng
    
    return w

def rover_dynamics(t, y, rover, planet, experiment):
    if not isinstance(t, (int, float)):
        raise Exception('rover_dynamics - Input time, t, should be a scalar')
    if not isinstance(y, np.ndarray) or y.shape!= (2,):
        raise Exception('rover_dynamics - Input state vector, y, should be a 1D array containing velocity and position')
    if type(rover) != dict:
        raise Exception('rover_dynamics - Input rover must be a dict')
    if type(planet) != dict:
            raise Exception('rover_dynamics - Input planet must be a dict')
    if type(experiment) != dict:
                raise Exception('rover_dynamics - Input experiment must be a dict')
    # def our state vector
    v = float(y[0])
    x = float(y[1])

    #interpolation function for terrain_angle
    alpha_fun = interp1d(experiment['alpha_dist'], experiment['alpha_deg'], kind = 'cubic', fill_value = 'extrapolate')
    terrain_angle = float(alpha_fun(x))

    #set up equation
    omega = motorW(v,rover)
    Fnet = F_net(omega, terrain_angle, rover, planet, experiment['Crr'])

    m = get_mass(rover)
    a = Fnet/m
    dydt = np.array([a,v])

    return dydt



def mechpower(v,rover):
    if not isinstance(v, (int, np.ndarray, float )):
        raise Exception('mechpower - Input v must be a scalar or a numpy array if vector')
    if isinstance(v, np.ndarray) and v.ndim != 1:
        raise Exception('mechpower - Input v must be a scalar or a 1D numpy array, not a matrix')
    if type(rover) != dict:
        raise Exception('mechpower - Input rover must be a dict')

    # P = tau * angular velo
    motor = rover['wheel_assembly']['motor']
    w = motorW(v,rover)
    tau = tau_dcmotor(w,motor)
    P = tau * w

    return P

def battenery(t, v, rover):
    """
    Computes the total electrical energy consumed from the rover battery pack
    over a simulation profile defined by time-velocity pairs, accounting for all 
    6 motors and for DC motor inefficiency.
    
    """


#input validation

    if not isinstance(t, np.ndarray) or t.ndim !=1:
        raise Exception('battenergy - Input t must be a 1D numpy array')
    if not isinstance(v, np.ndarray) or v.ndim !=1:
        raise Exception('battenergy-Input v must be a 1D numpy array')
    if len(t) !=len(v):
        raise Exception ('battenergy- t and v must be the same length')
    if type(rover)!=dict:
        raise Exception('battenergy - Input rover must be a dict')
    

    motor = rover['wheel_assembly']['motor']
    # omega is the motor speed, tau is the motor torque at each time sample [W]
    omega = motorW(v, rover)
    tau = tau_dcmotor(omega,motor)
    #mechanical power output of a single motor at each time sample [W]
    P_motor = mechpower(v,rover)
    #interpolation function for motor efficiency
    effcy_tau = motor['effcy_tau']
    effcy = motor['effcy']
    effcy_fun = interp1d(effcy_tau, effcy, kind = 'cubic')
    #clip torque into range of efficiency data before interpolating
    #because 170Nm is greater than 165Nm which is max bounds at zero rover velocity
    tau_clipped = np.clip(tau,np.min(effcy_tau),np.max(effcy_tau))

    #motor efficiency at each operating point [-]
    eta = effcy_fun(tau_clipped)
    #Samples set to nonzero to avoid dividing by zero

    P_batt_single = np.zeros(len(v), dtype=float)
    nonzero = eta >0
    P_batt_single[nonzero]=P_motor[nonzero] / eta[nonzero]

#   Multiply by number of wheels (6)
    P_batt_total = 6*P_batt_single

    # integrate battery power to get the total energy in joules
    E = np.trapezoid(P_batt_total, t)
           

    return E
