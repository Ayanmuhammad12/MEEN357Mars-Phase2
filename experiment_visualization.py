positions = np.linspace(min(alpha_dist), max(alpha_dist), 100)

angles = alpha_fun(positions)

plt.plot(alpha_dist, alpha_deg, 'r*', label = 'Original Data')

plt.plot(positions, angles, 'b-', label = 'Cubic Spline')

plt.xlabel('Position (m)')
plt.ylabel('Terrain Angles (degree)')
plt.title('Mars Rover Experimental Terrain')
plt.legend()
plt.gird(True)
