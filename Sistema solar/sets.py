from matplotlib import pyplot as plt
from matplotlib import animation
from math import cos, sin, sqrt, pi
import json
from scipy.integrate import ode
from scipy.optimize import fmin_slsqp
H_TO_S = 3600
S_TO_H = 1 / 3600
M_TO_KM = 1 / 1000
KM_TO_M = 1000
G_TO_KG = 1 / 1000
KG_TO_G = 1000
CM_TO_KM = 1 / 100000
KM_TO_CM = 100000
MICRO_TO_UNI = 1E-6
UNI_TO_MICRO = 1E6
DISTANCE_SUN_EARTH = 149600000 # Km
DISTANCE_SUN_MARS = 227940000 # Km
DISTANCE_SUN_JUPITER = 7786000000 # Km
EARTH_ORBITAL_VELOCITY = 29.8 / S_TO_H # Km/h
MARS_ORBITAL_VELOCITY = 24.1 / S_TO_H # Km/h
EARTH_ANGULAR_VELOCITY = EARTH_ORBITAL_VELOCITY / DISTANCE_SUN_EARTH # rad/h
MARS_ANGULAR_VELOCITY = MARS_ORBITAL_VELOCITY / DISTANCE_SUN_MARS # rad/h
SUN_RADIUS = 695700 # Km
SUN_MASS = 1988500E4 # 10^20Kg
GRAVITATIONAL_CONSTANT = 6.67408 * 1E9 * M_TO_KM ** 3 / S_TO_H ** 2 # 10^(-20) Km^3/(kg*h^2)
GM = GRAVITATIONAL_CONSTANT * SUN_MASS # Km^3/h^2
MAX_SPEED_DIFFER = 9 / S_TO_H # Km/h
MARS_RADIUS = 6792 / 2 # Km
MARS_EARTH_INITIAL_ANGLE_DIFF = 0 # rad
EARTH_ESCAPE_VELOCITY = 11.186 / S_TO_H # Km/h
ODE_SOLVER_STEP = 1 / 500
MAX_WEIGHT = 2E3 # Kg maximum weight we can carry with our sail
SAIL_DENSITY = 7 * G_TO_KG / CM_TO_KM ** 2 # Kg/Km^2
P0 = 4.563 * MICRO_TO_UNI / (M_TO_KM * S_TO_H ** 2) # (Kg/Km*h^2)
r02 = DISTANCE_SUN_EARTH ** 2 # Km
ac = 2 * P0 * r02 / MAX_WEIGHT # Km/h^2
b1 = 0.0864
b2 = 0.8277
b3 = -5.45e-3
N = 15 # Default number of partitions
AREA_RELEVANCE = 1. # default relevance given to area minimization
TIME_RELEVANCE = 1. # default relevance given to time minimization
R_init = DISTANCE_SUN_EARTH
R_fin = DISTANCE_SUN_MARS
# define a function which given a value s between 0 and 1 returns the angle that the sail should
be
# form with the position vector at the corresponding time.
def alpha(s, alpha_list):
       s = min(s, 0.99)
return alpha_list[int(s * len(alpha_list))]
# Define the ordinary differential equation we must solve in order to compute the final position
# and velocity of our solar wind propelled rocket first parameter is radius (from Sun), second is
# angle (respect to initial angle) third is radial speed and last is angular speed. The parameters
# consist on a list of N+3 doubles. First N doubles correspond to the angle the sail is forming with
# the position vector taking the Sun as coordinates origin. The fourth to last parameter is a
value
# between 0 and 1 which indicates the initial speed the sail has whith respect to the Earth
# (in proortion to escape velocity). The third to last parameter is the initial
# angle that the ship forms with Earth’s orbit. The second to last parameter is the sail area
# and the last parameter is the expected time the trip is going to last.
def f(t, y, params):
       r, theta, u, v = y
A = params[-2]
tf = params[-1]
r_deriv = tf * u
theta_deriv = tf * v / r
u_deriv = tf * (v ** 2 - GM / r) / r
v_deriv = -tf * u * v / r
a = alpha(t, params[:-4])
cosalpha = cos(a)
sinalpha = sin(a)
aux = tf * ac * A * cosalpha / r ** 2
u_deriv += aux * (b1 + cosalpha * (b2 * cosalpha + b3))
v_deriv += aux * sinalpha * (b2 * cosalpha + b3)
return [r_deriv, theta_deriv, u_deriv, v_deriv]
# define function that solves the differential equation and returns the final values of position and speed in
# polar coordinates.
def solve_ode(params):
# Initialize an object to solve the previous differential equation using Runge-Kutta 4(5) method.ode_solver = ode(f).set_integrator(’dopri5’, nsteps=10000)
# Set initial value
        theta_init = params[-3]
v_init = params[-4]
ode_solver.set_initial_value([DISTANCE_SUN_EARTH, 0, sin(theta_init) * v_init,
cos(theta_init) * v_init + EARTH_ORBITAL_VELOCITY],
0)
ode_solver.set_f_params(params)
while ode_solver.successful() and ode_solver.t < 1:
            ode_solver.integrate(ode_solver.t + ODE_SOLVER_STEP)
            return ode_solver.y
# define function that solves the differential equation and returns the final values of position and speed in
# polar coordinates.
def solve_ode_orbital(params):
            global R_init
# Initialize an object to solve the previous differential equation using Runge-Kutta 4(5) method.
ode_solver = ode(f).set_integrator('dopri5', nsteps=10000)
# Set initial value
theta_init = params[-3]
v_init = params[-4]
ode_solver.set_initial_value([R_init, 0, sin(theta_init) * v_init,
cos(theta_init) * v_init + sqrt(GM/R_init)],
0)
ode_solver.set_f_params(params)
while ode_solver.successful() and ode_solver.t < 1:
   ode_solver.integrate(ode_solver.t + ODE_SOLVER_STEP)
return ode_solver.y
# Our goal is that of maximizing the amount of material we can receive and return per time unit.
def objective(params):
    global AREA_RELEVANCE, TIME_RELEVANCE
r_fin, theta_fin, radial_speed, theta_speed = solve_ode(params)
tang_speed_diff = (theta_speed - DISTANCE_SUN_MARS * MARS_ANGULAR_VELOCITY)
speed_condition = sqrt(max(tang_speed_diff ** 2 + radial_speed ** 2 - MAX_SPEED_DIFFER ** 2,
0))
Mars_theta = MARS_EARTH_INITIAL_ANGLE_DIFF + MARS_ANGULAR_VELOCITY * params[-1]
x_diff = r_fin * cos(theta_fin) - DISTANCE_SUN_MARS * cos(Mars_theta)
y_diff = r_fin * sin(theta_fin) - DISTANCE_SUN_MARS * sin(Mars_theta)
pos_condition = sqrt(max(x_diff ** 2 + y_diff ** 2 - MARS_RADIUS ** 2, 0))
# we square the value to make sure it is not negative.
# this option results in error, the NLP solver learns to make the value params[-1]
#infinitely small.
# return -((MAX_WEIGHT - SAIL_DENSITY * params[-2]) / params[-1]) ** 2
#we add constraints (terms multiplied by 1e8) in order to achieve the desired results that
#being a small
# final distance to Mars (of the order of Mars radius) and a small relative speed to Mars
(smaller that 9 Km/s)
time_condition = TIME_RELEVANCE * abs(params[-1])
area_condition = AREA_RELEVANCE * params[-2] ** 2
return pos_condition * 1E-4 + speed_condition * 5E-3 + time_condition / 180 + area_condition
* 5
# Our goal is that of maximizing the amount of material we can receive and return per time unit.
def orbit_objective(params):
global R_fin
r_fin, theta_fin, radial_speed, theta_speed = solve_ode_orbital(params)
tang_speed_diff = abs(theta_speed - sqrt(GM/R_fin))
speed_condition = sqrt(max(tang_speed_diff**2 + radial_speed ** 2 - MAX_SPEED_DIFFER ** 2,
0))
orbit_condition = abs(r_fin - R_fin)
# we square the value to make sure it is not negative.
# we add constraints (terms multiplied by 1e8) in order to achieve the desired results that
being a small
# final distance to Mars (of the order of Mars radius) and a small relative speed to Mars
(smaller that 9 Km/s)
time_condition = TIME_RELEVANCE * abs(params[-1])
area_condition = AREA_RELEVANCE * params[-2] ** 2
return orbit_condition * 1E-4 + speed_condition * 5E-3 + time_condition / 180 +
area_condition * 5
def find_optimal_values(n=N, verbose=True, area_relevance=AREA_RELEVANCE,
time_relevance=TIME_RELEVANCE):
global AREA_RELEVANCE, TIME_RELEVANCE
AREA_RELEVANCE = area_relevance
TIME_RELEVANCE = time_relevance
# Define constraints for the angles values (angles must be in the range [-pi/2, pi/2])
bounds = [(-pi / 2, pi / 2) for i in range(n)]
# Define bounds for initial speed respect to Earth (from 0 to escape velocity).
bounds.append((0, EARTH_ESCAPE_VELOCITY))
# Define bounds for the initil angle formed with Earth.
bounds.append((-pi, pi))
# Define constraints for the Area of the sail (area must be between 0 and maximum
transportable weight.
bounds.append((0, 0.28))
# Define time constraints between 10 hours and 1.8 years
bounds.append((10, 24 * 1.8 * 365))
initial = [11 * pi / 180 for i in range(n)] # angles
initial.append(EARTH_ESCAPE_VELOCITY / 100) # Initial speed relative to Earth
initial.append(0) # initial angle with Earth orbit
initial.append(0.14) # Area sail.
initial.append(500 * 24) # Time.
if verbose:
results = fmin_slsqp(objective, initial, bounds=bounds, iprint=2, iter=1000)
else:
results = fmin_slsqp(objective, initial, bounds=bounds, iter=1000)
show_results(results)
return results
def show_results(params):
angles = params[:-4]
speed_init, angle_init, sail_size, time = params[-4:]
print("\nThe optimal angle values are:")
for angle in angles:
print("\t{:.2f} radians".format(angle))
print("\nThe optimal initial speed respect to Earth is:")
print("\t{:.2f} Km/s".format(speed_init / H_TO_S))
print("\nThe optimal initial angle formed with the Earth’s orbit is:")
print("\t{:.2f} radians".format(angle_init))
print("\nThe optimal sail size is:")
print("\t{:.2f} Km^2".format(sail_size))
print("\nThe optimal travel time is:")
print("\t{:.2f} days".format(time / 24))
r_fin, theta_fin, radial_speed, tangential_speed = solve_ode(params)
print("\nSpeed difference between solar sail and Mars is:")
tang_speed_diff = tangential_speed - DISTANCE_SUN_MARS * MARS_ANGULAR_VELOCITY
print("\t{:.2f} Km/s".format(sqrt(tang_speed_diff ** 2 + radial_speed ** 2) / H_TO_S))
print("\nDistance between the center of Mars (with radius 3396 Km) and the solar sail at
trajectory end is:")
theta_mars = MARS_EARTH_INITIAL_ANGLE_DIFF + MARS_ANGULAR_VELOCITY * params[-1]
x_diff = (r_fin * cos(theta_fin) - DISTANCE_SUN_MARS * cos(theta_mars))
y_diff = (r_fin * sin(theta_fin) - DISTANCE_SUN_MARS * sin(theta_mars))
print("\t{:.2f} Km".format(sqrt(x_diff ** 2 + y_diff ** 2)))
def show_orbital_results(params):
angles = params[:-4]
speed_init, angle_init, sail_size, time = params[-4:]
print("\nThe optimal angle values are:")
for angle in angles:
print("\t{:.2f} radians".format(angle))
print("\nThe optimal initial speed respect to Earth is:")
print("\t{:.2f} Km/s".format(speed_init / H_TO_S))
print("\nThe optimal initial angle formed with the Earth’s orbit is:")
print("\t{:.2f} radians".format(angle_init))
print("\nThe optimal sail size is:")
print("\t{:.2f} Km^2".format(sail_size))
print("\nThe optimal travel time is:")
print("\t{:.2f} days".format(time / 24))
r_fin, theta_fin, radial_speed, tangential_speed = solve_ode_orbital(params)
print("\nSpeed difference between solar sail Orbital speed and Mars is:")
tang_speed_diff = tangential_speed - sqrt(GM/R_fin)
print("\t{:.2f} Km/s".format(sqrt(tang_speed_diff ** 2 + radial_speed ** 2) / H_TO_S))
print("\nDistance between the solar sail and the orbit:")
print("\t{:.2f} Km".format(abs(r_fin-R_fin)))
print("\nThe Mars angle respect Earth should be:")
print("\t{:.2f} radians".format(theta_fin - MARS_ANGULAR_VELOCITY * params[-1]))
# Define function that solves the differential equation and .
def generate_plotting_data(params):
# Initialize an object to solve the previous differential equation using Runge-Kutta 4(5)
method.
ode_solver = ode(f).set_integrator(’dopri5’, nsteps=10000)
# Set initial value
theta_init = params[-3]
v_init = params[-4]
ode_solver.set_initial_value([DISTANCE_SUN_EARTH, 0, sin(theta_init) * v_init,
                              cos(theta_init) * v_init + EARTH_ORBITAL_VELOCITY],
0)
ode_solver.set_f_params(params)
positions = []
while ode_solver.successful() and ode_solver.t + ODE_SOLVER_STEP < 1:
ode_solver.integrate(ode_solver.t + ODE_SOLVER_STEP)
positions.append([ode_solver.y[0] * cos(ode_solver.y[1]), ode_solver.y[0] *
sin(ode_solver.y[1])])
return positions
# Define function that solves the differential equation and .
def generate_orbital_plotting_data(params):
# Initialize an object to solve the previous differential equation using Runge-Kutta 4(5)
method.
ode_solver = ode(f).set_integrator(’dopri5’, nsteps=10000)
# Set initial value
theta_init = params[-3]
v_init = params[-4]
ode_solver.set_initial_value([R_init, 0, sin(theta_init) * v_init,
cos(theta_init) * v_init + sqrt(GM/R_init)],
0)
ode_solver.set_f_params(params)
positions = []
while ode_solver.successful() and ode_solver.t + ODE_SOLVER_STEP < 1:
ode_solver.integrate(ode_solver.t + ODE_SOLVER_STEP)
positions.append([ode_solver.y[0] * cos(ode_solver.y[1]), ode_solver.y[0] *
sin(ode_solver.y[1])])
return positions
# Creating celestial_object class
class CelestialObject(object):
def __init__(self, x=0, y=0):
self.x = x
self.y = y
def move(self, newx, newy):
self.x = newx
self.y = newy
def save_params(params, file=’learned_parameters.json’):
angles = [angle for angle in params[:-4]]
results = {’angles’: angles,
’initial_speed’: params[-4],
’initial_angle’: params[-3],
’sail_size’: params[-2],
’time’: params[-1]}
with open(file, ’w’) as outfile:
json.dump(results, outfile)
def load_params(file=’learned_parameters.json’):
json_data = open("learned_parameters.json", "r")
results = json.load(json_data)
params = [results[’initial_speed’], results[’initial_angle’], results[’sail_size’],
results[’time’]]
params = results[’angles’] + params
return params
def show_animated_plot(params):
positions = generate_plotting_data(params)
tf = params[-1]
n = len(positions)
step = tf / n
# Initializing dots
Earth = CelestialObject(DISTANCE_SUN_EARTH, 0)
Ship = CelestialObject(DISTANCE_SUN_EARTH, 0)
Mars = CelestialObject(DISTANCE_SUN_MARS, 0)
plt.clf()
fig = plt.figure()
ax = plt.axes(xlim=(-2 * DISTANCE_SUN_MARS, 2 * DISTANCE_SUN_MARS),
ylim=(-2 * DISTANCE_SUN_MARS, 2 * DISTANCE_SUN_MARS))
d, = ax.plot([Earth.x, Mars.x, Ship.x],
[Earth.y, Mars.y, Ship.y], ’ro’, markersize=4)
circle = plt.Circle((0, 0), DISTANCE_SUN_MARS, color=’r’, fill=False)
ax.add_artist(circle)
circle = plt.Circle((0, 0), DISTANCE_SUN_EARTH, color=’b’, fill=False)
ax.add_artist(circle)
circle = plt.Circle((0, 0), 10 * SUN_RADIUS, color=’y’)
ax.add_artist(circle)
# animation function. This is called sequentially
def animate(i):
Earth_theta = i * step * EARTH_ANGULAR_VELOCITY
Earth.move(DISTANCE_SUN_EARTH * cos(Earth_theta), DISTANCE_SUN_EARTH * sin(Earth_theta))
Mars_theta = MARS_EARTH_INITIAL_ANGLE_DIFF + i * step * MARS_ANGULAR_VELOCITY
Mars.move(DISTANCE_SUN_MARS * cos(Mars_theta), DISTANCE_SUN_MARS * sin(Mars_theta))
x, y = positions[i]
Ship.move(x, y)
d.set_data([Earth.x, Mars.x, Ship.x],
[Earth.y, Mars.y, Ship.y])
plt.title("Time: {:.2f} days".format(i * step / 24))
return d,
# call the animator.
anim = animation.FuncAnimation(fig, animate, frames=n, interval=5)
plt.show()
def show_orbital_animated_plot(params):
positions = generate_orbital_plotting_data(params)
tf = params[-1]
n = len(positions)
step = tf / n
# Initializing dots
Earth = CelestialObject(DISTANCE_SUN_EARTH, 0)
Ship = CelestialObject(DISTANCE_SUN_EARTH, 0)
Mars = CelestialObject(DISTANCE_SUN_MARS, 0)
plt.clf()
fig = plt.figure()
ax = plt.axes(xlim=(-2 * DISTANCE_SUN_MARS, 2 * DISTANCE_SUN_MARS),
ylim=(-2 * DISTANCE_SUN_MARS, 2 * DISTANCE_SUN_MARS))
d, = ax.plot([Earth.x, Mars.x, Ship.x],
[Earth.y, Mars.y, Ship.y], ’ro’, markersize=4)
circle = plt.Circle((0, 0), DISTANCE_SUN_MARS, color=’r’, fill=False)
ax.add_artist(circle)
circle = plt.Circle((0, 0), DISTANCE_SUN_EARTH, color=’b’, fill=False)
ax.add_artist(circle)
circle = plt.Circle((0, 0), 10 * SUN_RADIUS, color=’y’)
ax.add_artist(circle)
# animation function. This is called sequentially
def animate(i):
Earth_theta = i * step * EARTH_ANGULAR_VELOCITY
Earth.move(DISTANCE_SUN_EARTH * cos(Earth_theta), DISTANCE_SUN_EARTH * sin(Earth_theta))
Mars_theta = MARS_EARTH_INITIAL_ANGLE_DIFF + i * step * MARS_ANGULAR_VELOCITY
Mars.move(DISTANCE_SUN_MARS * cos(Mars_theta), DISTANCE_SUN_MARS * sin(Mars_theta))
x, y = positions[i]
Ship.move(x, y)
d.set_data([Earth.x, Mars.x, Ship.x],
[Earth.y, Mars.y, Ship.y])
plt.title("Time: {:.2f} days".format(i * step / 24))
return d,
# call the animator.
anim = animation.FuncAnimation(fig, animate, frames=n, interval=5)
plt.show()
def save_animated_plot(params, file_name=’animated_plot.html’):
Writer = animation.writers[’ffmpeg’]
writer = Writer(fps=15, metadata=dict(artist=’Team 744’), bitrate=1800)
positions = generate_plotting_data(params)
tf = params[-1]
n = len(positions)
step = tf / n
# Initializing dots
Earth = CelestialObject(DISTANCE_SUN_EARTH, 0)
Ship = CelestialObject(DISTANCE_SUN_EARTH, 0)
Mars = CelestialObject(DISTANCE_SUN_MARS, 0)
fig = plt.figure()
ax = plt.axes(xlim=(-2 * DISTANCE_SUN_MARS, 2 * DISTANCE_SUN_MARS),
ylim=(-2 * DISTANCE_SUN_MARS, 2 * DISTANCE_SUN_MARS))
d, = ax.plot([Earth.x, Mars.x, Ship.x],
[Earth.y, Mars.y, Ship.y], ’ro’, markersize=4)
circle = plt.Circle((0, 0), DISTANCE_SUN_MARS, color=’r’, fill=False)
ax.add_artist(circle)
circle = plt.Circle((0, 0), DISTANCE_SUN_EARTH, color=’b’, fill=False)
ax.add_artist(circle)
circle = plt.Circle((0, 0), 10 * SUN_RADIUS, color=’y’)
ax.add_artist(circle)
# animation function. This is called sequentially
def animate(i):
Earth_theta = i * step * EARTH_ANGULAR_VELOCITY
Earth.move(DISTANCE_SUN_EARTH * cos(Earth_theta), DISTANCE_SUN_EARTH * sin(Earth_theta))
Mars_theta = i * step * MARS_ANGULAR_VELOCITY
Mars.move(DISTANCE_SUN_MARS * cos(Mars_theta), DISTANCE_SUN_MARS * sin(Mars_theta))
x, y = positions[i]
Ship.move(x, y)
d.set_data([Earth.x, Mars.x, Ship.x],
[Earth.y, Mars.y, Ship.y])
plt.title("Time: {:.2f} days".format(i * step / 24))
return d,
# call the animator.
anim = animation.FuncAnimation(fig, animate, frames=n, interval=5)
anim.save(filename=file_name, writer=writer)
def plot_trajectory(params):
positions = generate_plotting_data(params)
x = [pos[0] for pos in positions]
y = [pos[1] for pos in positions]
plt.clf()
plt.plot(x, y, ’r’, color=’black’)
N = 100
x = [DISTANCE_SUN_EARTH * cos(2 * pi * i / N) for i in range(N + 1)]
y = [DISTANCE_SUN_EARTH * sin(2 * pi * i / N) for i in range(N + 1)]
plt.plot(x, y, ’r’, color=’blue’)
x = [DISTANCE_SUN_MARS * cos(2 * pi * i / N) for i in range(N + 1)]
y = [DISTANCE_SUN_MARS * sin(2 * pi * i / N) for i in range(N + 1)]
plt.plot(x, y, ’r’, color=’red’)
plt.title(’{} angle variations’.format(len(params) - 4))
show_results(params)
plt.show()
def save_trajectory(params, file_name=’animated_plot.mp4’):
positions = generate_plotting_data(params)
x = [pos[0] for pos in positions]
y = [pos[1] for pos in positions]
plt.clf()
plt.plot(x, y, ’r’, color=’black’)
N = 100
x = [DISTANCE_SUN_EARTH * cos(2 * pi * i / N) for i in range(N + 1)]
y = [DISTANCE_SUN_EARTH * sin(2 * pi * i / N) for i in range(N + 1)]
plt.plot(x, y, ’r’, color=’blue’)
x = [DISTANCE_SUN_MARS * cos(2 * pi * i / N) for i in range(N + 1)]
y = [DISTANCE_SUN_MARS * sin(2 * pi * i / N) for i in range(N + 1)]
plt.plot(x, y, ’r’, color=’red’)
plt.title(’{} angle variations’.format(len(params) - 4))
plt.savefig(file_name)
def save_orbital_trajectory(params, file_name=’animated_plot.mp4’):
positions = generate_orbital_plotting_data(params)
x = [pos[0] for pos in positions]
y = [pos[1] for pos in positions]
plt.clf()
plt.plot(x, y, ’r’, color=’black’)
N = 100
x = [R_init * cos(2 * pi * i / N) for i in range(N + 1)]
y = [R_init * sin(2 * pi * i / N) for i in range(N + 1)]
plt.plot(x, y, ’r’, color=’blue’)
x = [R_fin * cos(2 * pi * i / N) for i in range(N + 1)]
y = [R_fin * sin(2 * pi * i / N) for i in range(N + 1)]
plt.plot(x, y, ’r’, color=’red’)
plt.title(’{} angle variations’.format(len(params) - 4))
plt.savefig(file_name)
def save_areatime(params, file=’learned_parameters.json’):
angles = [angle for angle in params[:-4]]
results = {’angles’: angles,
’initial_speed’: params[-4],
’initial_angle’: params[-3],
’sail_size’: params[-2],
’time’: params[-1]}
with open(file, ’w’) as outfile:
json.dump(results, outfile)
def plot_angles(params):
angles = params[:-4]
n = len(angles)
t = params[-1] / n
x = [(i + 1 / 2) * t / 24 for i in range(n)]
plt.clf()
plt.xlabel(’time (days)’)
plt.ylabel(’Angle (radians)’)
plt.plot(x, angles, ’r’)
plt.show()
def plot_multiple_angles(params_list):
angles_list = [params[:-4] for params in params_list]
x_list = [[(i + 1 / 2) * params_list[j][-1] / (len(angles_list[j])*24) for i in
range(len(angles_list[j]))] for j in
range(len(angles_list))]
plt.clf()
plt.xlabel(’time (days)’)
plt.ylabel(’Angle (radians)’)
lines = []
for i in range(len(angles_list)):
line, = plt.plot(x_list[i], angles_list[i], label=str(len(x_list[i])))
lines.append(line)
plt.legend(handles=lines)
plt.show()
def find_optimal_value_list(n_list, verbose=False, area_relevance=AREA_RELEVANCE,
time_relevance=TIME_RELEVANCE):
param_list = [find_optimal_values(n, verbose, area_relevance=area_relevance,
time_relevance=time_relevance) for n in
n_list]
return param_list
def optimizing_area_time(n=20):
x = []
y = []
d = []
s = []
for i in range(n):
params = find_optimal_values(15, True, area_relevance=i / n, time_relevance=(1 - i / n))
angles = params[:-4]
init_speed, init_angle, area, time = params[-4:]
r_fin, theta_fin, radial_speed, tangential_speed = solve_ode(params)
tang_speed_diff = tangential_speed - DISTANCE_SUN_MARS * MARS_ANGULAR_VELOCITY
s.append(sqrt(tang_speed_diff ** 2 + radial_speed ** 2) / H_TO_S)
theta_mars = MARS_EARTH_INITIAL_ANGLE_DIFF + MARS_ANGULAR_VELOCITY * params[-1]
x_diff = (r_fin * cos(theta_fin) - DISTANCE_SUN_MARS * cos(theta_mars))
y_diff = (r_fin * sin(theta_fin) - DISTANCE_SUN_MARS * sin(theta_mars))
d.append(sqrt(x_diff ** 2 + y_diff ** 2))
x.append(area)
y.append(time / 24)
print("n \t area (km^2)\t time (days) \t Dv (km/s)\t dist (km) \n")
for i in range(n):
print(i, "\t", x[i], "\t", y[i], "\t", s[i], "\t", d[i], "\n")
# apply the same algorithm bu to transfer from orbit to orbit and not from planet to planet
def inter_orbital_transport(n=N, verbose=True, initial_orbit_radius=DISTANCE_SUN_EARTH,
final_orbit_radius=DISTANCE_SUN_MARS):
global R_init, R_fin
R_init = initial_orbit_radius
R_fin = final_orbit_radius
# Define constraints for the angles values (angles must be in the range [-pi/2, pi/2])
bounds = [(-pi / 2, pi / 2) for i in range(n)]
# Define bounds for initial speed respect to Earth (from 0 to escape velocity).
bounds.append((0, EARTH_ESCAPE_VELOCITY))
# Define bounds for the initil angle formed with Earth.
bounds.append((-pi, pi))
# Define constraints for the Area of the sail (area must be between 0 and maximum
transportable weight.
bounds.append((0, 0.28))
# Define time constraints between 10 hours and 15 years
bounds.append((10, 24 * 15 * 365))
initial = [11 * pi / 180 for i in range(n)] # angles
initial.append(EARTH_ESCAPE_VELOCITY / 100) # Initial speed relative to Earth
initial.append(0) # initial angle with Earth orbit
initial.append(0.14) # Area sail.
initial.append(500 * 24) # Time.
if verbose:
results = fmin_slsqp(orbit_objective, initial, bounds=bounds, iprint=2, iter=1000)
else:
results = fmin_slsqp(orbit_objective, initial, bounds=bounds, iter=1000)
show_orbital_results(results)
return results
def test_model():
n_list = [10, 15, 20, 25, 30]
# n_list = [1, 2, 3]
param_list = find_optimal_value_list(n_list, False)
# param_list = [load_params(’learned_parameters_n={}.json’.format(n)) for n in n_list]
for params in param_list:
n = len(params) - 4
save_params(params, ’learned_parameters_n={}.json’.format(n))
# save_animated_plot(params, ’animated_plot_n={}.html’.format(n))
save_trajectory(params, ’trajectory_plot_n={}.png’.format(n))
n = input(’insert number of different angles in the simulation you want to see\n’)
params = load_params(’learned_parameters_n={}.json’.format(n))
show_animated_plot(params)
plot_multiple_angles(param_list)
def test_inter_orbital():
params = inter_orbital_transport(initial_orbit_radius=DISTANCE_SUN_MARS,
final_orbit_radius=DISTANCE_SUN_EARTH)
save_params(params, ’learned_orbital_parameters_earth_mars.json’)
save_orbital_trajectory(params, ’learned_orbital_parameters_mars_earth.png’)
plot_angles(params)
show_orbital_animated_plot(params)
print("")
if __name__ == ’__main__’:
# test_model()
test_inter_orbital()
