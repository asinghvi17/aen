# %% codecell
%matplotlib inline
# dorsal and ventral 6 neurons each side, one way connections, runs nicely
# HMH 20180722
# rechecked and run 20181115



#Imports
import numpy as np
import scipy
from scipy.integrate import odeint
import scipy.ndimage.filters as filt
from scipy.interpolate import CubicSpline

import matplotlib.pyplot as plt
import matplotlib.animation as animation

import julia

from julia import Makie
from julia import Base

import math
vent_diff = 0.9

#constants

vx0= 1
vy0= -0.49
dx0= 1.0
dy0=-0.51
vx1= 1
vy1= -0.5
dx1= 1.0
dy1=-0.5
vx2= 1
vy2= -0.5
dx2= 1
dy2=-0.5
vx3= 1
vy3= -0.49
dx3= 1.0
dy3=-0.51
vx4= 1
vy4= -0.5
dx4= 1.0
dy4=-0.5
vx5= 1
vy5= -0.5
dx5= 1
dy5=-0.5

#vector of x and y values as initial conditions
v=[vx0,vy0,dx0,dy0,vx1,vy1,dx1,dy1,vx2, vy2, dx2, dy2, vx3, vy3, dx3, dy3, vx4, vy4, dx4, dy4, vx5, vy5, dx5, dy5 ]

v

tr = (0, 2000)
t = np.linspace(tr[0], tr[1], 20001)

def f(v):
    return min(max(-2 - v, v), 2 - v)

def FHND(v, t):
    # calculates the derivatives in FHN with diffusion
    vx0, vy0, dx0, dy0, vx1, vy1, dx1, dy1, vx2, vy2, dx2, dy2, vx3, vy3, dx3, dy3, vx4, vy4, dx4, dy4, vx5, vy5, dx5, dy5 = v
    # v denotes vector (x0, y0, x1, y1, x2, y2, x3, y3) passed to the function FHND
    # the statement x0, y0, x1, y1, x2, y2, x3, y3 = v assigns components in the vector v to the
    # corresponding components.
    # Thus, do not use v = x0, y0, x1, y1, x2, y2, x3, y3
    # this allows for four neurons
    # for now, just use one-⁠way coupling from neuron 0 to neuron 1
# allow for different epsilons
    e0 = 0.08
    e1 = 0.08
    g = 0.8
    b0 = 0.46
    b1 = 0.47
 # diffusion constant
    Dhead = -0.2
    Drest= -0.02
    Dgap=0.03
    J = 0.5
    #
    dvdt=[f(vx0)- vy0+ Dhead*vent_diff*dx0                  + J, e0*(vx0-g*vy0+b0),   # neurons 0 (ventral and dorsal)
          f(dx0)- dy0+ Dhead*vx0                  + J, e0*(dx0-g*dy0+b0),   # coupled by 1-way diffusion
                                                                                    # negative diffusion constant Dhead
                                                                                    # simulated inhibitory synapse

          f(vx1)- vy1+ Drest*dx1 + Dgap*(vx0-vx1) + J, e1*(vx1-g*vy1+b1),   # neurons 1 (ventral and dorsal)
          f(dx1)- dy1+ Drest*vx1 + Dgap*(dx0-dx1) + J, e1*(dx1-g*dy1+b1),   # driven by neurons 0 through one-way
                                                                                    # diffusion (selective gap junction)
                                                                                    # Dgap
                                                                                    # Drest simulates inhibitary synapse
                                                                                    # coupling ventral and dorsam neurons 1

          f(vx2)- vy2+ Drest*dx2 + Dgap*(vx1-vx2) + J, e1*(vx2-g*vy2+b1),   # neurons 2
          f(dx2)- dy2+ Drest*vx2 + Dgap*(dx1-dx2) + J, e1*(dx2-g*dy2+b1),

          f(vx3)- vy3+ Drest*dx3 + Dgap*(vx2-vx3) + J, e1*(vx3-g*vy3+b1),   # neurons 3
          f(dx3)- dy3+ Drest*vx3 + Dgap*(dx2-dx3) + J, e1*(dx3-g*dy3+b1),

          f(vx4)- vy4+ Drest*dx4 + Dgap*(vx3-vx4) + J, e1*(vx4-g*vy4+b1),   # neurons 4
          f(dx4)- dy4+ Drest*vx4 + Dgap*(dx3-dx4) + J, e1*(dx4-g*dy4+b1),

          f(vx5)- vy5+ Drest*dx5 + Dgap*(vx4-vx5) + J, e1*(vx5-g*vy5+b1),   # neurons 5
          f(dx5)- dy5+ Drest*vx5 + Dgap*(dx4-dx5) + J, e1*(dx5-g*dy5+b1)]


    # diffusion is represented by the term D*(x1-⁠x0)
    return dvdt


def sol_vd(vd):
    global vent_diff
    vent_diff = vd
    return odeint(FHND, v, t)   #the solver


def eff_signal(sol):
    effective_signal=np.zeros((20001,24))
    for neuron_number in range(0,6):
        effective_signal[:,neuron_number]=  (sol[:,4*neuron_number]+abs(sol[:,4*neuron_number]))/2
        effective_signal[:,neuron_number]-=(sol[:,4*neuron_number+2]+abs(sol[:,4*neuron_number+2]))/2
        effective_signal[:,neuron_number]= filt.gaussian_filter1d(effective_signal[:,neuron_number],40)
        # Gaussian filter simulates a combination of muscle reponse to neural
        # stimulus; also effects of elastic properties of worm and interaction with
        # fluid

    return effective_signal

sol = sol_vd(1)

#graphs

#Comparison between neuron 0 and 1
plt.plot(t[19000:20000],(sol[19000:20000,0]), 'g', label='neuron v0 voltage')
plt.plot(t[19000:20000],(sol[19000:20000,2]), 'b', label='neuron d0 voltage')
plt.plot(t[19000:20000],(sol[19000:20000,4]), 'r', label='neuron v1 voltage')
plt.plot(t[19000:20000],(sol[19000:20000,6]), 'k', label='neuron d1 voltage')
plt.legend(loc='best')
plt.grid()
plt.show()

# ventral side
plt.plot(t[16000:20000],(sol[16000:20000,0]), 'm', label='neuron v0 voltage')
plt.plot(t[16000:20000],(sol[16000:20000,4]), 'b', label='neuron v1 voltage')
plt.plot(t[16000:20000],(sol[16000:20000,8]), 'g', label='neuron v2 voltage')
plt.plot(t[16000:20000],(sol[16000:20000,12]), 'y', label='neuron v3 voltage')
plt.plot(t[16000:20000],(sol[16000:20000,16]), 'orange', label='neuron v4 voltage')
plt.plot(t[16000:20000],(sol[16000:20000,20]), 'r', label='neuron v5 voltage')
plt.legend(loc='best')
plt.grid()
plt.show()

# dorsal side
plt.plot(t[16000:20000],(sol[16000:20000,2]), 'm', label='neuron d0 voltage')
plt.plot(t[16000:20000],(sol[16000:20000,6]), 'b', label='neuron d1 voltage')
plt.plot(t[16000:20000],(sol[16000:20000,10]), 'g',label='neuron d2 voltage')
plt.plot(t[16000:20000],(sol[16000:20000,14]), 'y', label='neuron d3 voltage')
plt.plot(t[16000:20000],(sol[16000:20000,18]), 'orange', label='neuron d4 voltage')
plt.plot(t[16000:20000],(sol[16000:20000,22]), 'r', label='neuron d5 voltage')
plt.legend(loc='best')
plt.grid()
plt.show()

# last 3 on each side
plt.plot(t[19000:20000],(sol[19000:20000,12]+1), 'm', label='neuron v3 voltage+1')
plt.plot(t[19000:20000],(sol[19000:20000,16]), 'b', label='neuron v4 voltage')
plt.plot(t[19000:20000],(sol[19000:20000,20]-1), 'g', label='neuron v5 voltage-1')
plt.plot(t[19000:20000],(sol[19000:20000,14]+1), 'y', label='neuron d3 voltage+1')
plt.plot(t[19000:20000],(sol[19000:20000,18]), 'orange', label='neuron d4 voltage')
plt.plot(t[19000:20000],(sol[19000:20000,22]-1), 'r', label='neuron d5 voltage-1')
plt.legend(loc='best')
plt.grid()
# %% codecell
effective_signal=np.zeros((20001,24))

effective_signal[:,neuron_number]=  (sol[:,4*neuron_number]+abs(sol[:,4*neuron_number]))/2

effective_signal[:, 1].shape

neuron_number = 1
for neuron_number in range(0,6):
    effective_signal[:,neuron_number]=  (sol[:,4*neuron_number]+abs(sol[:,4*neuron_number]))/2
    effective_signal[:,neuron_number]-=(sol[:,4*neuron_number+2]+abs(sol[:,4*neuron_number+2]))/2
    effective_signal[:,neuron_number]= filt.gaussian_filter1d(effective_signal[:,neuron_number],40) # Gaussian filter simulates
                                                                                                    # a combination of
                                                                                                    # muscle reponse to neural
                                                                                                    # stimulus
                                                                                                    # also effects of elastic
                                                                                                    # properties of worm and
                                                                                                    # also interaction with
                                                                                                    # fluid
    plt.plot(t[18000:19000],effective_signal[18000:19000,neuron_number], 'k', label='effective signal')
plt.legend(loc='best')
plt.grid()
plt.show()
# %% codecell
plt.plot(t[19000:20000],effective_signal[19000:20000,1], 'k', label='effective signal')
plt.legend(loc='best')
plt.grid()
plt.show()
# %% codecell
effective_signal_smooth_1 = filt.gaussian_filter1d(effective_signal[:,1],40)
plt.plot(t[19000:20000],effective_signal[19000:20000,1], 'k', label='effective signal')
plt.plot(t[19000:20000],effective_signal_smooth_1[19000:20000], 'r', label='smoothed')
plt.legend(loc='best')
plt.grid()
plt.show()
# %% codecell

number_of_bends = 6
for framenumber in range (0,200,5):
    s = 0     # arclength
    theta = 0
    svals = [s]
    thetavals= [theta]
    for j in range(0,number_of_bends):
        s += 1
        theta+= effective_signal[8000+framenumber,j]
        svals.append(s)
        thetavals.append(theta)
    print(svals)
    print(thetavals)
    # spline fit
    cstheta=CubicSpline(svals,thetavals)
    worm_coords=np.arange(0,number_of_bends+1.01,0.01)                         # I am incrementing the position in steps of
                                                                                # 0.01*(distance between nodes)
    smooththeta=cstheta(worm_coords)
    plt.plot(worm_coords,smooththeta,'k')
    plt.show()
# %% codecell
# worm_20180704 animates worms from neuronal inputs
# assumes effective_signal for 6 neurons, "time" in range 0 to 10000
# inspired by and uses https://matplotlib.org/examples/animation/double_pendulum_animated.html





number_of_joints = 6   # check #  # actually muscles #
tightness_of_bend=0.5

fig = plt.figure()#=dpi = 400=#
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-8, 8), ylim=(-8, 8))
ax.grid()
ax.set_aspect('equal')

line, = ax.plot([], [], 'o-', lw=2)


def init():
    line.set_data([], [])
    return line


def animate(framenumber):
    # begin outer loop over frames in video
    xpos=0
    ypos=0
    # s = 0         to be added
    xcoords=[xpos]
    ycoords=[ypos]
    # scoords=[s]    to be added
    # deltatheta # this will be vector of neuronal inputs
    phase=0
    theta=0



    for j in range(0,number_of_joints):
        deltatheta = effective_signal[8000+framenumber,j]
        theta+=tightness_of_bend*deltatheta
        xpos=xpos+math.cos(theta)
        ypos=ypos+math.sin(theta)
        xcoords.append(xpos)
        ycoords.append(ypos)
        # print(i,deltatheta, xpos, ypos)


    # spline fit
    nodes=np.arange(0, number_of_joints+1)
    csx=CubicSpline(nodes,xcoords)
    csy=CubicSpline(nodes,ycoords)
    worm_coords=np.arange(0,number_of_joints+1.01,0.01)                         # I am incrementing the position in steps of
                                                                                # 0.01*(distance between nodes)
    smoothwormx=csx(worm_coords)
    smoothwormy=csy(worm_coords)

    line.set_data(smoothwormx, smoothwormy)
    return line

ani = animation.FuncAnimation(fig, animate, np.arange(1, 200),
                              interval=25, blit=True, init_func=init)

# %% codecell
def save(ani, filename):
    writer = animation.FFMpegWriter(fps=15, metadata=dict(artist='Hastings, Harold and Singhvi, Anshul'), bitrate=3000)
    ani.save(filename, writer=writer)

save(ani, "smoothworm-0.8v.mp4")
# %% codecell
