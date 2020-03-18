# %%
#basic program
# dorsal and ventral 6 neurons each side, one way connections, runs nicely
# HMH 20180722

#Imports
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math
from scipy.signal import correlate
from scipy.interpolate import CubicSpline

import matplotlib

from matplotlib.backends.backend_pgf import FigureCanvasPgf
matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)

#constants
dx0= 1.0
dy0= -0.51
vx0= 1
vy0= -0.49
dx1= 1.0
dy1=-0.5
vx1= 1
vy1=-0.5
dx2= 1
dy2=-0.5
vx2= 1
vy2= -0.5
dx3= 1.0
dy3=-0.51
vx3= 1
vy3= -0.49
dx4= 1.0
dy4=-0.5
vx4= 1
vy4=-0.5
dx5= 1
dy5=-0.5
vx5= 1
vy5= -0.5


#vector of x and y values
v=[vx0,vy0,dx0,dy0,vx1,vy1,dx1,dy1,vx2, vy2, dx2, dy2, vx3, vy3, dx3, dy3, vx4, vy4, dx4, dy4, vx5, vy5, dx5, dy5 ]

#tests for times
#t = np.linspace(0, 20000, 100001)
t = np.linspace(0, 2000, 20001)

def FHND(v, t):
    # calculates the derivatives in FHN with diffusion
    vx0, vy0, dx0, dy0, vx1, vy1, dx1, dy1, vx2, vy2, dx2, dy2, vx3, vy3, dx3, dy3, vx4, vy4, dx4, dy4, vx5, vy5, dx5, dy5 = v
    # v denotes vector (x0, y0, x1, y1, x2, y2, x3, y3) passed to the function FHND
# allow for different epsilons
    e0 = 0.08
    e1 = 0.08
    g = 0.8
    b0 = 0.46
    b1 = 0.47
 # diffusion constant
    Dhead = -0.2
    Drest=-0.02
    Dgap=0.05
    J = 0

# sine wave drives neuron 0
    dvdt=dvdt=[vx0-(vx0**3/3)- vy0+ Dhead*max(dx0-vx0,0)    + J, e0*(vx0-g*vy0+b0),   # neurons 0 (ventral and dorsal)
          dx0-(dx0**3/3)- dy0+ Dhead*max(vx0-dx0,0)         + J, e0*(dx0-g*dy0+b0),   # coupled by 1-way diffusion
                                                                                    # negative diffusion constant Dhead
                                                                                    # simulated inhibitory synapse

          vx1-(vx1**3/3)- vy1+ Drest*max(dx1-vx1,0) + Dgap*max(vx0-vx1,0) + J, e1*(vx1-g*vy1+b1),   # neurons 1 (ventral and dorsal)
          dx1-(dx1**3/3)- dy1+ Drest*max(vx1-dx1,0) + Dgap*max(dx0-dx1,0) + J, e1*(dx1-g*dy1+b1),   # driven by neurons 0 through one-way
                                                                                    # diffusion (selective gap junction)
                                                                                    # Dgap
                                                                                    # Drest simulates inhibitary synapse
                                                                                    # coupling ventral and dorsam neurons 1

          vx2-(vx2**3/3)- vy2+ Drest*max(dx2-vx2,0) + Dgap*max(vx1-vx2,0) + J, e1*(vx2-g*vy2+b1),   # neurons 2
          dx2-(dx2**3/3)- dy2+ Drest*max(vx2-dx2,0) + Dgap*max(dx1-dx2,0) + J, e1*(dx2-g*dy2+b1),

          vx3-(vx3**3/3)- vy3+ Drest*max(dx3-vx3,0) + Dgap*max(vx2-vx3,0) + J, e1*(vx3-g*vy3+b1),   # neurons 3
          dx3-(dx3**3/3)- dy3+ Drest*max(vx3-dx3,0) + Dgap*max(dx2-dx3,0) + J, e1*(dx3-g*dy3+b1),

          vx4-(vx4**3/3)- vy4+ Drest*max(dx4-vx4,0) + Dgap*max(vx3-vx4,0) + J, e1*(vx4-g*vy4+b1),   # neurons 4
          dx4-(dx4**3/3)- dy4+ Drest*max(vx4-dx4,0) + Dgap*max(dx3-dx4,0) + J, e1*(dx4-g*dy4+b1),

          vx5-(vx5**3/3)- vy5+ Drest*max(dx5-vx5,0) + Dgap*max(vx4-vx5,0) + J, e1*(vx5-g*vy5+b1),   # neurons 5
          dx5-(dx5**3/3)- dy5+ Drest*max(vx5-dx5,0) + Dgap*max(dx4-dx5,0) + J, e1*(dx5-g*dy5+b1)]
    # diffusion is represented by terms like D*(x1-⁠x0)
    return dvdt
#the solver
sol = odeint(FHND, v, t)


# find period by cross-correlation
position_of_peak1 = np.argmax(correlate(sol[8000:10000,4],sol[8000:10000,4])[2010:])+2010
position_of_peak0 = np.argmax(correlate(sol[8000:10000,4],sol[8000:10000,4]))
position_of_peak_minus1 = np.argmax(correlate(sol[8000:10000,4],sol[8000:10000,4])[:1990])
period = (position_of_peak1 - position_of_peak_minus1)/2

# find timelag_2_to_3 by cross-correlation
lag_2_to_3 = np.argmax(correlate(sol[8000:10000,4],sol[8000:10000,6])) # potentially offset by constant
lag_3_to_2 = np.argmax(correlate(sol[8000:10000,6],sol[8000:10000,4])) # potentially offset by same constant
timelag_2_to_3 = (lag_3_to_2 - lag_2_to_3)/2

# find phaselag_2_to_3
phaselag_2_to_3 = timelag_2_to_3/ period  # range 0 to 1
phaselag_2_to_3_deg = 360*phaselag_2_to_3



xpos=0
ypos=0
xcoords=[xpos]
ycoords=[ypos]
# deltatheta # this will be vector of neuronal inputs
phase=0
theta=0                              # this will be head angle

number_of_joints = 5				# 6 neurons
total_phase_shift=1.5*math.pi
phase_shift_per_joint=total_phase_shift/number_of_joints
tightness_of_bend=0.25			# probably 0.25 if neurons drive it

import csv

for t in np.arange(8000,8400,1):
    xpos=0
    ypos=0
    xcoords=[xpos]
    ycoords=[ypos]
    v = [sol[t,0], sol[t,2], sol[t,4], sol[t,6], sol[t,8],sol[t,10], sol[t,12],sol[t,14],sol[t,16],
        sol[t,18],sol[t,20],sol[t,22]] # from neuron
    phase=0
    theta=0

    for i in range(0,number_of_joints+1):
        deltatheta = phase_shift_per_joint*v[i]
        theta+=tightness_of_bend*deltatheta
        xpos=xpos+math.cos(theta)
        ypos=ypos+math.sin(theta)
        xcoords.append(xpos)
        ycoords.append(ypos)
       # print(i,deltatheta, xpos, ypos)

# spline fit
    nodes=np.arange(0, number_of_joints+2)
    csx=CubicSpline(nodes,xcoords)
    csy=CubicSpline(nodes,ycoords)
    worm_coords=np.arange(0,number_of_joints+1.01,0.01)                         # I am incrementing the position in steps of
                                                                            # 0.01*(distance between nodes)
    smoothwormx=csx(worm_coords)
    smoothwormy=csy(worm_coords)

    with open("CSVworms/worm-{}.csv".format(t), "w") as csvfile:

        writer = csv.writer(csvfile)

        writer.writerows([[xcoords[i], ycoords[i]] for i in range(len(xcoords))])

    #plt.plot(xcoords,ycoords,'b', label='stickworm', marker='o')

    plt.ylim(-6,6)
    plt.xlim(-6,6)# draw worm
    plt.plot(smoothwormx,smoothwormy,'k', label='smoothworm',color='k')
    plt.plot(xcoords,ycoords, linestyle=' ', marker='o', color='red')# draw worm
    plt.axis('off')
    #plt.legend(loc='best')
    #plt.grid()
    plt.show()
# %%

# %%
