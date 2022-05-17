# -*- coding: utf-8 -*-
"""
PH322: Celestial Mechanics 
Final Project: Non-impulsive orbital maneuver from Earth to Saturn orbit
May 17th 2022, Rose-Hulman Institute of Technology
Authors: Ishaan Mishra and Lee Trent
"""
import math
import matplotlib.pyplot as pl
def printbreak():
    print("---------------------------------")
    
#--------------------------initial_calculations--------------------------#

mu=1.327e20#for Solar orbit (m^3/s^2)
rp=149.6e9#m
ra=1432e9#m


#T=2*235e-3#N <<Using data from NASA NEXT-C Ion thruster
T=77.7e-3#<<maximum thrust for TRANSFER ORBIT
Isp=4220#s

m0=1425.6 # k
g0=9.81#m/s^2

a=(mu/ra)**(1/2)-(mu/rp)**(1/2)
t=m0*g0*Isp/T*(1-math.exp(a/(Isp*g0)))


mdot=T/(Isp*g0)
mf=m0-t*mdot
printbreak()
print("INITIAL CALCULATIONS (TRANSFER ORBIT)")
print("Initial mass:",m0,"kg")
print("Final mass:", mf, "kg")
print("Time taken to reach Saturn orbit from Earth orbit: ", t/60/60/24/365, "years")

printbreak()

#--------------------------function_variables--------------------------#

dt= 100#s
#Assuming Earth orbit is circular, initial conditions 
a0=[mu/rp**2,0]
r0=[rp,0]
v0=[0,math.sqrt(mu/rp)]

a=a0
r=r0
v=v0
m=2108.3#spacecraft initial mass
AU=1.496e11

#Initial datasets for graphing
X=[r[0]/AU,]
Y=[r[1]/AU,]

Mx=[0,]
My=[m,]

Vx=[0,]
Vy=[v[1]/1000,]

Ag=[]
At=[]
Ax=[]

delta_v=0
at=0
ag=0

#-------------------------------functions------------------------------#

def unit_vector(v):
    vmag=v_mag(v)
    v=[v[0]/vmag,v[1]/vmag]
    return(v)

def v_mag(v):
    vmag=math.sqrt(v[0]**2+v[1]**2)
    return vmag

def update_a():
    global a, dt, r, v, t_flight, delta_v, at, ag
    rh=unit_vector(r)
    ag=[-mu/((v_mag(r))**2)*rh[0], -mu/((v_mag(r))**2)*rh[1]]
    Ag.append(v_mag(ag))
    vh=unit_vector(v)
    '''
    if t_flight>365*24*60*60: 
        a2=[0,0,0]
    else:
        a2=[T/m*vh[0],T/m*vh[1]]
    '''
    at=[T/m*vh[0],T/m*vh[1]]
    At.append(v_mag(at))
    a=[ag[0]+at[0],ag[1]+at[1]]
    delta_v+=v_mag(at)*dt
    
def update_m():
    global m, mdot, dt, Mx, My
    m-=mdot*dt
    Mx.append(t_flight/60/60/24)
    My.append(m)
    
def update_v(): 
    global a, dt, v, t_flight
    v=[v[0]+a[0]*dt,v[1]+a[1]*dt]
    Vx.append(t_flight/60/60/24)
    Vy.append(v_mag(v)/1000)

def update_r():
    global v, dt, r 
    r= [r[0]+v[0]*dt, r[1]+v[1]*dt]
    X.append(r[0]/AU)
    Y.append(r[1]/AU)
    Ax.append(v_mag(r)/AU)
    
def plot_trajectory(): 
    ax = pl.gca()
    ax.plot(X,Y, color='purple')
    circle1 = pl.Circle((0, 0), ra/AU, color='b', fill=False)
    circle2 = pl.Circle((0, 0), rp/AU, color='r', fill=False)
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.set_xlabel('x (AU)')
    ax.set_ylabel('y (AU)')
    ax.set_title('Trajectory of NASA NEXT-C ion thruster-powered Saturn flyby mission')
    ax.set_aspect('equal', adjustable='datalim')

    pl.show()

def plot_speed(): 
    pl.plot(Vx, Vy)
    pl.xlabel('Time (days)')
    pl.ylabel('Speed (km/s)')
    pl.title('Speed of spacecraft with respect to time')
    pl.show()

def plot_mass(): 
    pl.plot(Mx, My)
    pl.xlabel('Time (days)')
    pl.ylabel('Mass (kg)')
    pl.title('Mass of spacecraft with respect to time')
    pl.show()
    
def plot_acc(): 
    pl.plot(Ax, Ag, color='r')
    pl.plot(Ax,At, color='b')
    pl.xlabel('Distance of spacecraft from Sun (AU)')
    pl.ylabel('Acceleration magnitude (m/s^2)')
    pl.title('Acceleration magnitude of spacecraft with respect to distance')
    pl.show()
#-------------------------------main_loop------------------------------#

print("NUMERICAL SIMULATION ")
print("Time step:", dt, "s")
print("Initial mass:",m, "kg")
t_flight=0
while v_mag(r)<ra:
    t_flight+=dt
    update_a()
    update_m()
    update_v()
    update_r()
print("Final mass:",m, "kg")
print("Initial position:",r0, "m")
print("Final position:", r, "m")
print("Initial velocity:",v0, "m/s")
print("Final velocity: ",v, "m/s")
print("Time taken to reach Saturn orbit from Earth orbit: ", t_flight/60/60/24/365, "years")
print("Delta-v:", delta_v/1000, "km/s")
if v_mag(at)>v_mag(ag): 
    print("NOT A TRANSFER ORBIT")
else: 
    print("TRANSFER ORBIT")
printbreak()
#print(math.sqrt(mu/r), v_mag(v))
#print(math.sqrt(mu/ra)/1000, v_mag(v)/1000)
plot_trajectory()
plot_speed()
plot_mass()
plot_acc()

#----------------------------------------------------------------------#
