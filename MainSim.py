#%%
"""
Created by Adam C 12:42 14/11/2022

"""

import matplotlib.pyplot as plt
import numpy as np

pos = np.array([[0,0,5],
                [100,0.1,0],
                [0,300,0]])

vel = np.array([[0,-0.02,0],
                [0,0.2,0],
                [0.1,0,0]])

mass= np.array([10,1,0.1])
N = np.size(mass)
r=np.zeros((N,3))
acc = r
accarr = np.zeros((N,N,3))

dt = 0.1
tsteps = 50000
posTime = np.zeros((N,3,tsteps))


#Calculation of position magnitudes and normalised unit vectors
def rHatCalc(pos,body):
    r = pos[body,:]-pos[:,:]
    r[body,:] = r[body,:]+10**40
    #replaces the self distance with a 
    #very large one so it is neglegible for 1/f(r) calculations
    rnorm=np.linalg.norm(r,axis = 1)
    rhat = r/rnorm[:,None]
    #finds normalised vecot
    rsqruared = np.sum(r**2,axis=1)
    #finds magnitude squared of displacement vectors
    return(rhat,rsqruared)

#Energy calculations
def particleEnergy(pos,vel,mass,size):
    bodyEnergy = np.zeros(size)
    for ie in range(size): #for each particle
        KE = 0.5*mass[ie]*np.sum(vel[ie]**2) #find kinetic energy from velocity
        r = rHatCalc(pos,ie)[1]**(1/2) #find distance to all other particles distance to self is very large
        GPE = np.sum(mass[ie]*(mass/r))
        bodyEnergy[ie] = KE + GPE
    return(bodyEnergy)

energy =np.zeros(tsteps)

#calculation of energy
for b in range(tsteps):
    energy[b] = np.sum(particleEnergy(pos,vel,mass,N))#find energy of current system
    for a in range(N):#for each body
        rdir,r2 = rHatCalc(pos,a) #calculate distance (distance to self is very large)
        for c in range(N):#for each other body
            accarr[a,c,:]= -((mass[c]/r2[c])*rdir[c,:])#calculate acceleration from position and massses of the 2 particles
        
        acc = np.sum(accarr,axis=1) #net acceleration
        vel = vel + acc*dt #basic ODE solve
        pos = pos + vel*dt 
        posTime[a,:,b] = pos[a] #append new positions to array
        



fig = plt.figure()
ax = plt.axes(projection='3d')
for i in range(N):
    ax.plot3D(posTime[i,0,:],posTime[i,1,:],posTime[i,2,:])
        
        
t = np.arange(0,tsteps)




plt.show()
#%%
plt.plot(range(tsteps),energy)