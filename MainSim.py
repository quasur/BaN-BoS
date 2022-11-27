#%%
"""
Created by Adam C 12:42 14/11/2022

"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as sci

#initial position in AU
posnat = np.array([[0   ,0,0]#sum
                ,[1     ,0,0]#earth
                ,[0.307 ,0,0]#mercury
                ,[0.723 ,0,0]#venus
                ,[1.523 ,0,0]#mars
                ,[5.2038,0,0]#jupiter
                ,[9.5826,0,0]#saturn
                ,[19.191,0,0]#uranus
                ,[30.07 ,0,0]#neptune
                ])
pos = posnat*10000 #scale to improve accuracy of the sim
#sim loses precision when vel mass and pos are all very different
#orders of magnitude

#vel in earth velocities
velnat = np.array([[0,0    ,0]#sun
                ,[0,1      ,0]#earth
                ,[0,1.9667 ,0]#mercury
                ,[0,1.1667 ,0]#venus
                ,[0,0.8    ,0]#mars
                ,[0,0.43567,0]#jupiter
                ,[0,0.32267,0]#saturn
                ,[0,0.22667,0]#uranus
                ,[0,0.181  ,0]#neptune
                ])
                #   sun   e  mer   ven   mar   jup    sat   ura    nep
massnat= np.array([330000,1,0.056,0.815,0.107,317.8,95.159,14.536,17.147])#in earth masses
#scaling mass and velocity
mass = massnat/100
vel = velnat * np.sqrt(mass[0]/pos[1,0])#scaling factor derived from equating gravitational and centripetal forces

#initialising arrays
N = np.size(mass)
r=np.zeros((N,3))
acc = r
accarr = np.zeros((N,N,3))

dt = 1
tsteps = 1000000
posTime = np.zeros((N,3,tsteps))


#Calculation of position magnitudes and normalised unit vectors
def rHatCalc(pos,body):
    r = pos[body,:]-pos[:,:]
    r[body,:] = r[body,:]+10**10
    #replaces the self distance with a 
    #very large one so it is neglegible for 1/f(r) calculations
    rnorm=np.linalg.norm(r,axis = 1)
    rhat = -r/rnorm[:,None]
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
    #energy[b] = np.sum(particleEnergy(pos,vel,mass,N))#find energy of current system
    for a in range(N):#for each body
        rdir,r2 = rHatCalc(pos,a) #calculate distance (distance to self is very large)
        for c in range(N):#for each other body
            if a == c:
                accarr[a,c,:] = [0,0,0]#acceleration to self is 0
            else:
                accarr[a,c,:]= ((mass[c]/r2[c])*rdir[c,:])#calculate acceleration from position and massses of the 2 particles
        
        acc = np.sum(accarr,axis=1) #net acceleration
        vel = vel + acc*dt #basic ODE solve
        pos = pos + vel*dt 
        posTime[a,:,b] = pos[a]/10000 #append new positions to array
        acc = np.sum(accarr,axis=1)

#%%

fig = plt.figure()
ax = plt.axes(projection='3d')
for i in range(N):
    ax.plot3D(posTime[i,0,:],posTime[i,1,:],posTime[i,2,:])

plt.show()
#%%
plt.plot(range(tsteps),energy)

#%%

ax = plt.axes()
for i in range(N):
    ax.plot(posTime[i,0,:],posTime[i,1,:])