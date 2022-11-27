#%%
"""
Created by Adam C 12:42 14/11/2022
Basic case of 2 body system
"""

import matplotlib.pyplot as plt
import numpy as np

pos = np.array([[10,0,0],
                [-10,0,0]])

vel = np.array([[0,0.2,0],
                [0,-0.2,0]])

mass= np.array([1,1])
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
    r[body,:] = r[body,:]+10**10
    #replaces the self distance with a 
    #very large one so it is neglegible for 1/f(r) calculations
    rnorm=np.linalg.norm(r,axis = 1)
    rhat = r/rnorm[:,None]
    #finds normalised vecot
    rsqruared = np.sum(r**2,axis=1)
    #finds magnitude squared of displacement vectors
    return(rhat,rsqruared)

#Initial energy calculations
def particleEnergy(pos,vel,mass,size):
    bodyEnergy = np.zeros(size)
    for ie in range(size):
        KE = 0.5*mass[ie]*np.sum(vel[ie]**2)
        r = rHatCalc(pos,ie)[1]**(1/2)
        GPE = np.sum(mass[ie]*(mass/r))
        bodyEnergy[ie] = KE + GPE
    return(bodyEnergy)

initialEnergy = np.sum(particleEnergy(pos,vel,mass,N))

for b in range(tsteps):
    for a in range(N):
        rdir,r2 = rHatCalc(pos,a)
        for c in range(N):
            if c == a:
                accarr[a,c,:] = 0
            else:
                accarr[a,c,:]= -((mass[c]/r2[c])*rdir[c,:])
        
        acc = np.sum(accarr,axis=1)
        vel = vel + acc*dt
        #print(vel)
        pos = pos + vel*dt
        posTime[a,:,b] = pos[a]
        

finalEnergy = np.sum(particleEnergy(pos,vel,mass,N))

print(finalEnergy-initialEnergy)


fig = plt.figure()
ax = plt.axes()
for i in range(N):
    ax.plot(posTime[i,0,:],posTime[i,1,:])
        
        
t = np.arange(0,tsteps)




plt.show()
#%%
