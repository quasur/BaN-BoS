#%%
"""
Created by Adam C 12:42 14/11/2022

"""

import numpy as np
import matplotlib.pyplot as plt

pos = np.array([[0,0,0],[50,0,0]])
vel = np.array([[0,0.05,0],[0,-0.5,0]])
mass= np.array([10,1])
N = np.size(mass)
r=np.zeros((N,3))
acc = r
accarr = np.zeros((N,N,3))

dt = 0.1
tsteps = 50000
posTime = np.zeros((N,3,tsteps))
fig = plt.figure()
ax = plt.axes()#projection='3d')


for b in range(tsteps):
    for a in range(N):
        r = pos[a,:]-pos[:,:]
        r[a,:] = r[a,:]+10**15
        rnorm=np.linalg.norm(r,axis = 1)
        rdir = r/rnorm[:,None]
        r2 = np.sum(r**2,axis=1)
        for c in range(N):
            accarr[a,c,:]= -((mass[c]/r2[c])*rdir[c,:])
        
        acc = np.sum(accarr,axis=1)
        vel = vel + acc*dt
        #print(vel)
        pos = pos + vel*dt
        posTime[a,:,b] = pos[a]
        

plt.plot(posTime[0,0,:],posTime[0,1,:],'b.')
plt.plot(posTime[1,0,:],posTime[1,1,:],'r.')

        
        
        
t = np.arange(0,tsteps)




plt.show()
# %%
