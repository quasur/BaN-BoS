#%%
"""
Created by Adam C 12:42 14/11/2022

"""

import numpy as np
import matplotlib.pyplot as plt

pos = np.array([[1,0,0],[-1,0,0]])
vel = np.array([[0,1,0],[0,-1,0]])
mass= np.array([   1   ,    1   ])
N = np.size(mass)
r=np.zeros((N,3))
acc = r
vel = acc

dt = 1
tsteps = 100
posTime = np.zeros((N,3,tsteps))

for b in range(tsteps):
    for a in range(N):
        r = pos[a,:]-pos[:,:]
        r[a,:] = r[a,:]+10000
        rnorm=np.linalg.norm(r,axis = 1)
        rdir = r/rnorm[:,None]
        r2 = np.sum(r**2,axis=1)
        acc[a,:]= np.dot((mass/r2),(rdir))
        
        vel = vel + acc*dt
        pos = pos + vel*dt
        posTime[a,:,b] = pos[a]
        
        
t = np.arange(0,tsteps)



# %%
