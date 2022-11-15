"""
Created by Adam C 12:42 14/11/2022

"""
#%%
import numpy as np
import matplotlib.pyplot as plt

pos = [[1,0,0],[-1,0,0]]
vel = [[0,1,0],[0,-1,0]]
mass= [   1   ,    1   ]
N = np.size(mass)
acc = [[1,2,3],[4,5,6 ]]
dt  = 1
r2=np.zeros(N)

tsteps = 100
print(acc[:,1])
for a in range(N):
    r2[a] = np.sum((np.abs(pos - pos[:,a].reshape((3,1))))**2,axis=0)
    #acc[:,a] = 