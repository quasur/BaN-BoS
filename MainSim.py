#%%
"""
Created by Adam C 12:42 14/11/2022
I have no idea why
but this code is fucked beyond repair

"""

import matplotlib.pyplot as plt
import numpy as np
import astropy.constants as const
#import scipy.integrate as sci

"""
poskm = np.array([[-1.358225530315482e6,5.877370504082298e4,3.115443511303229e4]#sun
                ,[1.376549883343749e7,-6.597839966211624e7 ,6.752602566463448e06]#mercury
                ,[-4.552006624849914e6,-1.086282818036603e8,-1.276688031704240e6]#venus
                ,[5.768299035133754e7,1.353196107909931e8,2.425963145492971e4]#earth
                ,[7.280632018212381e7,2.153756053377162e8,2.724578490375057e6]#mars
                ,[7.297743868536659e8,1.180309283196366e8,-1.681667105980964e7]#jupiter
                ,[1.203381795947770e9,-8.461591498062431e8,-3.319943552438039e7]#saturn
                ,[2.013076660522243e9,2.145938180777864e9,-1.810968377389681e7]#uranus
                ,[4.449357700380400e9,-4.557354440041960e8,-9.315523521869305e7]#neptune
                ])
velkm = np.array([[8.930963297742810e-4,-1.566187437452747E-02,1.099907037042883E-4]
                ,[3.771675113483448e1,1.334760193963047e1,-2.367346535682127]#mercury
                ,[3.477141440886862e1,-1.174603245348857,-2.022144898327964]#venus
                ,[-2.779298020411351e1,1.178045149785370e1,-3.681005562539141e4]#earth
                ,[-2.198932096459520E+01,9.936623918073089,7.481016310711777e-1]#mars
                ,[-2.236133842333254,1.351047045029083e1,-6.075031230644790e-3]#jupiter
                ,[5.016422322394626,7.882308975831709,-3.364741968945446e-1]#saturn
                ,[-5.016542103919061,4.341932264280020,8.117851592618530e-2]#uranus
                ,[5.177294270500546e-1,5.439258628469002,-1.238187905959987e-1]#neptune
                ])/1000


#initial position in AU
posnat = np.array([[0   ,0,0]#sun
                ,[1     ,0,0]#earth
                ,[0.313 ,0,0]#mercury
                ,[0.723 ,0,0]#venus
                ,[1.523 ,0,0]#mars
                ,[-5.2038,0,0]#jupiter
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
                ,[0,1.947,0]#mercury
                ,[0,1.1667 ,0]#venus
                ,[0,0.8    ,0]#mars
                ,[0,-0.43567,0]#jupiter
                ,[0,0.32267,0]#saturn
                ,[0,0.22667,0]#uranus
                ,[0,0.181  ,0]#neptune
                ])

                #   sun   mer   ven   e  mar   jup    sat   ura    nep
massnat= np.array([330000,0.056,0.815,1,0.107,317.8,95.159,14.536,17.147])#in earth masses

N = np.size(massnat)

earthmass = 5.972e24 #kg
G = 6.6743e-11 #m3kg-1s-2
AU = 1.495979e8#km
#meanearthvel=2*np.pi*AU/(369.25*24*60*60)

mass = massnat#*earthmass
#pos = poskm*1000
print(pos[1])
vel = velnat*np.sqrt(mass[0]/np.linalg.norm(pos[1]))#km*1000
"""
pos = np.array([[10,0,0],
                [-10,0,0]])

vel = np.array([[0,0.2,0],
                [0,-0.2,0]])

mass= np.array([1,1])
N =2
G=1
day = 60*60*24
dt =0.1#*day
tsteps = 1000

#initialising arrays
r=np.zeros((N,3))
acc = r
accarr = np.zeros((N,N,3))
posTime = np.zeros((N,3,tsteps))
energy =np.zeros(tsteps)

#Energy calculations
def particleEnergy(rArray,vel,mass,size):
    bodyEnergy = np.zeros(size)
    for ie in range(size-1): #for each particle
        ie2 = ie+1
        r=rArray[ie,:,:]
        KE = 0.5*mass[ie]*np.sum(vel[ie]**2) #find kinetic energy from velocity
        #find distance to all other particles distance to self is very large
        GPE = -np.sum(G*mass[ie]*(mass[ie2:]/np.linalg.norm(r[ie2:])))
        bodyEnergy[ie] = KE + GPE
    return(bodyEnergy)

def movementeq(R,t): 
    R[1],-M*R[0]*((R[0]*R[0]).sum(1))

rArray = np.zeros([N,N,3])#array of each distance between each respective body (N,M) is the distance from body N to body M
#calculation of energy
for b in range(tsteps):
    for a in range(N):
        rArray[a,:] = pos[a,:]-pos[:,:]
        rArray[a,a]=10e30 #self distance is huge to avoid a /0

    energy[b] = np.sum(particleEnergy(rArray,vel,mass,N))#find energy of current system

    for a in range(N):#for each body
        posTime[a,:,b] = pos[a] #append new positions to array
        rset = rArray[a,:,:]
        rCalculation = np.linalg.norm(rset)**(-3/2)*rset
        accarr[a,:,:]= np.dot(mass.T,rCalculation)#calculate acceleration from position and massses of the 2 particles
        accarr[a,a,:] = 0 #self acceleration = 0
        
        vel = vel + np.sum(G*accarr[a,:,:],axis=0)*dt #basic ODE solve
        pos = pos + vel*dt 
        
sunPosTime = posTime[0,:,:]
relPosTime = np.zeros((N,3,tsteps))
for a in range(N):
    relPosTime[a,:,:] = posTime[a,:,:] - sunPosTime

#%%general 3d plot

fig = plt.figure()
ax = plt.axes(projection='3d')
for i in range(N):
    ax.plot3D(relPosTime[i,0,:],relPosTime[i,1,:],relPosTime[i,2,:])
plt.show()

#%%energy over time plot
plt.plot(range(tsteps),energy)
print(max(energy-energy[0]))
#%%2d plot

ax = plt.axes()
for i in range(N):
    ax.plot(relPosTime[i,0,:],relPosTime[i,1,:])

#%%2d plot for first few points

ax = plt.axes()
for i in range(N):
    ax.plot(relPosTime[i,0,0:5],relPosTime[i,1,0:5])

#%% Rudamentary animation
axanim = plt.axes()
animstep = 50
framenum = int(tsteps/animstep)
j2 = 0
j=0
k=0
while tsteps>(j2+animstep):
    k= k +1
    j2 = (k +10)*animstep
    j = k*animstep
    axanim.clear()
    for i in range(N):  
        axanim.plot(relPosTime[i,0,j:j2],relPosTime[i,1,j:j2],'-')
        axanim.plot(relPosTime[i,0,j2],relPosTime[i,1,j2],'o')
    plt.pause(0.1)
    
#%%
