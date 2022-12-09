#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
print("Select mode: \n 1 - Fast mode, Simulates only sun, jupiter and 2nd star to save time \n 2 - Normal mode \n 3 - Basic mode, the basic N = 2 mode of the sim")
mode=input()
print("mode: ",mode)
#data from nasa horizons nov 29th
posm = np.array([[-1.358225530315482e6,5.877370504082298e4,3.115443511303229e4]#sun
                ,[1.376549883343749e7,-6.597839966211624e7 ,6.752602566463448e06]#mercury
                ,[-4.552006624849914e6,-1.086282818036603e8,-1.276688031704240e6]#venus
                ,[5.768299035133754e7,1.353196107909931e8,2.425963145492971e+4]#earth
                ,[7.280632018212381e7,2.153756053377162e8,2.724578490375057e6]#mars
                ,[7.297743868536659e8,1.180309283196366e8,-1.681667105980964e7]#jupiter
                ,[1.203381795947770e9,-8.461591498062431e8,-3.319943552438039e7]#saturn
                ,[2.013076660522243e9,2.145938180777864e9,-1.810968377389681e7]#uranus
                ,[4.449357700380400e9,-4.557354440041960e8,-9.315523521869305e7]#neptune
                ,[0,0,1e11]#interloper
                ])*1000#converting from km to m


velms = np.array([[8.930963297742810e-4,-1.566187437452747E-02,1.099907037042883E-4]
                ,[3.771675113483448e1,1.334760193963047e1,-2.367346535682127]#mercury
                ,[3.477141440886862e1,-1.174603245348857,-2.022144898327964]#venus
                ,[-2.779298020411351e1,1.178045149785370e+1,-3.681005562539141e-4]#earth
                ,[-2.198932096459520E+01,9.936623918073089,7.481016310711777e-1]#mars
                ,[-2.236133842333254,1.351047045029083e1,-6.075031230644790e-3]#jupiter
                ,[5.016422322394626,7.882308975831709,-3.364741968945446e-1]#saturn
                ,[-5.016542103919061,4.341932264280020,8.117851592618530e-2]#uranus
                ,[5.177294270500546e-1,5.439258628469002,-1.238187905959987e-1]#neptune
                ,[0,0.125*2**(-0.5),-8*2**(-0.5)]#interloper
                ])*1000


                   #sol       #mercury #venus      #earth      #mars      #jupiter    #saturn   #uranus   #neptune  #other star
masskg = np.array([1.988500e30,3.302e23,4.8685e24,5.97219e24,6.4171e23,1.89818722e27,5.6834e26,8.6813e25,1.02409e26,2e28])


if mode == '3':
    pos = np.array([[-1,0,0],[1,0,0]])
    mass = np.array([1,1])
    vel = np.array([[0,-0.5,0],[0,0.5,0]])
    print("simple case")
else:
    pos = posm
    vel = velms
    G = 6.6743e-11
    massfactor = np.linalg.norm(posm[2,:])*np.linalg.norm(velms[2,:]**2)/G
    mass = (masskg/masskg[0])*massfactor
    escapevel = np.sqrt(2*G*mass[0]/np.linalg.norm(pos[-1]))/100
    vel[-1] = vel[-1]*escapevel


if mode == '1':
    pos = np.array([pos[0],pos[5],pos[9]])
    vel = np.array([vel[0],vel[5],vel[9]])
    mass = np.array([mass[0],mass[5],mass[9]])
elif mode == '2':
    print("standard mode")

N = np.size(mass)


day = 24*60**2
if mode == '1':
    dt = day
    timestep=round(366*100)
elif mode=='2':
    dt = day/10
    timestep = 60225*10
elif mode == '3':
    dt = 1
    timestep = 2
else:
    dt=1
    timestep=0

print(timestep)

posTime = np.zeros([N,timestep,3])
accArray = np.zeros([N,N,3])
massmult = np.transpose(np.array([mass,mass,mass]))
GPEArray = np.zeros([N])
KEArray = np.zeros([N])
energyTime = np.zeros(timestep)
accTime = np.zeros(timestep)
rCalc = np.zeros([N,N,3])

for t in range(timestep):
    bigpos = pos*np.ones([N,N,3])
    swapbigpos = np.swapaxes(bigpos,0,1)
    rArray = bigpos-swapbigpos

    posTime[:,t,:] = pos[:,:]
    for i in range(N):
        rArray[i,i,:]=3**(-0.5)

    rScale = np.linalg.norm(rArray,axis=2)**(3)
    for i in range(N):
        rScale[i,i] = 1e100
        rArray[i,i,:] = 1e100
    rCalc[:,:,0] = rArray[:,:,0]/rScale
    rCalc[:,:,1] = rArray[:,:,1]/rScale
    rCalc[:,:,2] = rArray[:,:,2]/rScale

    for i in range(N):
        accArray[i,:,:]=G*massmult*rCalc[i,:,:]
        accArray[i,i,:]=0#self acceleration = 0

    #energy calculation
    for i in range(N-1):
        GPEArray[i] = -G*mass[i]*np.sum(mass[i+1:]/np.linalg.norm(rArray[i,i+1:,:],axis=1))
    KEArray = 0.5*mass*np.linalg.norm(vel[:,:],axis=1)**2
    energyTime[t] = np.sum(GPEArray+KEArray)
    acc = np.sum(accArray,axis=1)
    vel = vel + acc*dt
    pos = pos + vel*dt
    accTime[t] = np.linalg.norm(acc[-1,:])
    
sunPosTime = posTime[0,:,:]
relPosTime = (posTime-sunPosTime)/1.496e+11

if mode == '3':
    plt.plot(posTime[0,:,0],posTime[0,:,1],'r-')
    plt.plot(posTime[1,:,0],posTime[1,:,1])
else:
    fig = plt.figimage
    ax = plt.axes(projection='3d')
    for i in range(N):
        ax.plot3D(relPosTime[i,:,0],relPosTime[i,:,1],relPosTime[i,:,2])

#%%

fig = plt.figimage
ax = plt.axes(projection='3d')
ax.plot3D(relPosTime[0,:,0],relPosTime[0,:,1],relPosTime[0,:,2],'r-')
ax.plot3D(relPosTime[1,:,0],relPosTime[1,:,1],relPosTime[1,:,2],'b-')
ax.plot3D(relPosTime[-1,:,0],relPosTime[-1,:,1],relPosTime[-1,:,2],'g-')

#%%
for i in range(N-1):
    plt.plot(relPosTime[i,:,0],relPosTime[i,:,1],'-')
    maxacc = np.argmax(accTime)
    plt.plot(relPosTime[2,maxacc,0],relPosTime[2,maxacc,1],'go')
    plt.plot(0,0,'ro')

#%%
#%%2D animation

axanim = plt.axes()
animstep = 50
framenum = int(timestep/animstep)
j2=0
j=0
k=0
bodies = 1+4
while timestep>(j2+animstep):
    k= k +1
    j2 = (k +10)*animstep #time scale of animation
    j = k*animstep
    axanim.clear()
    axesLimit = np.abs(np.max(relPosTime[0:bodies,:,:]))*1.1
    axanim.set_xlim(-axesLimit,axesLimit)
    axanim.set_ylim(-axesLimit,axesLimit)
    for i in range(bodies):  
        axanim.plot(relPosTime[i,j:j2,0],relPosTime[i,j:j2,1],'-')
        axanim.plot(relPosTime[i,j2,0],relPosTime[i,j2,1],'o')
    plt.pause(0.1)
    
#%%energy graph
plt.plot(range(timestep),energyTime/energyTime[0]-1)

#%%
plt.plot(range(timestep),accTime)