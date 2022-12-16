#BaN-BoS, a forever WIP by Adam Corness 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
print("Select mode: \n 1 - Fast mode, Simulates only sun, jupiter and 2nd star to save time \n 2 - Normal mode \n 3 - Basic mode, the basic N = 2 mode of the sim")
mode=input()
print("mode: ",mode)
matplotlib.use("TkAgg") #this mode allowed me to display multiple figures in vscode without it crashing
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
                ,[778.479e6/2,0,1e11]#interloper
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
                ,[0.000001,0,-0.1]#interloper
                ])*1000

                   #sol       #mercury #venus      #earth      #mars      #jupiter    #saturn   #uranus   #neptune  #other star
masskg = np.array([1.988500e30,3.302e23,4.8685e24,5.97219e24,6.4171e23,1.89818722e27,5.6834e26,8.6813e25,1.02409e26,2e30])

if mode == '3':#overrides for basic case
    pos = np.array([[-1,0,0],[1,0,0]])
    mass = np.array([1,1])
    vel = np.array([[0,-0.4,0],[0,0.4,0]])
    print("simple case")
    G=1
else:#scaling initial positions to make the simulation actually have a model solar system rather than multiple unbound and eccentric planets
    pos = posm
    vel = velms
    G = 6.6743e-11
    massfactor = np.linalg.norm(posm[2,:])*np.linalg.norm(velms[2,:]**2)/G
    mass = (masskg/masskg[0])*massfactor
    escapevel = np.sqrt(2*G*mass[0]/np.linalg.norm(pos[-1]))
    vel[-1] = vel[-1]*escapevel


if mode == '1': #cut down arrays for timesave
    pos = np.array([pos[0],pos[5],pos[9]])
    vel = np.array([vel[0],vel[5],vel[9]])
    mass = np.array([mass[0],mass[5],mass[9]])
elif mode == '2':
    print("standard mode")

N = np.size(mass)

#specify time parameters
day = 24*60**2
if mode == '1':
    dt = day/10
    timestep=round(366*1000)
elif mode=='2':
    dt = day/40
    timestep = 60225*40
elif mode == '3':
    dt = 0.0001
    timestep = 200000
else:
    dt=1
    timestep=0

print("timesteps: " ,timestep)
print("number of bodies: ",N)
#initialise arrays for data storage
posTime = np.zeros([N,timestep,3])
accArray = np.zeros([N,N,3])
massmult = np.transpose(np.array([mass,mass,mass]))
GPEArray = np.zeros([N])
KEArray = np.zeros([N])
energyTime = np.zeros(timestep)
accTime = np.zeros(timestep)
rCalc = np.zeros([N,N,3])


###=================###
### MAIN SIMULATION ###
###=================###

for t in range(timestep):
    bigpos = pos*np.ones([N,N,3])
    swapbigpos = np.swapaxes(bigpos,0,1) #array of distances from body a to b stored for exaple as ([b,a,[x,y,z]])
    rArray = bigpos-swapbigpos

    posTime[:,t,:] = pos[:,:]
    for i in range(N):
        rArray[i,i,:]=3**(-0.5) #calculate r part in a = GM/r^2 * rhat

    rScale = np.linalg.norm(rArray,axis=2)**(3)
    for i in range(N):
        rScale[i,i] = 1e100
        rArray[i,i,:] = 1e100
    rCalc[:,:,0] = rArray[:,:,0]/rScale
    rCalc[:,:,1] = rArray[:,:,1]/rScale #there is probably a more efficient way of doing this calculation with numpy
    rCalc[:,:,2] = rArray[:,:,2]/rScale #but the time saved in simulation runtime will probably be less than the time it takes to learn how to implement it

    for i in range(N):
        accArray[i,:,:]=G*massmult*rCalc[i,:,:]
        accArray[i,i,:]=0#self acceleration = 0 

    #energy calculation
    for i in range(N-1):
        GPEArray[i] = -G*mass[i]*np.sum(mass[i+1:]/np.linalg.norm(rArray[i,i+1:,:],axis=1))
    KEArray = 0.5*mass*np.linalg.norm(vel[:,:],axis=1)**2
    energyTime[t] = np.sum(GPEArray+KEArray)

    #update vel and pos based on the acceleration measured
    acc = np.sum(accArray,axis=1)
    vel = vel + acc*dt
    pos = pos + vel*dt
    accTime[t] = np.linalg.norm(acc[-1,:])
    
sunPosTime = posTime[0,:,:]
relPosTime = (posTime-sunPosTime)/1.496e+11 #convert to AU centered on the suns position

###==================###
### PLOTTING SECTION ###
###==================###

if mode == '3': #plot for simple case
    plt.figure(figsize=(9,9),dpi=160)
    axsimple = plt.axes([0.1,0.1,0.8,0.8])
    axsimple.set_xlim(-1.05,1.05)
    axsimple.set_ylim(-1.05,1.05)
    axsimple.set_title("Paths of bodies in simple case",fontsize=18)
    axsimple.set_xlabel("x position",fontsize=16)
    axsimple.set_ylabel("y position",fontsize=16)
    axsimple.plot(posTime[0,:,0],posTime[0,:,1],'r-')
    axsimple.plot(posTime[1,:,0],posTime[1,:,1],'b-')    
    
i=0
if mode == "2": #3 different zooms of the long simulation uses a loop to iterate through them
    while i < 3:
        plt.cla()
        maxacc = np.argmax(accTime)
        plt.figure(figsize=(16,9),dpi=160)
        axxy = plt.axes([0.06,0.07,0.75,0.9])
        axxz = plt.axes([0.82,0.07,0.1,0.9])
        axxz.yaxis.tick_right()
        axxz.yaxis.set_label_position("right")
        axxy.cla()
        axxz.cla()
        axxy.set_facecolor("midnightblue")
        axxz.set_facecolor("midnightblue")
        if i ==0:
            axxy.set_xlim(-150,35)
            axxy.set_ylim(-89.6,30)
            axxz.set_xlim(-15,15)
            axxz.set_ylim(-800,800)
        elif i ==1:
            axxy.set_xlim(-9.28,9.28)
            axxy.set_ylim(-6,6)
            axxz.set_xlim(-15,15)
            axxz.set_ylim(-800,800)
        elif i ==2:
            axxy.set_xlim(2,3)
            axxy.set_ylim(-0.32325,0.32325)
            axxz.set_xlim(-15,15)
            axxz.set_ylim(-800,800)
        axxy.set_xlabel("x/AU",fontsize=16)
        axxy.set_ylabel("y/AU",fontsize=16)
        axxz.set_xlabel("x/AU",fontsize=16)
        axxz.set_ylabel("z/AU",fontsize=16)
        axxy.set_title("Paths in xy plane",fontsize=16)
        axxz.set_title("Paths in xz plane",fontsize=16)
        if i <2:
            axxy.plot(0,0,'ro',label="Sun")
            axxz.plot(0,0,'ro')
            axxy.plot(relPosTime[1,:,0],relPosTime[1,:,1],color="slategrey",label="Mercury")
            axxz.plot(relPosTime[1,:,0],relPosTime[1,:,2],color="slategrey",label="Mercury")
            axxy.plot(relPosTime[2,:,0],relPosTime[2,:,1],color="darkorange",label="Venus")
            axxz.plot(relPosTime[2,:,0],relPosTime[2,:,2],color="darkorange",label="Venus")
            axxy.plot(relPosTime[3,:,0],relPosTime[3,:,1],color="turquoise",label="Earth")
            axxz.plot(relPosTime[3,:,0],relPosTime[3,:,2],color="turquoise",label="Earth")
            axxy.plot(relPosTime[4,:,0],relPosTime[4,:,1],color="red",label="Mars")
            axxz.plot(relPosTime[4,:,0],relPosTime[4,:,2],color="red",label="Mars") 
        if i < 3:
            axxy.plot(relPosTime[5,:,0],relPosTime[5,:,1],color="chocolate",label="Jupiter")
            axxz.plot(relPosTime[5,:,0],relPosTime[5,:,2],color="chocolate",label="Jupiter")
            axxy.plot(relPosTime[9,:,0],relPosTime[9,:,1],color="deeppink",label="Interloper")
            axxz.plot(relPosTime[9,:,0],relPosTime[9,:,2],color="deeppink",label="Interloper")
            axxy.plot(relPosTime[-1,maxacc,0],relPosTime[-1,maxacc,1],'rx',label="Position of closest approach")
            axxz.plot(relPosTime[-1,maxacc,0],relPosTime[-1,maxacc,2],'rx')
        if i <2:
            axxy.plot(relPosTime[6,:,0],relPosTime[6,:,1],color="gold",label="Saturn")
            axxz.plot(relPosTime[6,:,0],relPosTime[6,:,2],color="gold",label="Saturn")
        if i <1:
            axxy.plot(relPosTime[7,:,0],relPosTime[7,:,1],color="skyblue",label="Uranus")
            axxz.plot(relPosTime[7,:,0],relPosTime[7,:,2],color="skyblue",label="Uranus")
            axxy.plot(relPosTime[8,:,0],relPosTime[8,:,1],color="royalblue",label="Neptune")
            axxz.plot(relPosTime[8,:,0],relPosTime[8,:,2],color="royalblue",label="Neptune")

        i=i+1
        axxy.legend(loc="lower right",fontsize=14)
        plt.show()
        input("press enter to show next plot")

if mode == "1":   #plotting for fast run mode
    maxacc = np.argmax(accTime)
    plt.figure(figsize=(16,9),dpi=160)
    axxy = plt.axes([0.06,0.07,0.75,0.9])
    axxz = plt.axes([0.82,0.07,0.1,0.9])
    axxy.set_facecolor("midnightblue")
    axxz.set_facecolor("midnightblue")
    axxz.yaxis.tick_right()
    axxz.yaxis.set_label_position("right")
    axxy.set_xlim(-9.318,9.318)
    axxz.yaxis.tick_right()
    axxy.set_ylim(-6,6)
    axxz.set_xlim(-6,6)
    axxz.set_ylim(-0.2,0.2)
    axxy.set_xlabel("x/AU",fontsize=16)
    axxy.set_ylabel("y/AU",fontsize=16)
    axxz.set_xlabel("x/AU",fontsize=16)
    axxz.set_ylabel("z/AU",fontsize=16)
    axxy.set_title("Paths in xy plane",fontsize=16)
    axxz.set_title("Paths in xz plane",fontsize=16)
    axxy.plot(0,0,'ro',label="Sun")
    axxz.plot(0,0,'ro')
    axxy.plot(relPosTime[1,:,0],relPosTime[1,:,1],color="chocolate",label="Jupiter")
    axxz.plot(relPosTime[1,:,0],relPosTime[1,:,2],color="chocolate",label="Jupiter")
    axxy.plot(relPosTime[2,:,0],relPosTime[2,:,1],color="deeppink",label="Interloper")
    axxz.plot(relPosTime[2,:,0],relPosTime[2,:,2],color="deeppink",label="Interloper")
    axxy.plot(relPosTime[-1,maxacc,0],relPosTime[-1,maxacc,1],'rx',label="Position of closest approach")
    axxz.plot(relPosTime[-1,maxacc,0],relPosTime[-1,maxacc,2],'rx')
    axxy.legend(loc="lower right",fontsize=14)
    i=3
    plt.show()
    input("Press enter to show next graph")
    
plt.show()
#plotting for relative energy
plt.cla()
plt.figure(figsize=(16,9),dpi=160)
Eax = plt.axes([0.07,0.07,0.9,0.87])
Eax.set_xlabel("timestep/0.05*days",fontsize=16)
Eax.set_ylabel("$E/E_i -1$",fontsize=16)
Eax.set_title("Relative energy over time",fontsize=18)
Eax.plot(range(timestep),energyTime/energyTime[0]-1)
plt.show()

#The following is deprocated code that I would of liked to adapt to render an mp4 animation of the simulation had I enough time
"""
maxx = np.max(np.abs(relPosTime[-2,:,0]))
maxy = np.max(np.abs(relPosTime[-2,:,1]))
maxz = np.max(np.abs(relPosTime[-2,:,2]))
maxxy = np.min([maxy,maxx])
fig = plt.figimage
ax = plt.axes(projection='3d')

ax.set_xlim(-maxxy,maxxy)
ax.set_ylim(-maxxy,maxxy)
ax.set_zlim(-maxz,maxz)

ax.plot3D(relPosTime[1,:,0],relPosTime[1,:,1],relPosTime[1,:,2],color="slategrey")
ax.plot3D(relPosTime[2,:,0],relPosTime[2,:,1],relPosTime[2,:,2],color="darkorange")
if mode == "2":
    ax.plot3D(relPosTime[3,:,0],relPosTime[3,:,1],relPosTime[3,:,2],color="turquoise")
    ax.plot3D(relPosTime[4,:,0],relPosTime[4,:,1],relPosTime[4,:,2],color="red")
    ax.plot3D(relPosTime[5,:,0],relPosTime[5,:,1],relPosTime[5,:,2],color="chocolate")
    ax.plot3D(relPosTime[6,:,0],relPosTime[6,:,1],relPosTime[6,:,2],color="gold")
    ax.plot3D(relPosTime[7,:,0],relPosTime[7,:,1],relPosTime[7,:,2],color="skyblue")
    ax.plot3D(relPosTime[8,:,0],relPosTime[8,:,1],relPosTime[8,:,2],color="darkblue")
    ax.plot3D(relPosTime[9,:,0],relPosTime[9,:,1],relPosTime[9,:,2],color="deeppink")
    ax.plot3D(0,0,0,'ro')
else:
    fig = plt.figimage
    ax = plt.axes(projection='3d')
    ax.set_zlim(-50,50)
    ax.plot3D(relPosTime[0,:,0],relPosTime[0,:,1],relPosTime[0,:,2],'r-')
    ax.plot3D(relPosTime[5,:,0],relPosTime[5,:,1],relPosTime[5,:,2],'b-')
    ax.plot3D(relPosTime[-1,:,0],relPosTime[-1,:,1],relPosTime[-1,:,2],'g-')


####2D animation

axanim = plt.axes()
animstep = 500
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
    axesLimit = np.abs(np.max(relPosTime[0:bodies+4,:,:]))*1.1
    axanim.set_xlim(-axesLimit,axesLimit)
    axanim.set_ylim(-axesLimit,axesLimit)
    for i in range(bodies+1):  
        axanim.plot(relPosTime[i+4,j:j2,0],relPosTime[i+4,j:j2,1],'-')
        axanim.plot(relPosTime[i+4,j2,0],relPosTime[i+4,j2,1],'o')
    plt.pause(0.1)
"""