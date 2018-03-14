# ====================================================================
# Author                : swc21
# Date                  : 2018-03-14 09:42:27
# Project               : ClusterFiles
# File Name             : Solersystem2
# Last Modified by      : swc21
# Last Modified time    : 2018-03-14 10:57:43
# ====================================================================
# 
# ====================================================================
# Author 				: swc21
# Date 					: 2018-03-14 09:42:27
# Project 				: ClusterFiles
# File Name 			: Solersystem2
# Last Modified by 		: swc21
# Last Modified time 	: 2018-03-14 10:28:56
# ====================================================================
#
# Filename: SolerSystem.py
import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import sys
N_Stinkers = 3
DT = 3e-2*3.154e7  # Time year in s
C = 2.99*1e8  # m/s
C_squared = C**2  # (m/s)**2
G = 6.67e-11  # N*(m**2)/kg**2_
m_e = 9.109e-31  # kg_mass of electron
m_p = 1.6e-27  # kg_protonmass
u_0 = 4.0e-7*np.pi  # N/amp**2_permeability
E_0 = 1.0/(C_squared*u_0)  # Farad/m_permitivity
K_e = 1.0/4*np.pi*E_0  # N*(m**2)/Coulombs**2_Coulombsconstant
scale_m = 1e-1*1.988e30  # Mass kg
scale_x = 3e16  # Length m
scale_v = 4e2  # m/s_Velocity
scale_q = 1.6e-19  # protoncharge_Coulombs
#[Bring Mass, xyz position, Velocity and ID together in a 2D array]


def add_particles(N):
    alpha = 2.18
    M_min = 1.0
    M_max = 10.0
    V_avg = 1.0

    def coin_flip(a, b):
        if np.random.rand(1) >= 0.5:
            return a
        else:
            return b
    #[Initial Mass Function for (N)stars]

    def sampleFromSalpeter(N, alpha, M_min, M_max):
        #<import>
        import random as random
        import math
        #<Salpeter power law function>
        #--Random seed
        random.seed(13)
        #--Convert limits to LogM
        log_M_Min = math.log(M_min)
        log_M_Max = math.log(M_max)
        #--Likelihood @ M_min
        maxlik = math.pow(M_min, 1.0 - alpha)
        #<Output array for particle mass>
        #--Mass output list
        Masses = []
        #--Build the list
        while (len(Masses) < N):
            # for i in range(0,N):
            #--Draw candidate from logM interval
            logM = random.uniform(log_M_Min, log_M_Max)
            M = math.exp(logM)
            #--Compute likelihood of candidate from Salpeter SMF.
            likelihood = math.pow(M, 1.0 - alpha)
            #--Accept randomly.
            u = random.uniform(0.0, maxlik)
            if (u < likelihood):
                #--Append mass for particle[i/N] to the list
                mass = M*scale_m
                Masses.append([mass])
        #--Return populated list of masses
        for i in range(N_Stinkers):
            if not i:
                Masses[i][0] = Masses[i][0]*30
            else:
                Masses[i][0] = Masses[i][0]*10
        return Masses
    #[get (N)number of random sets of x,y,z coords as a function of (R)radius]

    def get_angular_xyz_vxvyvz(N):
        phi = (2*np.pi-0)*np.random.random_sample((N,))+0
        costh = (2)*np.random.random_sample((N,))-1
        theta = np.arccos(costh)
        u = np.random.random_sample((N,))
        #<Convert to x,y,z coords>
        r = np.power(u, 1/3)*(1.0-np.random.rand(1)[0]**3)*scale_x
        x = r*np.sin(theta)*np.cos(phi)
        y = r*np.sin(theta)*np.sin(phi)
        z = r*np.cos(theta)
        #<Output array for x,y,z position>
        #--Position output list
        Velocity_xyz = []
        Positions_xyz = []
        #--Build the list
        for i in range(0, N):
            #--Append position for particle[i/N] to the list
            Positions_xyz.append([x[i], y[i], z[i]])
            Velocity_xyz.append(
                [-np.sign(y[i])*scale_v*1.5, np.sign(x[i])*scale_v*1.5, -np.sign(z[i])*scale_v*2])
        #--Return populated list of positions
        for i in range(N_Stinkers):
            if not i:
                Positions_xyz[i] = [0., 0., 0., ]
                Velocity_xyz[i] = [0., 0., 0., ]
        return Positions_xyz, Velocity_xyz
    #[get (N)number of random sets of x,y,z coords as a function of (R)radius]

    def get_random_xyz(N):
        phi = (2*np.pi-0)*np.random.random_sample((N,))+0
        costh = (2)*np.random.random_sample((N,))-1
        theta = np.arccos(costh)
        u = np.random.random_sample((N,))
        #<Convert to x,y,z coords>
        r = np.power(u, 1/3)*(1.0-np.random.rand(1)[0]**2)*scale_x
        x = r*np.sin(theta)*np.cos(phi)
        y = r*np.sin(theta)*np.sin(phi)
        z = r*np.cos(theta)
        #<Output array for x,y,z position>
        #--Position output list
        Positions_xyz = []
        #--Build the list
        for i in range(0, N):
            #--Append position for particle[i/N] to the list
            Positions_xyz.append([x[i], y[i], z[i]])
        #--Return populated list of positions
        return Positions_xyz
    #[get (N)number of random sets of Vx,Vy,Vz velocity based on (V_avg)avarage velocity]

    def get_random_VxVyVz(N):
        #<imports>
        import random
        #<Output array for vx,vy,vz velocity>
        #--Velocity output list
        Velocity_xyz = []
        #--Build the list
        for i in range(0, N):
            #--Multiply the array by avarage velocity paramater V_avg
            vx = scale_v*random.uniform(-1, 1)
            vy = scale_v*random.uniform(-1, 1)
            vz = scale_v*random.uniform(-1, 1)
            #--Append velocity for particle[i/N] to the list
            Velocity_xyz.append([vx, vy, vz])
        #--Return populated list of velocities
        return Velocity_xyz
    #<CALL Initial Mass Function>
    Masses = sampleFromSalpeter(N, alpha, M_min, M_max)
    #<Continuous uniform distribution as function of radius(R)>
    #Positions_xyz = get_random_xyz(N)
    #<Avarage Velocity Distribution>
    #Velocity_xyz = get_random_VxVyVz(N)
    Positions_xyz, Velocity_xyz = get_angular_xyz_vxvyvz(N)
    #Accelerations = [0.,0.,0.]*(N+1)
    # Positions_xyz.append([0.,0.,0.])
    # Velocity_xyz.append([0.,0.,0.])
    # Masses.append([(scale_m*100)])
    # return zip(Positions_xyz,Velocity_xyz,Accelerations,Masses)
    #List = []
    r = np.ones((N, 5))
    v = np.ones((N, 3))
    for i in range(0, N):
        #--Variables for Mass, Position, Velocity and Charge per particle
        partmass = Masses[i]
        partpos = Positions_xyz[i]
        partvel = Velocity_xyz[i]
        partcharge = [scale_q*partmass[0]/(coin_flip(m_p, -m_e))]
        r[i] = np.asarray(
            [partpos[0], partpos[1], partpos[2], partmass[0], partcharge[0]])
        v[i] = np.asarray([partvel[0], partvel[1], partvel[2]])
        #a = [0., 0., 0.]
    return zip(r, v)


'''
def energy(Data):
    KE_tot = 0.0
    PE_tot = 0.0
    for particle_a in Data:
        PE_a = 0.0
        for particle_b in Data:
            if particle_a is not particle_b:
                x1,y1,z1 = particle_a[0]
                x2,y2,z2 = particle_b[0]
                vec = [(x2-x1),(y2-y1),(z2-z1)]
                magnitude = np.sqrt(np.square(vec).sum())
                PE_a += particle_b[3][0]/magnitude
        PE_a *= particle_a[3][0]
        PE_tot += PE_a
        KE_tot += particle_a[3][0]*(particle_a[1][0]**2+particle_a[1][1]**2+particle_a[1][2]**2)
    KE_tot /= 2.0
    PE_tot *= G/-2.0
    return PE_tot,KE_tot
    def plot(Data,T):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection='3d')
    plt.gca().patch.set_facecolor('white')
    masses = np.asarray([size[3] for size in Data])
    max_mass = masses[N_Stinkers:].max()
    masses /= max_mass
    #velocity = np.asarray([size[1] for size in Data])
    mymap = plt.cm.cool
    for i in range(len(Data)):
        particles = Data[i]
        #speeds = np.sqrt(np.square(particles[1]).sum())/scale_v*10
        x, y, z = particles[0] 
        if 0.0 < masses[i] <= 1.0: 
            ax.scatter(x,y,z,s=masses[i]**2*100, c=z/scale_x, cmap=mymap, alpha=.99)
        else:
            ax.scatter(x,y,z,s=550, c='r', alpha=.2)
            #print('STINKER ALERT')
        lim = scale_x/2.0
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])
        ax.set_zlim([-lim,lim])
        #ax.set_axis_off()
    fig.savefig('/home/sol/Documents/code/SolerSysytem/plots/plot_%d.png' % T)
    plt.close()
    #print('ploted',T)
def gamma(vector):
    total = vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]
    return 1.0/np.sqrt(1.0+total/C_squared)
def sr_helper(a_vec,v_vec):
    coeff = gamma(v_vec)
    return [a_vec[0]+coeff*v_vec[0], a_vec[1]+coeff*v_vec[1], a_vec[2]+coeff*v_vec[2]]
def update_all(Data):
    N = len(Data)
    for particle_a in Data:
        acceleration_on_ax_update = 0.0
        acceleration_on_ay_update = 0.0
        acceleration_on_az_update = 0.0
        for particle_b in Data:
            if particle_a[3] != particle_b[3]:
                x1,y1,z1 = particle_a[0]
                x2,y2,z2 = particle_b[0]
                vec = [(x2-x1),(y2-y1),(z2-z1)]
                magnitude = np.sqrt(np.square(vec).sum())
                if magnitude < scale_x*1e-5:
                    print('Close encounter')
                    continue
                acceleration_on_ax_update += particle_b[3]*(x2-x1)/magnitude**3
                acceleration_on_ay_update += particle_b[3]*(y2-y1)/magnitude**3
                acceleration_on_az_update += particle_b[3]*(z2-z1)/magnitude**3
        impulse = [(particle_a[2][0]+acceleration_on_ax_update)*G*DT/2.0, (particle_a[2][1]+acceleration_on_ay_update)*G*DT/2.0, (particle_a[2][2]+acceleration_on_az_update)*G*DT/2.0]
        new_velocity = sr_helper([0.0,0.0,0.0], sr_helper(impulse, particle_a[1]))
        particle_a[1][0] = new_velocity[0]
        particle_a[1][1] = new_velocity[1]
        particle_a[1][2] = new_velocity[2]
        if gamma(particle_a[1])>1.1:
            print('they gone plaid')
        #particle_a[1][0] += (particle_a[2][0]+acceleration_on_ax_update)*G*DT/2.0
        #particle_a[1][1] += (particle_a[2][1]+acceleration_on_ay_update)*G*DT/2.0
        #particle_a[1][2] += (particle_a[2][2]+acceleration_on_az_update)*G*DT/2.0
        particle_a[0][0] += particle_a[1][0]*DT
        particle_a[0][1] += particle_a[1][1]*DT
        particle_a[0][2] += particle_a[1][2]*DT
        particle_a[2][0] = acceleration_on_ax_update
        particle_a[2][1] = acceleration_on_ay_update
        particle_a[2][2] = acceleration_on_az_update
def step_time(Data, t=5000, save=20, error=None):
    PE_int,KE_int = energy(Data)
    E_int = KE_int+2*PE_int
    print(PE_int,KE_int,E_int)
    plot(Data,1)
    if error is None:
        PE,KE = energy(Data)
        for i in range(0,t):
            if i == t//2:
                global DT
                DT *= -1.0
            update_all(Data)
            my_str = '\r'+str(i)+'----------->PE:'+str(PE/PE_int)+'  KE:'+str(KE/KE_int)+'  Sum:'+str((KE+2*PE)/(E_int))
            sys.stdout.write(' '*len(my_str))
            sys.stdout.flush()
            sys.stdout.write(my_str)
            sys.stdout.flush()
            if not i%save:
                plot(Data,i)
                PE,KE = energy(Data)
                #sys.stdout.write('\r'+str(i)+'----------->PE:'+str(PE/PE_int)+'  KE:'+str(KE/KE_int)+'  Sum:'+str((KE+2*PE)/(E_int)))
                #sys.stdout.flush()
    else:
        update_all(Data)
        PE,KE = energy(Data)
        i=0
        while 1-error<(PE/PE_int)<1+error and 1-error<(KE/KE_int)<1+error:
            update_all(Data)
            i+=1
            if not i%save:
                plot(Data,i)
                PE,KE = energy(Data)
                #print('PE:',PE,'   KE:',KE)
                sys.stdout.write('\r'+str(i)+'----------->PE:'+str(PE/PE_int)+'  KE:'+str(KE/KE_int))
                sys.stdout.flush()
        print(i)
Data = add_particles(5)
step_time(Data)
'''
