# ====================================================================
# Author 				: swc21
# Date 					: 2018-03-14 09:47:35
# Project 				: ClusterFiles
# File Name 			: Solersystem
# Last Modified by 		: swc21
# Last Modified time 	: 2018-03-14 10:28:50
# ====================================================================
#
# Filename: SolerSystem.py
import numpy as np

from datetime import datetime


def pps(data):
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    print("PARTICLE LIST-------------------------------------------------------------")
    for i in range(len(data)):
        a = data[i]
        print(i, "X,Y,Z:", a[0], "Vx,Vy,Vz:", a[1],
              "Ax,Ay,Az:", a[2], "Mass:", a[3])
    print("END OF LIST---------------------------------------------------------------")
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    return
#[Bring Mass, xyz position, Velocity and Acceleration together in a list]


def add_particles(N, alpha, M_min, M_max, R, V_avg, A_avg):
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
                Masses.append([M])
        #--Return populated list of masses
        return Masses
    #[get (N)number of random sets of x,y,z coords as a function of (R)radius]

    def get_random_xyz(N, R):
        #<imports>
        import numpy as np
        from numpy import pi
        #<Continuous uniform distribution as function of radius(R)>
        # Spherical coods:
        #--Phi | [0, 2*pi]
        phi = (2*pi-0)*np.random.random_sample((N,))+0
        #--Theta | [0, pi]
        costh = (2)*np.random.random_sample((N,))-1
        theta = np.arccos(costh)
        #--U unit radius r**3 | [0, 1]
        u = np.random.random_sample((N,))
        #<Convert to x,y,z coords>
        r = R * np.power(u, 1/3)
        x = r * np.sin(theta)*np.cos(phi)
        y = r * np.sin(theta)*np.sin(phi)
        z = r * np.cos(theta)
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

    def get_random_VxVyVz(N, V_avg):
        #<imports>
        import random
        #<Output array for vx,vy,vz velocity>
        #--Velocity output list
        Velocity_xyz = []
        #--Build the list
        for i in range(0, N):
            #--Multiply the array by avarage velocity paramater V_avg
            vx = V_avg*random.uniform(-1, 1)
            vy = V_avg*random.uniform(-1, 1)
            vz = V_avg*random.uniform(-1, 1)
            #--Append velocity for particle[i/N] to the list
            Velocity_xyz.append([vx, vy, vz])
        #--Return populated list of velocities
        return Velocity_xyz
    #[get (N)number of random sets of Ax,Ay,Az acceleration based on (A_avg)avarage acceleration]

    def get_random_AxAyAz(N, A_avg):
        #<imports>
        import random
        #<Output array for ax,ay,az acceleration>
        #--Acceleration output list
        Acceleration_xyz = []
        #--Build the list
        for i in range(0, N):
            #--Multiply the array by avarage acc paramater A_avg
            ax = A_avg*random.uniform(-1, 1)
            ay = A_avg*random.uniform(-1, 1)
            az = A_avg*random.uniform(-1, 1)
            #--Append acceleration for particle[i/N] to the list
            Acceleration_xyz.append([ax, ay, az])
        #--Return populated list of velocities
        return Acceleration_xyz
    start_time = datetime.now()
    #<CALL Initial Mass Function>
    Masses = sampleFromSalpeter(N, alpha, M_min, M_max)
    #<Continuous uniform distribution as function of radius(R)>
    Positions_xyz = get_random_xyz(N, R)
    #<Avarage Velocity Distribution>
    Velocity_xyz = get_random_VxVyVz(N, V_avg)
    #<Avarage Acceleration Distribution>
    Acceleration_xyz = get_random_AxAyAz(N, A_avg)
    List = []
    #--Build the data-set
    for i in range(0, N):
        #--Variables for ID, Mass, Position, Velocity and Acceleration per particle
        partmas = Masses[i]
        partpos = Positions_xyz[i]
        partvel = Velocity_xyz[i]
        partacc = Acceleration_xyz[i]
        r = np.ones((N, 3))
        v = np.ones((N, 3))
        a = np.ones((N, 3))
        m = np.ones((N, 1))
        r[i] = [partpos[0], partpos[1], partpos[2]]
        v[i] = [partvel[0], partvel[1], partvel[2]]
        a[i] = [partacc[0], partacc[1], partacc[2]]
        m[i] = [partmas[0]]
        List.append([r[i], v[i], a[i], m[i]])
    end_time = datetime.now()
    print('Duration of particle creation: {}'.format(end_time - start_time))
    pps(List)
    return List


def last_gravity_calculation(b, a_vector):
    x1, y1, z1 = a[0]
    x2, y2, z2 = b[0]
    distance = np.sqrt((x2-x1)**2+(y2-y1)**2(z2-z1)**2)
    a1 = b
    PRINT("last_gravity_calculation_count finished loop",
          LAST_GRAVITY_CALCULATION_COUNT)
    return -np.concatenate([acceleration_on_a[0], acceleration_on_a[1], acceleration_on_a[2]])


def acceleration_xyz(a, Data):
    for j in range(N):
        a = Data[j]
        x1, y1, z1 = a[0]
        a_pos = x1**2 + y1**2 + z1**2
        a_vector = np.sqrt(a_pos)
        for i in range(N):
            b = Data[i]
            x2, y2, z2 = b[0]
            if x1+y1+z1 != x2+y2+z2:
                acceleration += last_gravity_calculation(b, a_vector)
            print("acceleration_xyz_outside_loop finished loop")


def update_particles(N, dt, Data, number):
    start_update_particles = datetime.now()
    print("updating particles for the", number, "time. N=", N,
          ". dt=", dt, ". Data is", len(Data), "lines long.")
    acceleration_updates = acceleration_xyz(Data, N)
    updated_particles = []
    for i in range(0, N):
        particle = Data[i]
        x, y, z = particle[0]
        vx, vy, vz = particle[1]
        m = particle[3]
        for j in range(0, N):
            try:
                ax2, axy2, az2 = acceleration_updates[j]
                vx += ax2*dt/m
                vy += ay2*dt/m
                vz += az2*dt/m
                particle[1] = [vx, vy, vz]
                x += vx*dt
                y += vy*dt
                z += vz*dt
                particle[0] = [x, y, z]
                particles.append(particle)
            except Exception:
                pass  # print(acceleration_updates[j])
    return updated_particles


def step(T, dt, Data):
    from datetime import datetime
    print("[STEPPING", T, "STEPS: dt=", dt, ":", len(
        Data), "particles] : {}".format(datetime.now()))
    start_step = datetime.now()
    N = len(Data)
    i = 1
    while i <= T:
        print("[STARTING STEP", i, "NOW] : {}".format(datetime.now()))
        start_timestep = datetime.now()
        number = i
        print("------->calling [update_particles]")
        updates = update_particles(N, dt, Data, number)
        print("------->starting [update loop]")
        start_update_loop = datetime.now()
        for j in range(0, N):
            try:
                Data[j] = updates[j]
            except Exception:
                print(len(Data), len(updates))
        end_update_loop = datetime.now()
        print('------->update loop finished. Duration: {}'.format(end_update_loop - start_update_loop))
        end_timestep = datetime.now()
        print("------->[Time Step %d finished] - " %
              i, 'Duration: {}'.format(end_timestep - start_timestep))
        i += 1
    end_step = datetime.now()
    print('completed all steps. Duration: {}'.format(end_step - start_step))
