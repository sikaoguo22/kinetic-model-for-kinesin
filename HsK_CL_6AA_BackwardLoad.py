# -*- coding: utf-8 -*-
"""

Simulating the movement of Kinesin at different concentrations of ATP（MC Simulation）

"""

#Import the required packages
import numpy as np

#Set random number seed
myseed = 91

#The state of kinesin is represented by K
K = 113 
#K is a three-digit number, and first digit represents the state of the trailing head.（1-ATP;2-ADP;3-free）
#Second digit the status of Kinesin（1-two head bound;2-Intermediate state）
#Third digit represents the state of the leading head（1-ATP;2-ADP;3-free）

#Set parameter value
kplus = 107 #trailing head ATP->ADP rate；s-1
kminus = 4 #leading head ATP->ADP rate；s-1
kb = 3.3 #ATP bind second order rate；s-1uM-1
kb_slow = 1 #leading head ATP bind slow times
ATP_conc = 10 #ATP concentrations；uM
F = 7 #load；pN
F_all = [0,1,2] #all loads；pN
r0 = 600 #0pN step ratio
Fs = 5 #Stop Force;pN
r = r0 ** (1-F/Fs) #step ratio
Pe = r / (kplus/kminus + r) #Forward probability
P0 = 0.33 #0pN, The escape probability after ATP hydrolysis
deltad = 1 #nm, distance parameter
kBT = 4.11#pN.nm, kBT
P0F = 1 - (1 - P0)*np.exp(-F*deltad/kBT)#The escape probability of leading head after ATP hydrolysis when F>0
d = 8.2 #step length;nm
kminus2 = 30 #leading head ATP release rate; s-1
tdwell = 0.001 #Intermediate state judgment, s

#The time of kinesin movement is represented by T
T = 0.0

#The position of trailing head is represented by X
X = 0
P = 0.0#position of Kinesin

#Number of steps
Nstep = 0

#The number of ATP consumed
NofATP = 0

#Number of simulation experiments N
N = 1000

#Record the velocity of each simulation V
V = []

#Record the number of steps per simulation
NATP = []

#Record randomness
Randomness = []

#steps per simulation
Number = []

#Simulation start

for F_index in range(len(F_all)):

    #parameters initialize
    F = F_all[F_index]
    r = r0 ** (1-F/Fs) #step ratio
    Pe = r / (kplus/kminus + r) #Forward probability
    P0F = 1 - (1 - P0)*np.exp(-F*deltad/kBT)#The escape probability of leading head after ATP hydrolysis when F>0
    print("Force: ", F, "pN")
    print("[ATP]: ", ATP_conc, "uM")
    V = []
    NATP = []
    Randomness = []
    Number = []
    Ratio = []
    Frecord = []
    Brecord = []
    Fstep = 0
    Bstep = 0
    np.random.seed(myseed)
    
    for index in range(N):
        K = 113 #Initial state：T--0
        T = 0.0
        X = 0
        Time = []
        Position = []
        NofATP = 0
        Nstep = 0
        Fstep = 0
        Bstep = 0
        Time.append(T)
        P = X+d/2
        Position.append(P)
        while 1:        
            K_now = K #Record current status
            if K_now == 113: #T--0
                kT = kplus + kb*ATP_conc
                T = T + (np.random.exponential(1.0/kT))
                decision_ch = np.random.random() #choice
                if decision_ch <= kplus / kT: #trailing head hydr
                    if np.random.random() <= P0:
                        K = 3221 #0D;last digit 1 represents forward
                        X = X+d
                        P = X
                        NofATP = NofATP+1
                    else:
                        K = 313 #not escape after hydr
                        X = X
                        P = X
                        NofATP = NofATP+1    
                else: #leading head bind ATP
                    K = 111 #T--T
                    X = X
                    P = X+d/2
                Time.append(T)
                Position.append(P)
                    
            if K_now == 3221: #0D;last digit 1 represents forward
                kT = kb*ATP_conc
                T = T + (np.random.exponential(1.0/kT))
                decision_ch = np.random.random() #choice
                if decision_ch <= Pe: #Forward
                    K = 113 #T--0
                    X = X
                    P = X+d/2
                    Nstep = Nstep+1
                    Fstep = Fstep+1
                else: #Back
                    K = 311 #0--T
                    X = X-d
                    P = X+d/2
                Time.append(T)
                Position.append(P)
                    
            if K_now == 3222: #0D;last digit 2 represents backward
                kT = kb*ATP_conc
                dwellmid = (np.random.exponential(1.0/kT))
                T = T + dwellmid
                decision_ch = np.random.random() #choice
                if dwellmid > tdwell:
                    Bstep = Bstep+1
                if decision_ch <= Pe: #Forward
                    K = 113 #T--0
                    X = X
                    P = X+d/2
                    if dwellmid > tdwell:
                        Fstep = Fstep+1
                else: #Back
                    K = 311 #0--T
                    X = X-d
                    P = X+d/2
                    Nstep = Nstep+1
                Time.append(T)
                Position.append(P)
                    
            if K_now == 111: #T--T
                kT = kplus + kminus + kminus2
                T = T + (np.random.exponential(1.0/kT))
                decision_ch = np.random.random() #choice
                if decision_ch <= kplus / kT: #trailing head hydr
                    if np.random.random() <= P0:
                        K = 1221 #TD;last digit 1 represents forward
                        X = X+d
                        P = X
                        NofATP = NofATP+1
                    else:
                        K = 311#not escape after hydr
                        X = X
                        P = X
                        NofATP = NofATP+1
                elif decision_ch <= (kplus+kminus) / kT: #leading head hydr
                    if np.random.random() <= P0F:
                        K = 1222 #TD;last digit 2 represents backward
                        X = X
                        P = X
                        NofATP = NofATP+1
                    else:
                        K = 113#not escape after hydr
                        X = X
                        P = X
                        NofATP = NofATP+1    
                else : #leading head release ATP
                    K = 113 #T--0
                    X = X
                    P = X+d/2
                Time.append(T)
                Position.append(P)
                    
            if K_now == 311: #0--T
                kT = kminus + kb*ATP_conc + kminus2
                T = T + (np.random.exponential(1.0/kT))
                decision_ch = np.random.random() #choice
                if decision_ch <= kminus / kT: #leading head hydr
                    if np.random.random() <= P0F:
                        K = 3222 #0D;last digit 2 represents backward
                        X = X
                        P = X
                        NofATP = NofATP+1
                    else:
                        K = 313#not escape after hydr
                        X = X
                        P = X
                        NofATP = NofATP+1
                elif decision_ch <= (kminus+kb*ATP_conc) / kT: #trailing head bind ATP
                    K = 111 #T--T
                    X = X
                    P = X+d/2
                else : #leading head release ATP
                    K = 313 #0--0
                    X = X
                    P = X+d/2
                Time.append(T)
                Position.append(P)
                    
            if K_now == 1221: #TD
                decision_ch = np.random.random() #choice
                if decision_ch <= Pe: #Forward
                    K = 113 #T--0
                    X = X
                    P = X+d/2
                    Nstep = Nstep+1
                    Fstep = Fstep+1
                else: #Back
                    K = 311 #0--T
                    X = X-d
                    P = X+d/2
                Time.append(T)
                Position.append(P)
                    
            if K_now == 1222: #TD
                decision_ch = np.random.random() #choice
                if decision_ch <= Pe: #Forward
                    K = 113 #T--0
                    X = X
                    P = X+d/2
                else: #Back
                    K = 311 #0--T
                    X = X-d
                    P = X+d/2
                    Nstep = Nstep+1
                    Bstep = Bstep+1
                Time.append(T)
                Position.append(P)
                    
            if K_now == 313: #0--0
                kT = kb*ATP_conc + kb*ATP_conc
                T = T + (np.random.exponential(1.0/kT))
                decision_ch = np.random.random() #choice
                if decision_ch <= kb*ATP_conc / kT: #leading head bind ATP
                    K = 311 #0--T
                    X = X
                    P = X+d/2
                else : #trailing head bind ATP
                    K = 113 #T--0
                    X = X
                    P = X+d/2
                                        
            if T > 1: #duration time per simulation,s
                break
                
        V.append(P/T)
        NATP.append(P)
        if index % 200 == 199:
            Randomness.append(1.0/(np.var(NATP)/d/np.mean(NATP)))
            NATP = []
        Frecord.append(Fstep)
        Brecord.append(Bstep)
            
        if Nstep != 0:
            Number.append(NofATP/Nstep)
        if Bstep != 0:
            Ratio.append(Fstep/Bstep)

    print("Velocity: ",np.mean(V),"+-",np.std(V)/(len(V)**0.5),"nm/s,SEM")
    print("Randomness: ",np.mean(Randomness),"+-",np.std(Randomness)/(len(Randomness)**0.5),",SEM")
    print("ATPs per step: ",np.mean(Number),"+-",np.std(Number)/(len(Number)**0.5),",SEM")
    print("Forward step/Backward step: ",np.mean(Ratio),"+-",np.std(Ratio)/(len(Ratio)**0.5),",SEM")
    print("-------------------------------------------------------------------")
