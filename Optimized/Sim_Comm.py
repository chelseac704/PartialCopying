import SLNK as slnk
import numpy as np
import sys

trials = 2500
psize = 100
reps = 50
n = 15
solutions = 2**n

SHARE_RATE_PARTIAL = 0.5 # Partial share
SHARE_RATE_FULL = 1.0 # Full share
SHARE_RADIUS_GLOBAL = psize - 1 # Global
SHARE_RADIUS_LOCAL = 1 # Local
NO_MUT_RATE = 0.0 # Regular, no copying mutation
MUT_RATE = 0.5    # New, half genes mutated after copying

### Community experiment
ncomm = 5
conditions = 5
ig = np.array([0.2,0.8,0.95,0.9875,1.0])
og = (1-ig)/(ncomm-1)
sh = [0.5,1.0]

class SLSim:
    def __init__(self):
        pass

    def simulateReg(self,n,k,reps):
        '''For a single k, run the simulation many times for each condition'''
        avg = np.zeros((conditions,reps))
        for rep in range(reps):
            #print("\tR: ",rep)
            nk = slnk.NKLandscape(n,k)
            pop = slnk.Population(psize,n,nk,False)
            initial_genotypes = pop.genotypes.copy()
            for condition in range(conditions):
                #print("\t\tC: ",condition)
                pop.set_pop(initial_genotypes)
                pop.share_rate = exp_condition[condition][0]
                pop.share_radius = exp_condition[condition][1]
                pop.mut_rate = exp_condition[condition][2]
                for trial_num in range(1,trials+1):
                    pop.share()
                    pop.learn()
                avg[condition,rep] = pop.avgfit()
        return avg

    def simulateComm(self,n,k,reps):
        '''For a single k, run the simulation many times for each condition'''
        avg = np.zeros((int(len(ig)*len(sh)),reps))
        for rep in range(reps):
            print("\tR: ",rep)
            nk = slnk.NKLandscape(n,k)
            pop = slnk.Population(psize,n,nk,True)
            initial_genotypes = pop.genotypes.copy()
            condition = 0
            for comm in range(len(ig)):
                pop.set_community(ncomm, ig[comm], og[comm])
                for share in range(len(sh)):
                    pop.set_pop(initial_genotypes)
                    pop.share_rate = sh[share]
                    print("\t\tC: ", condition, ig[comm], og[comm], sh[share])
                    for trial_num in range(1,trials+1):
                        pop.share()
                        pop.learn()
                    avg[condition,rep] = pop.avgfit()
                    condition += 1
        return avg

    def exp(self,minK,maxK,step,id,community=False):
        '''Run the experiment for all values of K and all conditions and save the data'''
        for k in range(minK,maxK+1,step):
            print("K: ",k)
            if community:
                avg=self.simulateComm(n,k,reps)
            else:
                avg=self.simulateReg(n,k,reps)
            np.save("avg_"+str(k)+"_"+str(id)+".npy", avg)

id = int(sys.argv[1])
sim = SLSim()
sim.exp(0,14,1,id,True)
