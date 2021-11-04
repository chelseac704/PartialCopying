import SLNK as slnk
import numpy as np
import sys

trials = 2500
psize = 150
reps = 50
n = 15
solutions = 2**n

SHARE_RATE_PARTIAL = 0.5 # Partial share
SHARE_RATE_FULL = 1.0 # Full share
SHARE_RADIUS_GLOBAL = psize - 1 # Global
SHARE_RADIUS_LOCAL = 1 # Local
NO_MUT_RATE = 0.0 # Regular, no copying mutation
MUT_RATE = 0.5    # New, half genes mutated after copying

# Main four conditions: Full and Partial sharing, Global and Local populations.
conditions = 4
exp_condition = np.zeros((conditions,3))
exp_condition[0] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, NO_MUT_RATE]
exp_condition[1] = [SHARE_RATE_FULL, SHARE_RADIUS_LOCAL, NO_MUT_RATE]
exp_condition[2] = [SHARE_RATE_PARTIAL, SHARE_RADIUS_GLOBAL, NO_MUT_RATE]
exp_condition[3] = [SHARE_RATE_PARTIAL, SHARE_RADIUS_LOCAL, NO_MUT_RATE]

####
# conditions = 16
# exp_condition = np.zeros((conditions,4))
# exp_condition[0] = [0.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[1] = [1.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[2] = [2.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[3] = [3.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[4] = [4.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[5] = [5.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[6] = [6.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[7] = [7.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[8] = [8.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[9] = [9.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[10] = [10.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[11] = [11.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[12] = [12.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[13] = [13.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[14] = [14.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]
# exp_condition[15] = [15.0/n, SHARE_RADIUS_GLOBAL, NO_MUT_RATE, 3]

####
# conditions = 16
# exp_condition = np.zeros((conditions,4))
# exp_condition[0] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 15.0/n, 3]
# exp_condition[1] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 14.0/n, 3]
# exp_condition[2] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 13.0/n, 3]
# exp_condition[3] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 12.0/n, 3]
# exp_condition[4] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 11.0/n, 3]
# exp_condition[5] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 10.0/n, 3]
# exp_condition[6] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 9.0/n, 3]
# exp_condition[7] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 8.0/n, 3]
# exp_condition[8] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 7.0/n, 3]
# exp_condition[9] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 6.0/n, 3]
# exp_condition[10] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 5.0/n, 3]
# exp_condition[11] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 4.0/n, 3]
# exp_condition[12] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 3.0/n, 3]
# exp_condition[13] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 2.0/n, 3]
# exp_condition[14] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 1.0/n, 3]
# exp_condition[15] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, 0.0/n, 3]

# Varying Demes
# conditions = 22
# exp_condition = np.zeros((conditions, 4))
# exp_condition[0] = [SHARE_RATE_FULL, SHARE_RADIUS_LOCAL, NO_MUT_RATE,3]
# exp_condition[1] = [SHARE_RATE_FULL, 10, NO_MUT_RATE,3]
# exp_condition[2] = [SHARE_RATE_FULL, 20, NO_MUT_RATE,3]
# exp_condition[3] = [SHARE_RATE_FULL, 30, NO_MUT_RATE,3]
# exp_condition[4] = [SHARE_RATE_FULL, 40, NO_MUT_RATE,3]
# exp_condition[5] = [SHARE_RATE_FULL, 50, NO_MUT_RATE,3]
# exp_condition[6] = [SHARE_RATE_FULL, 60, NO_MUT_RATE,3]
# exp_condition[7] = [SHARE_RATE_FULL, 70, NO_MUT_RATE,3]
# exp_condition[8] = [SHARE_RATE_FULL, 80, NO_MUT_RATE,3]
# exp_condition[9] = [SHARE_RATE_FULL, 90, NO_MUT_RATE,3]
# exp_condition[10] = [SHARE_RATE_FULL, SHARE_RADIUS_GLOBAL, NO_MUT_RATE,3]
# exp_condition[11] = [SHARE_RATE_PARTIAL, SHARE_RADIUS_LOCAL, NO_MUT_RATE,3]
# exp_condition[12] = [SHARE_RATE_PARTIAL, 10, NO_MUT_RATE,3]
# exp_condition[13] = [SHARE_RATE_PARTIAL, 20, NO_MUT_RATE,3]
# exp_condition[14] = [SHARE_RATE_PARTIAL, 30, NO_MUT_RATE,3]
# exp_condition[15] = [SHARE_RATE_PARTIAL, 40, NO_MUT_RATE,3]
# exp_condition[16] = [SHARE_RATE_PARTIAL, 50, NO_MUT_RATE,3]
# exp_condition[17] = [SHARE_RATE_PARTIAL, 60, NO_MUT_RATE,3]
# exp_condition[18] = [SHARE_RATE_PARTIAL, 70, NO_MUT_RATE,3]
# exp_condition[19] = [SHARE_RATE_PARTIAL, 80, NO_MUT_RATE,3]
# exp_condition[20] = [SHARE_RATE_PARTIAL, 90, NO_MUT_RATE,3]
# exp_condition[21] = [SHARE_RATE_PARTIAL, SHARE_RADIUS_GLOBAL, NO_MUT_RATE,3]

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
            np.save("avg_"+str(k)+"_"+str(id)+"_150.npy", avg)

id = int(sys.argv[1])
sim = SLSim()
sim.exp(0,14,1,id,False)
