import numpy as np

import networkx as nx

import matplotlib.pyplot as plt



def ind2gen(index,n):

    '''Return a genotype of length n that encodes the index in binary'''

    # For example, ind2gen(255,8) = [1, 1, 1, 1, 1, 1, 1, 1]

    genotype = np.zeros(n)

    while n > 0:

        n = n - 1

        if index % 2 == 0:

            genotype[n] = 0

        else:

            genotype[n] = 1

        index = index // 2

    return genotype



def gen2ind(genotype, n):

    '''Return the index encoded in the genotype'''

    # For example, gen2ind([1,1,1,1,1,1,1,1]) = 255

    i = 0

    index = 0

    while i < n:

        index += genotype[i]*(2**(n-i-1))

        i += 1

    return int(index)



class NKLandscape:

    '''Create a tunably rugged landscape with N dimensions and K epistatic interactions'''



    def __init__(self, n, k):

        self.n = n

        self.k = k

        self.maxfit = 0.0

        self.minfit = 1000000

        # Create random matrix for interactions

        # For example, when N = 5 and K = 1, there are 5 rows each with 4 columns

        # to score each combination of an individual and its neighbor (i.e. 01, 10, 11, 00)

        self.interactions = np.random.rand(n*(2**(k+1))).reshape(n,2**(k+1))

        self.fit_table = np.zeros(2**n)

        self.best = 0

        for solution in range(2**n):      # Figure out fitness for each possible solution

            genotype = ind2gen(solution,n)

            fit = self.calcFit(genotype)

            self.fit_table[solution] = fit

            if fit > self.maxfit:

                self.maxfit = fit

                self.best = genotype

            if fit < self.minfit:

                self.minfit = fit

        self.fit_table = (self.fit_table - self.minfit)/(self.maxfit-self.minfit) # Normalize

        self.fit_table = self.fit_table**8    # Scale



    def calcFit(self, genotype):

        fit = 0.0

        for gene in range(self.n):  # Calculate contribution of each gene in the current solution

            subgen = []

            for nbr in range(self.k + 1):  # Identify neighbors

                nbr_ind = (gene + nbr) % self.n

                subgen.append(genotype[nbr_ind])

            ind = gen2ind(subgen, self.k + 1)  # Calculate epistatic interactions with each neighbor

            fit += self.interactions[gene][ind]

        return fit



    def fitness(self, genotype):

        '''Return the fitness of a solution.'''

        index = gen2ind(genotype,self.n)

        return self.fit_table[index]



class Population:

    ''''''



    def __init__(self, popsize, n, landscape, community=False):

        self.popsize = popsize  # population size

        self.ng = n              # number of genes

        self.share_rate = 1.0   # recombination rate

        self.share_radius = popsize - 1  # how big your "neighborhood" is

        self.mut_rate = 0.0 #(E:May14) Percentage of the solution that is mutated during a copy/share

        self.landscape = landscape          # NK landscape

        self.genotypes = np.random.randint(2,size=popsize*n).reshape(popsize,n)

        self.learned = np.zeros(popsize,dtype=int) # Keeps track of individuals who just learned

        self.shared = np.zeros(popsize,dtype=int)  # and who just shared

        #self.dist_list = np.zeros(int(((popsize*popsize)-popsize)/2))

        self.community = community



    def set_community(self, ncomm, in_group, out_group):

        pmatrix = np.eye(ncomm) * in_group + np.where(np.eye(ncomm) == 1, 0, 1) * out_group

        csizes = [int(self.popsize / ncomm)] * ncomm

        graph = nx.stochastic_block_model(sizes=csizes, p=pmatrix)

        self.adj_matrix = nx.to_numpy_array(graph)



    def set_pop(self, genotypes):

        '''Set the population genotypes to the given genotypes'''

        self.genotypes = genotypes.copy()



    def avgfit(self):

        '''Return the average fitness of the population '''

        avg = 0

        for i in range(self.popsize):

            avg += self.landscape.fitness(self.genotypes[i])

        return avg / self.popsize



    def stats(self):

        '''Return the average fitness of the population and the fitness of the best individual'''

        #Calculate Avg and Best Fitness

        avg = 0

        best = 0

        for i in range(self.popsize):

            fit = self.landscape.fitness(self.genotypes[i])

            self.landscape.visited(self.genotypes[i])

            avg += fit

            if fit > best:

                best = fit

        # Calculate Avg Hamming Distance

        # Also calculate Spread (i.e., number of unique solutions)

        k=0

        unique = np.ones(self.popsize)

        for i in range(self.popsize):

            for j in range(i+1,self.popsize):

                self.dist_list[k] = np.mean(np.abs(self.genotypes[i] - self.genotypes[j]))

                if (self.dist_list[k]==0):

                    unique[i]=0.0

                    unique[j]=0.0

                k+=1

        return avg/self.popsize, np.mean(self.dist_list), np.mean(unique)



    def learn(self): # When pref=0, learning happens regardless of if sharing didn't happen

        '''Update the genotypes of the population by having them learn'''

        self.learned = np.zeros(self.popsize,dtype=int)

        for i in range(self.popsize):

            if self.shared[i]==0:

                original_fitness = self.landscape.fitness(self.genotypes[i])

                new_genotype = self.genotypes[i].copy()

                j = np.random.randint(self.ng)  # choose a random gene and "flip" it

                if new_genotype[j] == 0:

                    new_genotype[j] = 1

                else:

                    new_genotype[j] = 0

                new_fitness = self.landscape.fitness(new_genotype)

                if new_fitness > original_fitness:   # If the change is better, we keep it.

                    self.learned[i] = 1

                    self.genotypes[i] = new_genotype.copy()



    def share(self):

        '''Update the genotypes of the population by having them share'''

        self.shared = np.zeros(self.popsize,dtype=int)

        new_genotypes = self.genotypes.copy()

        for i in range(self.popsize):

            if self.community:

                neighbors = [j for j in range(self.popsize) if self.adj_matrix[i, j] == 1]

                j = np.random.choice(neighbors)

            else:

                # Pick a neighbor in your radius

                j = np.random.randint(i-self.share_radius,i+self.share_radius+1)

                while (j == i):

                    j = np.random.randint(i-self.share_radius,i+self.share_radius+1)

                j = j % self.popsize;

            # Compare fitness with neighbor

            if self.landscape.fitness(self.genotypes[j]) > self.landscape.fitness(self.genotypes[i]):

                self.shared[i] = 1

                # If neighbor is better, get some of their answers

                for g in range(self.ng):

                    if np.random.rand() <= self.share_rate:

                        new_genotypes[i][g] = self.genotypes[j][g]

                        if np.random.rand() <= self.mut_rate:

                            new_genotypes[i][g] = np.random.randint(2)

        self.genotypes = new_genotypes.copy()

