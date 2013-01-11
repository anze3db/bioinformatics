from random import randint, sample, random, choice
from itertools import izip

class Genetic:

    def __init__(self):
        
        self.weights = [randint(12,50) for _ in range(100)]
        #self.weights = [1,2,3,1,2,3,1,2,3,1,2,3]
        self.GENOME_SIZE = len(self.weights)
        self.RATE = 0.01
        self.EPOCHS = 100
        
        self.target_weight = 100
        self.n = 100
        self.pop = self.initial_population(self.n)
        
    def initial_population(self, n):
        return [[randint(0,1) for _ in range(self.GENOME_SIZE)] for _ in range(self.n)]
    
    def weight(self, genome):
        return sum(w for g, w in izip(genome, self.weights) if g)
    
    def fitness(self, genome):
        if self.weight(genome) > self.target_weight:
            return (self.weight(genome) - self.target_weight)*10
        return (self.target_weight - self.weight(genome))
    
    def crossover(self, g1, g2):
        site = randint(0, self.GENOME_SIZE)
        return g1[:site] + g2[site:], g2[:site] + g1[site:]
    
    def mutation(self, genome):
        return [g if random()  > self.RATE else randint(0,1) for g in genome]

    
    def evolve(self):

        for t in range(self.EPOCHS):
            pop_fitness = sorted([(self.fitness(g), g) for g in self.pop])
            self.alpha = pop_fitness[0][1]
            next_pop = [g for _, g in pop_fitness[:int(self.n*.1)]]
            kids = []
            for _ in range(int(self.n*.6*.5)):
               kids.extend(self.crossover(*sample(next_pop, 2)))
            
            next_pop.extend(kids)   
            next_pop = [self.mutation(k) for k in next_pop]
            self.pop = next_pop
            
            self.pop = next_pop     
        
    def print_final_alpha(self):
        print "Final alpha:"
        print sum(w for w, g in izip(self.weights, self.alpha) if g) 
        print " ".join("%d" % w for w, g in izip(self.weights, self.alpha) if g) 
        

if __name__ == '__main__':
    
    g = Genetic()
    g.evolve()
    g.print_final_alpha()
