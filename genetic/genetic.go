package genetic

import (
	"fmt";
	"math/rand";
	"time";
)

type Solver struct {
	RandSeed int64
	MaxSecondsToRunWithoutImprovement float64
}

func (solver *Solver) GetBest(getFitness func(string) int, display func(string), geneSet string, numberOfChromosomes, numberOfGenesPerChromosome int) string {
	if solver.RandSeed > 0 {
		rand.Seed(solver.RandSeed)
	} else {
		rand.Seed(time.Now().UnixNano())
	}
	if solver.MaxSecondsToRunWithoutImprovement == 0 {
		solver.MaxSecondsToRunWithoutImprovement = 20
	}
	fmt.Printf("\tSolver will run at most %v second(s) without improvement.\n", solver.MaxSecondsToRunWithoutImprovement)

	nextGene := make(chan string)
	go generateGene(nextGene, geneSet)
	
	nextChromosome := make(chan string)
	go generateChromosome(nextChromosome, nextGene, geneSet, numberOfGenesPerChromosome)
	
	start := time.Now()
	
	var bestParent = generateParent(nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome)
	fitness := getFitness(bestParent)
	var bestFitness = fitness
	
	strategies := []func (parent, geneSet string, nextGene chan string) string { mutate, swap }
	
	for time.Since(start).Seconds() < solver.MaxSecondsToRunWithoutImprovement {
		strategyIndex := rand.Intn(len(strategies))
		strategy := strategies[strategyIndex]
		current := strategy(bestParent, geneSet, nextGene)
		fitness := getFitness(current)
		if fitness >= bestFitness {
			if fitness > bestFitness {
				display(current)
				start = time.Now()
				bestFitness = fitness
			}
			bestParent = current
		}
	}

	return bestParent
}

func swap(parent, geneSet string, nextGene chan string) string {
	parentIndexA := rand.Intn(len(parent))
	parentIndexB := rand.Intn(len(parent))
	for ; parentIndexA == parentIndexB; parentIndexB = rand.Intn(len(parent)) {}
	
	parentIndexA, parentIndexB = sort(parentIndexA, parentIndexB)
	
	child := ""
	if parentIndexA > 0 {
		child += parent[:parentIndexA]
	}
	
	child += parent[parentIndexB:parentIndexB+1]
	
	if parentIndexB - parentIndexA > 1 {
		child += parent[parentIndexA+1:parentIndexB]
	}

	child += parent[parentIndexA:parentIndexA+1]

	if parentIndexB+1 < len(parent) {
		child += parent[parentIndexB+1:]
	}

	return child
}

func mutate(parent, geneSet string, nextGene chan string) string {
	parentIndex := rand.Intn(len(parent))
	child := ""
	if parentIndex > 0 {
		child += parent[:parentIndex]
	}
	
	newGene := <- nextGene
	currentGene := parent[parentIndex:parentIndex+1]
	for ; newGene == currentGene ; newGene = <- nextGene {}
	child += newGene
	
	if parentIndex+1 < len(parent) {
		child += parent[parentIndex+1:]
	}
	return child
}

func generateParent(nextChromosome chan string, geneSet string, numberOfChromosomes, numberOfGenesPerChromosome int) string {
	
	s := ""
	for i := 0; i < numberOfChromosomes; i++ {
		s += <- nextChromosome
	}
	return s
}

func generateChromosome(nextChromosome, nextGene chan string, geneSet string, numberOfGenesPerChromosome int) {
	for {
		c := ""
		for i := 0; i < numberOfGenesPerChromosome; i++ {
			c += <- nextGene
		}
		nextChromosome <- c
	}
}

func generateGene(nextGene chan string, geneSet string) {
	for {
		index := rand.Intn(len(geneSet))
		nextGene <- geneSet[index:index+1]	
	}
}

func sort(a, b int) (int, int) {
	if a < b {
		return a,b
	}
	return b,a
}
