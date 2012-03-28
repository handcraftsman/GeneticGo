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
	
	for time.Since(start).Seconds() < solver.MaxSecondsToRunWithoutImprovement {
		current := mutateParent(bestParent, geneSet, nextGene)
		fitness := getFitness(current)
		if fitness > bestFitness {
			display(current)
			bestFitness = fitness
			bestParent = current
			start = time.Now()
		}
	}

	return bestParent
}

func mutateParent(parent, geneSet string, nextGene chan string) string {
	parentIndex := rand.Intn(len(parent))
	child := ""
	if parentIndex > 0 {
		child += parent[:parentIndex]
	}
	child += <- nextGene
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
