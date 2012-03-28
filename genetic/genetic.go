package genetic

import (
	"math/rand";
	"time";
)

type Solver struct {
	randSeed int64
}

func (s *Solver) GetBest(getFitness func(string) int, display func(string), geneSet string, numberOfChromosomes, numberOfGenesPerChromosome int) string {
	if s.randSeed > 0 {
		rand.Seed(s.randSeed)
	} else {
		rand.Seed(time.Now().UnixNano())
	}

	nextGene := make(chan string)
	go generateGene(nextGene, geneSet)
	
	nextChromosome := make(chan string)
	go generateChromosome(nextChromosome, nextGene, geneSet, numberOfGenesPerChromosome)
	
	var bestGenes = generateParent(nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome)
	value := getFitness(bestGenes)
	var bestValue = value
	
	for bestValue < numberOfChromosomes {
		current := mutateParent(bestGenes, geneSet, nextGene)
		value := getFitness(current)
		if value > bestValue {
			display(current)
			bestValue = value
			bestGenes = current
		}
	}

	return bestGenes
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
