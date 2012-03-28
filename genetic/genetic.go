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
	
	cs := make(chan string)
	go generateChromosome(cs, geneSet, numberOfGenesPerChromosome)
	
	var bestGenes = generateParent(cs, geneSet, numberOfChromosomes, numberOfGenesPerChromosome)
	value := getFitness(bestGenes)
	var bestValue = value
	
	for bestValue < numberOfChromosomes {
		current := mutateParent(bestGenes, geneSet)
		value := getFitness(current)
		if value > bestValue {
			display(current)
			bestValue = value
			bestGenes = current
		}
	}

	return bestGenes
}

func mutateParent(parent, geneSet string) string {
	geneSetIndex := rand.Intn(len(geneSet))
	parentIndex := rand.Intn(len(parent))
	current := ""
	if parentIndex > 0 {
		current += parent[:parentIndex]
	}
	current += geneSet[geneSetIndex:1+geneSetIndex]
	if parentIndex+1 < len(parent) {
		current += parent[parentIndex+1:]
	}
	return current
}

func generateParent(cs chan string, geneSet string, numberOfChromosomes, numberOfGenesPerChromosome int) string {
	
	s := ""
	for i := 0; i < numberOfChromosomes; i++ {
		chromosome := <- cs
		s += chromosome
	}
	return s
}

func generateChromosome(cs chan string, geneSet string, numberOfGenesPerChromosome int) {
	for {
		c := ""
		for i := 0; i < numberOfGenesPerChromosome; i++ {
			index := rand.Intn(len(geneSet))
			c += geneSet[index:1+index]
		}
		cs <- c
	}
}

