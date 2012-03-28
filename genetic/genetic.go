package genetic

import (
	"math/rand";
	"time";
)

type Solver struct {
}

func (s *Solver) GetBest(calculate func(string) int, disp func(string), genes string, length int) string {
	rand.Seed(time.Now().UnixNano())
	var bestGenes = generateParent(genes, length)
	value := calculate(bestGenes)
	var bestValue = value
	
	for bestValue < length {
		current := mutateParent(bestGenes, genes)
		value := calculate(current)
		if value > bestValue {
			disp(current)
			bestValue = value
			bestGenes = current
		}
	}

	return bestGenes
}

func mutateParent(parent, genes string) string {
	geneIndex := rand.Intn(len(genes))
	parentIndex := rand.Intn(len(parent))
	current := ""
	if parentIndex > 0 {
		current += parent[:parentIndex]
	}
	current += genes[geneIndex:1+geneIndex]
	if parentIndex+1 < len(parent) {
		current += parent[parentIndex+1:]
	}
	return current
}

func generateParent(genes string, length int) string {
	s := ""
	for i := 0; i < length; i++ {
		index := rand.Intn(len(genes))
		s += genes[index:1+index]
	}
	return s
}
