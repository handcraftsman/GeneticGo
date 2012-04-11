package genetic

import (
	"math/rand"
)

func generateChromosome(nextChromosome, nextGene chan string, geneSet string, numberOfGenesPerChromosome int, quit *bool) {
	for {
		c := ""
		for i := 0; i < numberOfGenesPerChromosome && !*quit; i++ {
			c += <-nextGene
		}
		if *quit {
			break
		}
		nextChromosome <- c
	}
	close(nextChromosome)
}

func generateGene(nextGene chan string, geneSet string, quit *bool) {
	for {
		index := rand.Intn(len(geneSet))
		if *quit {
			break
		}
		nextGene <- geneSet[index : index+1]
	}
	close(nextGene)
}

func generateParent(nextChromosome chan string, geneSet string, numberOfChromosomes, numberOfGenesPerChromosome int) string {
	s := ""
	for i := 0; i < numberOfChromosomes; i++ {
		s += <-nextChromosome
	}
	return s
}

func populatePool(pool []sequenceInfo, nextChromosome chan string, geneSet string, numberOfChromosomes, numberOfGenesPerChromosome int, compareFitnesses func(sequenceInfo, sequenceInfo) bool, getFitness func(string) int) map[string]bool {
	distinctPool := make(map[string]bool, len(pool))
	itemGenes := generateParent(nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome)
	pool[0] = sequenceInfo{itemGenes, getFitness(itemGenes)}

	
	for i := 1; i < len(pool); {
		itemGenes = generateParent(nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome)

		if distinctPool[itemGenes] {
			continue
		}
		distinctPool[itemGenes] = true

		pool[i] = sequenceInfo{itemGenes, getFitness(itemGenes)}

		insertionSort(pool, compareFitnesses, i)
		i++
	}

	return distinctPool
}
