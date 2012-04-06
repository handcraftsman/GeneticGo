package genetic

import (
	"math/rand"
)

func generateChromosome(nextChromosome, nextGene chan string, geneSet string, numberOfGenesPerChromosome int) {
	for {
		c := ""
		for i := 0; i < numberOfGenesPerChromosome; i++ {
			c += <-nextGene
		}
		nextChromosome <- c
	}
}

func generateGene(nextGene chan string, geneSet string) {
	for {
		index := rand.Intn(len(geneSet))
		nextGene <- geneSet[index : index+1]
	}
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

	i := 1
	for i < len(pool) {
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
