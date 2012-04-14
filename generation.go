package genetic

import (
	"bytes"
)

func generateChromosome(nextChromosome, nextGene chan string, geneSet string, numberOfGenesPerChromosome int, quit chan bool) {
	defer func() { close(nextChromosome) }()

	for {
		c := bytes.NewBuffer(make([]byte, 0, numberOfGenesPerChromosome))
		for i := 0; i < numberOfGenesPerChromosome; i++ {
			select {
			case <-quit:
				quit <- true
				return
			default:
				c.WriteString(<-nextGene)
			}
		}
		select {
		case <-quit:
			quit <- true
			return
		default:
			nextChromosome <- c.String()
		}
	}
}

func generateGene(nextGene chan string, geneSet string, quit chan bool, randSeed int64) {
	localRand := createRandomNumberGenerator(randSeed)
	defer func() { close(nextGene) }()
	for {
		index := localRand.Intn(len(geneSet))
		select {
		case <-quit:
			quit <- true
			return
		default:
			nextGene <- geneSet[index : index+1]
		}
	}
}

func generateParent(nextChromosome chan string, geneSet string, numberOfChromosomes, numberOfGenesPerChromosome int) string {
	s := bytes.NewBuffer(make([]byte, 0, numberOfChromosomes*numberOfGenesPerChromosome))
	for i := 0; i < numberOfChromosomes; i++ {
		s.WriteString(<-nextChromosome)
	}
	return s.String()
}

func populatePool(pool []sequenceInfo, nextChromosome chan string, geneSet string, numberOfChromosomes, numberOfGenesPerChromosome int, compareFitnesses func(sequenceInfo, sequenceInfo) bool, getFitness func(string) int) map[string]bool {
	distinctPool := make(map[string]bool, len(pool))
	itemGenes := generateParent(nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome)
	initialStrategy := strategyInfo{name: "initial   "}
	pool[0] = sequenceInfo{itemGenes, getFitness(itemGenes), initialStrategy}

	for i := 1; i < len(pool); {
		itemGenes = generateParent(nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome)

		if distinctPool[itemGenes] {
			continue
		}
		distinctPool[itemGenes] = true

		pool[i] = sequenceInfo{itemGenes, getFitness(itemGenes), initialStrategy}

		insertionSort(pool, compareFitnesses, i)
		i++
	}

	return distinctPool
}
