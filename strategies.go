package genetic

import (
	"math/rand"
)

func crossover(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string, useBestParent bool) string {
	if len(parentA) == 1 || len(parentB) == 1 {
		return mutate(parentA, parentB, geneSet, numberOfGenesPerChromosome, nextGene, useBestParent)
	}

	sourceStart := rand.Intn((len(parentB)-1)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
	destinationStart := rand.Intn((len(parentA)-1)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
	maxLength := min(len(parentA)-destinationStart, len(parentB)-sourceStart)
	length := (1 + rand.Intn(maxLength/numberOfGenesPerChromosome-1)) * numberOfGenesPerChromosome

	child := ""

	if destinationStart > 0 {
		child += parentA[0:destinationStart]
	}

	child += parentB[sourceStart : sourceStart+length]

	if len(child) < len(parentA) {
		child += parentA[len(child):len(parentA)]
	}

	return child
}

func mutate(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string, useBestParent bool) string {
	parent := parentA
	if useBestParent {
		parent = parentB
	}

	parentIndex := rand.Intn(len(parent))
	child := ""
	if parentIndex > 0 {
		child += parent[:parentIndex]
	}

	newGene := <-nextGene
	currentGene := parent[parentIndex : parentIndex+1]
	for ; newGene == currentGene; newGene = <-nextGene {
	}
	child += newGene

	if parentIndex+1 < len(parent) {
		child += parent[parentIndex+1:]
	}
	return child
}

func reverse(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string, useBestParent bool) string {
	parent := parentA
	if useBestParent {
		parent = parentB
	}

	if len(parent) == 1 {
		return mutate(parentA, parentB, geneSet, numberOfGenesPerChromosome, nextGene, useBestParent)
	}

	reversePointA := rand.Intn(len(parent)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
	reversePointB := rand.Intn(len(parent)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
	for ; reversePointA == reversePointB; reversePointB = rand.Intn(len(parent)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome {
	}

	min, max := sort(reversePointA, reversePointB)

	fragments := make([]string, max-min)
	for i := min; i < max; i += numberOfGenesPerChromosome {
		fragments[i-min] = parent[i : i+numberOfGenesPerChromosome]
	}

	child := ""
	if min > 0 {
		child = parent[0:min]
	}

	reverseArray(fragments)
	for _, fragment := range fragments {
		child += fragment
	}

	if len(child) < len(parent) {
		child += parent[len(child):len(parent)]
	}

	return child
}

func shift(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string, useBestParent bool) string {
	parent := parentA
	if useBestParent {
		parent = parentB
	}

	if len(parent) < numberOfGenesPerChromosome+1 {
		return mutate(parent, parentB, geneSet, numberOfGenesPerChromosome, nextGene, useBestParent)
	}
	shiftRight := rand.Intn(2) == 1

	segmentStart := rand.Intn((len(parent)-numberOfGenesPerChromosome)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
	segmentCount := 1 + rand.Intn((len(parent)-segmentStart+numberOfGenesPerChromosome)/numberOfGenesPerChromosome-1)

	// +2 because first and last will be empty to leave room
	fragments := make([]string, segmentCount+2)
	for i := 0; i < segmentCount; i++ {
		start := segmentStart + i*numberOfGenesPerChromosome
		fragments[i+1] = parent[start : start+numberOfGenesPerChromosome]
	}
	start := 1
	end := segmentCount

	if shiftRight {
		start = 0
		fragments[0] = fragments[end]
		end--
	} else {
		end++
		fragments[end] = fragments[0]
		start++
	}

	child := ""
	if segmentStart > 0 {
		child += parent[0:segmentStart]
	}

	for i := start; i <= end; i++ {
		child += fragments[i]
	}

	if len(child) < len(parent) {
		child += parent[len(child):len(parent)]
	}

	return child
}

func swap(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string, useBestParent bool) string {
	if len(parentA) == 1 || len(parentB) == 1 {
		return mutate(parentA, parentB, geneSet, numberOfGenesPerChromosome, nextGene, useBestParent)
	}

	parent := parentA
	if useBestParent {
		parent = parentB
	}

	parentIndexA := rand.Intn(len(parent))
	parentIndexB := rand.Intn(len(parent))
	for ; parentIndexA == parentIndexB; parentIndexB = rand.Intn(len(parent)) {
	}

	parentIndexA, parentIndexB = sort(parentIndexA, parentIndexB)

	child := ""
	if parentIndexA > 0 {
		child += parent[:parentIndexA]
	}

	child += parent[parentIndexB : parentIndexB+1]

	if parentIndexB-parentIndexA > 1 {
		child += parent[parentIndexA+1 : parentIndexB]
	}

	child += parent[parentIndexA : parentIndexA+1]

	if parentIndexB+1 < len(parent) {
		child += parent[parentIndexB+1:]
	}

	return child
}
