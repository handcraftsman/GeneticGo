package genetic

import (
	"fmt"
	"math/rand"
	"time"
)

type Solver struct {
	RandSeed                          int64
	MaxSecondsToRunWithoutImprovement float64
	LowerFitnessesAreBetter           bool
	PrintStrategyUsage                bool
}

func (solver *Solver) GetBest(getFitness func(string) int, display func(string), geneSet string, numberOfChromosomes, numberOfGenesPerChromosome int) string {
	if solver.RandSeed > 0 {
		rand.Seed(solver.RandSeed)
	} else {
		rand.Seed(time.Now().UnixNano())
	}
	if solver.MaxSecondsToRunWithoutImprovement == 0 {
		solver.MaxSecondsToRunWithoutImprovement = 20
		fmt.Printf("\tSolver will run at most %v second(s) without improvement.\n", solver.MaxSecondsToRunWithoutImprovement)
	}

	childFitnessIsBetter := func(childFitness, parentFitness int) bool {
		return childFitness > parentFitness
	}
	if solver.LowerFitnessesAreBetter {
		childFitnessIsBetter = func(childFitness, parentFitness int) bool {
			return childFitness < parentFitness
		}
	}

	childFitnessIsSameOrBetter := func(childFitness, parentFitness int) bool {
		return childFitness >= parentFitness
	}
	if solver.LowerFitnessesAreBetter {
		childFitnessIsSameOrBetter = func(childFitness, parentFitness int) bool {
			return childFitness <= parentFitness
		}
	}

	nextGene := make(chan string)
	go generateGene(nextGene, geneSet)

	nextChromosome := make(chan string)
	go generateChromosome(nextChromosome, nextGene, geneSet, numberOfGenesPerChromosome)

	start := time.Now()

	strategies := [...]strategyInfo{
		{"crossover ", crossover, 1, false},
		{"mutate    ", mutate, 1, true},
		{"reverse   ", reverse, 1, true},
		{"shift     ", shift, 1, true},
		{"swap      ", swap, 1, true},
	}

	strategySuccessSum := len(strategies)

	//pool := make([]sequenceInfo, 3*numberOfChromosomes*numberOfGenesPerChromosome)
	maxPoolSize := 3 * numberOfChromosomes * numberOfGenesPerChromosome
	pool := make([]sequenceInfo, 10)
	bestParent, bestFitness, distinctPool := populatePool(pool, nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome, childFitnessIsBetter, getFitness)
	defer func() { pool = nil }()
	defer func() { distinctPool = nil }()
	nextDistinctPool := make(map[string]bool, len(pool))
	nextDistinctPoolFitnesses := make(map[int]bool, len(pool))
	nextDistinctPoolFitnesses[pool[0].fitness] = true
	nextDistinctPoolFitnesses[pool[len(pool)/2].fitness] = true
	nextDistinctPoolFitnesses[pool[len(pool)-1].fitness] = true

	defer func() { nextDistinctPool = nil }()
	poolIndex := 0

	usedPoolParentSuccessCount := 80
	usedBestParentSuccessCount := 20
	usedParentSuccessSum := 100

	expandedPool := false

	for time.Since(start).Seconds() < solver.MaxSecondsToRunWithoutImprovement {
		strategyIndex := rand.Intn(strategySuccessSum)
		var strategy = strategies[0]
		for i, potentialStrategy := range strategies {
			if strategyIndex < potentialStrategy.count {
				strategyIndex = i
				strategy = strategies[strategyIndex]
				break
			}
			strategyIndex -= potentialStrategy.count
		}

		useBestParent := false
		whichParent := rand.Intn(usedParentSuccessSum)
		if whichParent > usedPoolParentSuccessCount {
			useBestParent = true
		}

		poolIndex = rand.Intn(len(pool))
		parent := pool[poolIndex]
		child := strategy.function(parent.genes, bestParent, geneSet, numberOfGenesPerChromosome, nextGene, useBestParent)
		if distinctPool[child] {
			continue
		}
		distinctPool[child] = true

		childFitness := getFitness(child)
		if childFitnessIsSameOrBetter(childFitness, pool[len(pool)-1].fitness) {
			if childFitnessIsBetter(childFitness, pool[len(pool)-1].fitness) {
				pool[len(pool)-1] = sequenceInfo{child, childFitness}
				insertionSort(pool, childFitnessIsBetter, len(pool)-1)
				nextDistinctPool[child] = true
				nextDistinctPoolFitnesses[childFitness] = true
			} else if len(pool) < maxPoolSize && len(nextDistinctPoolFitnesses) < 4 {
				pool = append(pool, sequenceInfo{child, childFitness})
				insertionSort(pool, childFitnessIsBetter, len(pool)-1)
				nextDistinctPoolFitnesses[childFitness] = true
				if solver.PrintStrategyUsage {
					print(".")
				}
				expandedPool = true
			}
			if childFitnessIsBetter(childFitness, bestFitness) {
				if solver.PrintStrategyUsage {
					if expandedPool {
						println()
						expandedPool = false
					}
					print(strategy.name)
				}
				display(child)
				start = time.Now()
				bestFitness = childFitness
				bestParent = child
				strategies[strategyIndex].count++
				strategySuccessSum++
			}
		}

		if len(nextDistinctPool) == len(pool) && len(nextDistinctPoolFitnesses) > 3 {
			distinctPool = nextDistinctPool
			nextDistinctPool = make(map[string]bool, len(pool))
			nextDistinctPool[bestParent] = true
			nextDistinctPoolFitnesses = make(map[int]bool, len(pool))
			nextDistinctPoolFitnesses[pool[0].fitness] = true
			nextDistinctPoolFitnesses[pool[len(pool)/2].fitness] = true
			nextDistinctPoolFitnesses[pool[len(pool)-1].fitness] = true
		}
	}

	if solver.PrintStrategyUsage {
		printStrategyUsage(strategies, strategySuccessSum, usedBestParentSuccessCount, usedPoolParentSuccessCount, usedParentSuccessSum)
	}

	return bestParent
}

func populatePool(pool []sequenceInfo, nextChromosome chan string, geneSet string, numberOfChromosomes, numberOfGenesPerChromosome int, compareFitnesses func(int, int) bool, getFitness func(string) int) (string, int, map[string]bool) {
	distinctPool := make(map[string]bool, len(pool))
	item := generateParent(nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome)
	bestFitness := getFitness(item)
	pool[0] = sequenceInfo{item, bestFitness}
	bestIndex := 0

	i := 1
	for i < len(pool) {
		item = generateParent(nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome)

		if distinctPool[item] {
			continue
		}
		distinctPool[item] = true

		itemFitness := getFitness(item)

		pool[i] = sequenceInfo{item, bestFitness}
		insertionSort(pool, compareFitnesses, i)
		if compareFitnesses(itemFitness, bestFitness) {
			bestFitness = itemFitness
			bestIndex = i
		}
		i++
	}

	return pool[bestIndex].genes, bestFitness, distinctPool
}

func insertionSort(items []sequenceInfo, compareFitnesses func(int, int) bool, index int) {
	if index < 1 || index > len(items) {
		return
	}
	for i := index; i > 0; i-- {
		if compareFitnesses(items[i].fitness, items[i-1].fitness) {
			items[i], items[i-1] = items[i-1], items[i]
			continue
		}
		break
	}
}

func printStrategyUsage(strategies [5]strategyInfo, strategySuccessSum, usedBestParentSuccessCount, usedPoolParentSuccessCount, usedParentSuccessSum int) {
	println("\nstrategy usage:")
	for _, strategy := range strategies {
		println(fmt.Sprint(
			strategy.name, "\t",
			strategy.count, "\t",
			100.0*strategy.count/strategySuccessSum, "%"))
	}
	println("\nparent selection in strategies:")
	println(fmt.Sprint(
		"used best parent\t",
		usedBestParentSuccessCount, "\t",
		100.0*usedBestParentSuccessCount/usedParentSuccessSum, "%"))
	println(fmt.Sprint(
		"used pool parent\t",
		usedPoolParentSuccessCount, "\t",
		100.0*usedPoolParentSuccessCount/usedParentSuccessSum, "%"))
	println()
}

func crossover(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string, useBestParent bool) string {
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

func reverse(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string, useBestParent bool) string {
	parent := parentA
	if useBestParent {
		parent = parentB
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

func generateParent(nextChromosome chan string, geneSet string, numberOfChromosomes, numberOfGenesPerChromosome int) string {
	s := ""
	for i := 0; i < numberOfChromosomes; i++ {
		s += <-nextChromosome
	}
	return s
}

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

func sort(a, b int) (int, int) {
	if a < b {
		return a, b
	}
	return b, a
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func reverseArray(a []string) {
	for i, j := 0, len(a)-1; i < j; i, j = i+1, j-1 {
		a[i], a[j] = a[j], a[i]
	}
}

type strategyInfo struct {
	name                         string
	function                     func(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string, useBestParent bool) string
	count                        int
	incrementParentSuccessCounts bool
}

type sequenceInfo struct {
	genes   string
	fitness int
}
