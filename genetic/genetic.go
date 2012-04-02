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
	}
	fmt.Printf("\tSolver will run at most %v second(s) without improvement.\n", solver.MaxSecondsToRunWithoutImprovement)

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

	var bestParent = generateParent(nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome)
	fitness := getFitness(bestParent)
	bestFitness := fitness
	parent := bestParent
	parentFitness := fitness

	strategies := [...]strategyInfo{
		{"crossover ", crossover, 1},
		{"mutate    ", mutate, 1},
		{"reverse   ", reverse, 1},
		{"shift     ", shift, 1},
		{"swap      ", swap, 1},
	}

	strategySuccessSum := len(strategies)

	for time.Since(start).Seconds() < solver.MaxSecondsToRunWithoutImprovement {
		strategyIndex := rand.Intn(strategySuccessSum)
		var strategy = strategies[0]
		for i, potentialStrategy := range strategies {
			if strategyIndex < potentialStrategy.Count {
				strategyIndex = i
				strategy = strategies[strategyIndex]
				break
			}
			strategyIndex -= potentialStrategy.Count
		}

		child := strategy.function(parent, bestParent, geneSet, numberOfGenesPerChromosome, nextGene)
		fitness := getFitness(child)
		if childFitnessIsSameOrBetter(fitness, parentFitness) {
			if childFitnessIsBetter(fitness, bestFitness) {
				if solver.PrintStrategyUsage {
					fmt.Print(strategy.name)
				}
				display(child)
				start = time.Now()
				bestFitness = fitness
				bestParent = child
				strategies[strategyIndex].Count++
				strategySuccessSum++
			} else {
				parent = child
				parentFitness = fitness
			}
		}
	}

	if solver.PrintStrategyUsage {
		printStrategyUsage(strategies, strategySuccessSum)
	}

	return bestParent
}

func printStrategyUsage(strategies [5]strategyInfo, strategySuccessSum int) {
	println("\nstrategy usage:")
	for _, strategy := range strategies {
		println(fmt.Sprint(
			strategy.name, "\t",
			strategy.Count, "\t",
			100.0*strategy.Count/strategySuccessSum, "%"))
	}
}

func crossover(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string) string {
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

func reverse(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string) string {
	parent := parentB
	useParentA := rand.Intn(2) == 1
	if useParentA {
		parent = parentA
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

func shift(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string) string {
	parent := parentB
	useParentA := rand.Intn(2) == 1
	if useParentA {
		parent = parentA
	}

	if len(parent) < numberOfGenesPerChromosome+1 {
		return mutate(parent, parentB, geneSet, numberOfGenesPerChromosome, nextGene)
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

func swap(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string) string {
	parent := parentB
	useParentA := rand.Intn(2) == 1
	if useParentA {
		parent = parentA
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

func mutate(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string) string {
	parent := parentB
	useParentA := rand.Intn(2) == 1
	if useParentA {
		parent = parentA
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
	name     string
	function func(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string) string
	Count    int
}
