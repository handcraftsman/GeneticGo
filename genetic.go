package genetic

import (
	"fmt"
	"math/rand"
	"time"
)

func (solver *Solver) GetBest(getFitness func(string) int, display func(string), geneSet string, numberOfChromosomes, numberOfGenesPerChromosome int) string {
	seedRandomNumberGenerator(solver.RandSeed)
	solver.ensureMaxSecondsToRunIsValid()
	childFitnessIsBetter, childFitnessIsSameOrBetter := solver.createFitnessComparisonFunctions()

	nextGene := make(chan string)
	go generateGene(nextGene, geneSet)

	nextChromosome := make(chan string)
	go generateChromosome(nextChromosome, nextGene, geneSet, numberOfGenesPerChromosome)

	strategies := []strategyInfo{
		{"crossover ", crossover, 1, false},
		{"mutate    ", mutate, 1, true},
		{"reverse   ", reverse, 1, true},
		{"shift     ", shift, 1, true},
		{"swap      ", swap, 1, true},
	}
	strategySuccessSum := len(strategies)

	start := time.Now()

	maxPoolSize := 3 * numberOfChromosomes * numberOfGenesPerChromosome

	pool := make([]sequenceInfo, maxPoolSize, maxPoolSize)
	defer func() { pool = nil }()

	distinctPool := populatePool(pool, nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome, childFitnessIsBetter, getFitness)
	defer func() { distinctPool = nil }()

	bestParent := pool[0]

	children := make([]sequenceInfo, 1, maxPoolSize)
	defer func() { children = nil }()
	children[0] = bestParent

	distinctChildren := make(map[string]bool, len(pool))
	defer func() { distinctChildren = nil }()

	distinctChildrenFitnesses := populateDistinctPoolFitnessesMap(pool)
	defer func() { distinctChildrenFitnesses = nil }()

	needNewlineBeforeDisplay := false

	for time.Since(start).Seconds() < solver.MaxSecondsToRunWithoutImprovement {
		strategyIndex := getNextStrategyIndex(strategies, strategySuccessSum)
		var strategy = strategies[strategyIndex]

		parent := pool[rand.Intn(len(pool))]
		useBestParent := rand.Intn(100) == 0

		childGenes := strategy.function(parent.genes, bestParent.genes, geneSet, numberOfGenesPerChromosome, nextGene, useBestParent)

		if distinctPool[childGenes] {
			continue
		}
		distinctPool[childGenes] = true

		child := sequenceInfo{childGenes, getFitness(childGenes)}

		if !childFitnessIsSameOrBetter(child, pool[len(pool)-1]) {
			continue
		}

		if child.fitness == pool[len(pool)-1].fitness {
			if len(pool) < maxPoolSize {
				pool = append(pool, child)
			} else {
				pool[len(pool)-1] = child
				insertionSort(pool, childFitnessIsSameOrBetter, len(pool)-1)
			}

			continue
		}

		if len(children) < maxPoolSize &&
			(len(distinctChildrenFitnesses) < 4 ||
				child.fitness == children[len(children)-1].fitness) {

			children = append(children, child)

			if solver.PrintDiagnosticInfo {
				print(".")
				needNewlineBeforeDisplay = true
			}
		} else {
			children[len(children)-1] = child
		}

		insertionSort(children, childFitnessIsSameOrBetter, len(children)-1)

		distinctChildren[child.genes] = true
		distinctChildrenFitnesses[child.fitness] = true

		if childFitnessIsBetter(child, bestParent) {
			if needNewlineBeforeDisplay {
				println()
				needNewlineBeforeDisplay = false
			}
			if solver.PrintStrategyUsage {
				print(strategy.name)
			}
			display(child.genes)
			start = time.Now()
			bestParent = child
			strategies[strategyIndex].count++
			strategySuccessSum++
			if solver.PrintDiagnosticInfo {
				print("+")
				needNewlineBeforeDisplay = true
			}

			pool[len(pool)-1] = child
			insertionSort(pool, childFitnessIsSameOrBetter, len(pool)-1)
		}

		if len(distinctChildren) >= maxPoolSize &&
			len(distinctChildrenFitnesses) > 3 {
			if solver.PrintDiagnosticInfo {
				print(">")
				needNewlineBeforeDisplay = true
			}

			pool = children
			children = make([]sequenceInfo, 1, maxPoolSize)
			children[0] = bestParent

			distinctPool = distinctChildren
			distinctChildren = make(map[string]bool, len(children))
			distinctChildren[bestParent.genes] = true

			distinctChildrenFitnesses = make(map[int]bool, len(pool))
			distinctChildrenFitnesses[bestParent.fitness] = true
		}
	}

	if needNewlineBeforeDisplay {
		println()
	}
	if solver.PrintStrategyUsage {
		printStrategyUsage(strategies, strategySuccessSum)
	}

	return bestParent.genes
}

func (solver *Solver) createFitnessComparisonFunctions() (func(child, other sequenceInfo) bool, func(child, other sequenceInfo) bool) {
	childFitnessIsBetter := func(child, other sequenceInfo) bool {
		return child.fitness < other.fitness
	}

	childFitnessIsSameOrBetter := func(child, other sequenceInfo) bool {
		return child.fitness <= other.fitness
	}

	if !solver.LowerFitnessesAreBetter {
		childFitnessIsBetter = func(child, other sequenceInfo) bool {
			return child.fitness > other.fitness
		}

		childFitnessIsSameOrBetter = func(child, other sequenceInfo) bool {
			return child.fitness >= other.fitness
		}
	}

	return childFitnessIsBetter, childFitnessIsSameOrBetter
}

func (solver *Solver) ensureMaxSecondsToRunIsValid() {
	if solver.MaxSecondsToRunWithoutImprovement == 0 {
		solver.MaxSecondsToRunWithoutImprovement = 20
		fmt.Printf("\tSolver will run at most %v second(s) without improvement.\n", solver.MaxSecondsToRunWithoutImprovement)
	}
}

func getNextStrategyIndex(strategies []strategyInfo, strategySuccessSum int) int {
	strategyIndex := rand.Intn(strategySuccessSum)
	for i, potentialStrategy := range strategies {
		if strategyIndex < potentialStrategy.count {
			strategyIndex = i
			break
		}
		strategyIndex -= potentialStrategy.count
	}
	return strategyIndex
}

func populateDistinctPoolFitnessesMap(pool []sequenceInfo) map[int]bool {
	distinctChildrenFitnesses := make(map[int]bool, len(pool))
	distinctChildrenFitnesses[pool[0].fitness] = true
	distinctChildrenFitnesses[pool[len(pool)/2].fitness] = true
	distinctChildrenFitnesses[pool[len(pool)-1].fitness] = true
	return distinctChildrenFitnesses
}

func printStrategyUsage(strategies []strategyInfo, strategySuccessSum int) {
	println("\nstrategy usage:")
	for _, strategy := range strategies {
		println(fmt.Sprint(
			strategy.name, "\t",
			strategy.count, "\t",
			100.0*strategy.count/strategySuccessSum, "%"))
	}
	println()
}
