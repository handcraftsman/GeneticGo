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
	pool := make([]sequenceInfo, 10)
	defer func() { pool = nil }()

	distinctPool := populatePool(pool, nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome, childFitnessIsBetter, getFitness)
	defer func() { distinctPool = nil }()

	bestParent := pool[0]

	nextDistinctPool := make(map[string]bool, len(pool))
	defer func() { nextDistinctPool = nil }()

	nextDistinctPoolFitnesses := populateDistinctPoolFitnessesMap(pool)
	defer func() { nextDistinctPoolFitnesses = nil }()

	expandedPool := false

	for time.Since(start).Seconds() < solver.MaxSecondsToRunWithoutImprovement {
		strategyIndex := getNextStrategyIndex(strategies, strategySuccessSum)
		var strategy = strategies[strategyIndex]

		parent := pool[rand.Intn(len(pool))]
		useBestParent := rand.Intn(100) > 80

		childGenes := strategy.function(parent.genes, bestParent.genes, geneSet, numberOfGenesPerChromosome, nextGene, useBestParent)

		if distinctPool[childGenes] {
			continue
		}
		distinctPool[childGenes] = true

		child := sequenceInfo{childGenes, getFitness(childGenes)}

		if childFitnessIsSameOrBetter(child, pool[len(pool)-1]) {
			if childFitnessIsBetter(child, pool[len(pool)-1]) {
				pool[len(pool)-1] = child
				insertionSort(pool, childFitnessIsBetter, len(pool)-1)
				nextDistinctPool[child.genes] = true
				nextDistinctPoolFitnesses[child.fitness] = true
			} else if len(pool) < maxPoolSize && len(nextDistinctPoolFitnesses) < 4 {
				pool = append(pool, child)
				insertionSort(pool, childFitnessIsBetter, len(pool)-1)
				nextDistinctPoolFitnesses[child.fitness] = true
				if solver.PrintStrategyUsage {
					print(".")
				}
				expandedPool = true
			}
			if childFitnessIsBetter(child, bestParent) {
				if solver.PrintStrategyUsage {
					if expandedPool {
						println()
						expandedPool = false
					}
					print(strategy.name)
				}
				display(child.genes)
				start = time.Now()
				bestParent = child
				strategies[strategyIndex].count++
				strategySuccessSum++
			}
		}

		if len(nextDistinctPool) == len(pool) && len(nextDistinctPoolFitnesses) > 3 {
			distinctPool = nextDistinctPool
			nextDistinctPool = make(map[string]bool, len(pool))
			nextDistinctPool[bestParent.genes] = true

			nextDistinctPoolFitnesses = populateDistinctPoolFitnessesMap(pool)
		}
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
	nextDistinctPoolFitnesses := make(map[int]bool, len(pool))
	nextDistinctPoolFitnesses[pool[0].fitness] = true
	nextDistinctPoolFitnesses[pool[len(pool)/2].fitness] = true
	nextDistinctPoolFitnesses[pool[len(pool)-1].fitness] = true
	return nextDistinctPoolFitnesses
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
