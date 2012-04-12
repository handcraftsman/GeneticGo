package genetic

import (
	"fmt"
	"math"
	"time"
)

func (solver *Solver) GetBestUsingHillClimbing(getFitness func(string) int,
	display func(string),
	geneSet string,
	maxRoundsWithoutImprovement, numberOfGenesPerChromosome int,
	bestPossibleFitness int) string {

	solver.initialize(geneSet, numberOfGenesPerChromosome)

	roundsSinceLastImprovement := 0
	numberOfGenesToUse := numberOfGenesPerChromosome
	maxLength := numberOfGenesPerChromosome * maxRoundsWithoutImprovement
	generationCount := 1

	bestEver := sequenceInfo{genes: generateParent(solver.nextChromosome, geneSet, generationCount, numberOfGenesPerChromosome)}
	bestEver.fitness = getFitness(bestEver.genes)
	initialParent := bestEver
	display(bestEver.genes)

	solver.initializePool(generationCount, numberOfGenesPerChromosome, geneSet, initialParent, getFitness)

	defer func() {
		solver.quit = true
		<-solver.nextChromosome
		<-solver.nextGene
		solver.pool = nil
		solver.distinctPool = nil
	}()

	filteredDisplay := func(item sequenceInfo) {
		if solver.childFitnessIsBetter(item, bestEver) {
			display(item.genes)
		}
	}

	for ; numberOfGenesToUse < maxLength &&
		roundsSinceLastImprovement < maxRoundsWithoutImprovement &&
		bestEver.fitness != bestPossibleFitness; generationCount++ {

		solver.printNewlineIfNecessary()
		println("-- starting round", generationCount)

		result := solver.getBestWithInitialParent(getFitness,
			filteredDisplay,
			geneSet,
			generationCount,
			numberOfGenesPerChromosome)

		if solver.childFitnessIsBetter(result, bestEver) {
			roundsSinceLastImprovement = 0
			bestEver = result
		} else {
			roundsSinceLastImprovement++
			if roundsSinceLastImprovement >= maxRoundsWithoutImprovement {
				break
			}
		}

		newPool := make([]sequenceInfo, 0, solver.maxPoolSize)
		distinctPool := make(map[string]bool, solver.maxPoolSize)

		improved := false

		for round := 0; round < 100 && !improved; round++ {
			for _, parent := range solver.pool {
				childGenes := parent.genes + <-solver.nextChromosome
				if distinctPool[childGenes] {
					continue
				}
				distinctPool[childGenes] = true

				fitness := getFitness(childGenes)
				child := sequenceInfo{childGenes, fitness}
				if len(newPool) < solver.maxPoolSize {
					newPool = append(newPool, child)
				} else {
					newPool[len(newPool)-1] = child
				}
				insertionSort(newPool, solver.childFitnessIsSameOrBetter, len(newPool)-1)

				if solver.childFitnessIsBetter(child, bestEver) {
					roundsSinceLastImprovement = 0
					solver.printNewlineIfNecessary()
					if solver.PrintStrategyUsage {
						print("climb    ")
					}
					filteredDisplay(child)
					bestEver = child
					improved = true
				}
			}
		}

		solver.pool = newPool
		solver.distinctPool = distinctPool
	}

	solver.printNewlineIfNecessary()
	solver.printStrategyUsage()

	return bestEver.genes
}

func (solver *Solver) GetBest(getFitness func(string) int,
	display func(string),
	geneSet string,
	numberOfChromosomes, numberOfGenesPerChromosome int) string {

	solver.initialize(geneSet, numberOfGenesPerChromosome)

	defer func() {
		solver.quit = true
		<-solver.nextChromosome
		<-solver.nextGene
		solver.pool = nil
		solver.distinctPool = nil
	}()

	initialParent := sequenceInfo{genes: generateParent(solver.nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome)}
	initialParent.fitness = getFitness(initialParent.genes)

	solver.initializePool(numberOfChromosomes, numberOfGenesPerChromosome, geneSet, initialParent, getFitness)

	best := solver.getBestWithInitialParent(getFitness,
		func(sequence sequenceInfo) { display(sequence.genes) },
		geneSet,
		numberOfChromosomes,
		numberOfGenesPerChromosome)

	solver.printNewlineIfNecessary()
	solver.printStrategyUsage()

	return best.genes
}

func (solver *Solver) getBestWithInitialParent(getFitness func(string) int,
	display func(sequenceInfo),
	geneSet string,
	numberOfChromosomes, numberOfGenesPerChromosome int) sequenceInfo {

	start := time.Now()

	bestParent := solver.pool[0]

	children := make([]sequenceInfo, 1, solver.maxPoolSize)
	children[0] = bestParent

	distinctChildren := make(map[string]bool, len(solver.pool))
	distinctChildrenFitnesses := populateDistinctPoolFitnessesMap(solver.pool)

	quit := false
	nextStrategyIndex := make(chan int, 10)
	go func() {
		for {
			strategyIndex := solver.getNextStrategyIndex()
			if quit {
				break
			}
			nextStrategyIndex <- strategyIndex
		}
		close(nextStrategyIndex)
	}()

	defer func() {
		quit = true
		<-nextStrategyIndex
	}()

	maxPossiblePermutations := math.Pow(float64(len(geneSet)), float64(numberOfChromosomes*numberOfGenesPerChromosome))

	for time.Since(start).Seconds() < solver.MaxSecondsToRunWithoutImprovement {
		strategyIndex := <-nextStrategyIndex
		var strategy = solver.strategies[strategyIndex]

		parent := solver.pool[solver.nextRand(len(solver.pool))]
		useBestParent := solver.nextRand(100) == 0

		childGenes := strategy.generate(strategy, solver.mutationStrategy, parent.genes, bestParent.genes, geneSet, numberOfGenesPerChromosome, solver.nextGene, useBestParent)

		if solver.distinctPool[childGenes] {
			if float64(len(solver.distinctPool)) == maxPossiblePermutations {
				solver.printNewlineIfNecessary()
				if solver.PrintDiagnosticInfo {
					println("tried all permutations")
				}
				break
			}
			continue
		}
		solver.distinctPool[childGenes] = true

		child := sequenceInfo{childGenes, getFitness(childGenes)}

		if !solver.childFitnessIsSameOrBetter(child, solver.pool[len(solver.pool)-1]) {
			continue
		}

		if child.fitness == solver.pool[len(solver.pool)-1].fitness {
			if len(solver.pool) < solver.maxPoolSize {
				solver.pool = append(solver.pool, child)
			} else {
				solver.pool[len(solver.pool)-1] = child
				insertionSort(solver.pool, solver.childFitnessIsSameOrBetter, len(solver.pool)-1)
			}

			continue
		}

		if len(children) < solver.maxPoolSize &&
			(len(distinctChildrenFitnesses) < 4 ||
				child.fitness == children[len(children)-1].fitness) {

			children = append(children, child)

			if solver.PrintDiagnosticInfo {
				print(".")
				solver.needNewlineBeforeDisplay = true
			}
		} else {
			children[len(children)-1] = child
		}

		insertionSort(children, solver.childFitnessIsSameOrBetter, len(children)-1)

		distinctChildren[child.genes] = true
		distinctChildrenFitnesses[child.fitness] = true

		if solver.childFitnessIsBetter(child, bestParent) {
			solver.printNewlineIfNecessary()
			if solver.PrintStrategyUsage {
				print(strategy.name)
			}
			display(child)
			start = time.Now()
			bestParent = child
			solver.incrementStrategyUseCount(strategyIndex)
			if solver.PrintDiagnosticInfo {
				print("+")
				solver.needNewlineBeforeDisplay = true
			}

			solver.pool[len(solver.pool)-1] = child
			insertionSort(solver.pool, solver.childFitnessIsSameOrBetter, len(solver.pool)-1)
		}

		if len(distinctChildren) >= solver.maxPoolSize &&
			len(distinctChildrenFitnesses) > 3 {
			if solver.PrintDiagnosticInfo {
				print(">")
				solver.needNewlineBeforeDisplay = true
			}

			solver.pool = children
			children = make([]sequenceInfo, 1, solver.maxPoolSize)
			children[0] = bestParent

			solver.distinctPool = distinctChildren
			distinctChildren = make(map[string]bool, len(children))
			distinctChildren[bestParent.genes] = true

			distinctChildrenFitnesses = make(map[int]bool, len(solver.pool))
			distinctChildrenFitnesses[bestParent.fitness] = true
		}
	}

	for _, child := range children {
		solver.pool = append(solver.pool, child)
	}

	return bestParent
}

func (solver *Solver) createFitnessComparisonFunctions() {
	if solver.LowerFitnessesAreBetter {
		solver.childFitnessIsBetter = func(child, other sequenceInfo) bool {
			return child.fitness < other.fitness
		}

		solver.childFitnessIsSameOrBetter = func(child, other sequenceInfo) bool {
			return child.fitness <= other.fitness
		}
	} else {
		solver.childFitnessIsBetter = func(child, other sequenceInfo) bool {
			return child.fitness > other.fitness
		}

		solver.childFitnessIsSameOrBetter = func(child, other sequenceInfo) bool {
			return child.fitness >= other.fitness
		}
	}
}

func (solver *Solver) ensureMaxSecondsToRunIsValid() {
	if solver.MaxSecondsToRunWithoutImprovement == 0 {
		solver.MaxSecondsToRunWithoutImprovement = 20
		fmt.Printf("\tSolver will run at most %v second(s) without improvement.\n", solver.MaxSecondsToRunWithoutImprovement)
	}
}

func (solver *Solver) getNextStrategyIndex() int {
	strategyIndex := solver.nextRand(solver.strategySuccessSum)
	for i, potentialStrategy := range solver.strategies {
		if strategyIndex < potentialStrategy.count {
			strategyIndex = i
			break
		}
		strategyIndex -= potentialStrategy.count
	}
	return strategyIndex
}

func (solver *Solver) incrementStrategyUseCount(strategyIndex int) {
	solver.strategies[strategyIndex].count++
	solver.strategySuccessSum++
}

func (solver *Solver) initialize(geneSet string, numberOfGenesPerChromosome int) {
	if solver.RandSeed == 0 {
		solver.RandSeed = time.Now().UnixNano()
	}
	solver.rand = createRandomNumberGenerator(solver.RandSeed)
	solver.ensureMaxSecondsToRunIsValid()
	solver.createFitnessComparisonFunctions()
	solver.initializeStrategies()
	solver.initializeChannels(geneSet, numberOfGenesPerChromosome)
	solver.needNewlineBeforeDisplay = false
}

func (solver *Solver) initializeChannels(geneSet string, numberOfGenesPerChromosome int) {
	solver.quit = false

	solver.nextGene = make(chan string, 1+numberOfGenesPerChromosome)
	go generateGene(solver.nextGene, geneSet, &solver.quit, solver.RandSeed)

	solver.nextChromosome = make(chan string, 1)
	go generateChromosome(solver.nextChromosome, solver.nextGene, geneSet, numberOfGenesPerChromosome, &solver.quit)
}

func (solver *Solver) initializePool(numberOfChromosomes, numberOfGenesPerChromosome int, geneSet string, initialParent sequenceInfo, getFitness func(string) int) {
	solver.maxPoolSize = max(len(geneSet), 3*numberOfChromosomes*numberOfGenesPerChromosome)
	solver.pool = make([]sequenceInfo, solver.maxPoolSize, solver.maxPoolSize)
	solver.pool[0] = initialParent
	solver.distinctPool = populatePool(solver.pool, solver.nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome, solver.childFitnessIsBetter, getFitness)
}

func (solver *Solver) initializeStrategies() {
	solver.strategies = []strategyInfo{
		{"crossover ", crossover, 1, false, createRandomNumberGenerator(solver.RandSeed)},
		{"mutate    ", mutate, 1, true, createRandomNumberGenerator(solver.RandSeed)},
		{"reverse   ", reverse, 1, true, createRandomNumberGenerator(solver.RandSeed)},
		{"shift     ", shift, 1, true, createRandomNumberGenerator(solver.RandSeed)},
		{"swap      ", swap, 1, true, createRandomNumberGenerator(solver.RandSeed)},
	}
	solver.strategySuccessSum = len(solver.strategies)
	solver.mutationStrategy = solver.strategies[1]
}

func (solver *Solver) nextRand(limit int) int {
	return solver.rand.Intn(limit)
}

func (solver *Solver) printNewlineIfNecessary() {
	if solver.needNewlineBeforeDisplay {
		solver.needNewlineBeforeDisplay = false
		println()
	}
}

func (solver *Solver) printStrategyUsage() {
	if !solver.PrintStrategyUsage {
		return
	}
	println("\nstrategy usage:")
	for _, strategy := range solver.strategies {
		println(fmt.Sprint(
			strategy.name, "\t",
			strategy.count, "\t",
			100.0*strategy.count/solver.strategySuccessSum, "%"))
	}
	println()
}

func populateDistinctPoolFitnessesMap(pool []sequenceInfo) map[int]bool {
	distinctChildrenFitnesses := make(map[int]bool, len(pool))
	distinctChildrenFitnesses[pool[0].fitness] = true
	distinctChildrenFitnesses[pool[len(pool)/2].fitness] = true
	distinctChildrenFitnesses[pool[len(pool)-1].fitness] = true
	return distinctChildrenFitnesses
}
