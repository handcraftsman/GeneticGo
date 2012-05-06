package genetic

import (
	"fmt"
	"time"
)

func (solver *Solver) GetBestUsingHillClimbing(getFitness func(string) int,
	display func(string),
	geneSet string,
	maxNumberOfChromosomes, numberOfGenesPerChromosome int,
	bestPossibleFitness int) string {
	solver.isHillClimbing = true
	solver.geneSet = geneSet

	solver.initialize(geneSet, numberOfGenesPerChromosome, getFitness, bestPossibleFitness)

	roundsSinceLastImprovement := 0
	generationCount := 1

	bestEver := sequenceInfo{genes: solver.initialParent}
	if len(solver.initialParent) == 0 {
		bestEver = sequenceInfo{genes: generateParent(solver.nextChromosome, geneSet, generationCount, numberOfGenesPerChromosome)}
	}
	bestEver.fitness = getFitness(bestEver.genes)

	filteredDisplay := make(chan *sequenceInfo)

	go func() {
		for {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case candidate := <-filteredDisplay:
				if solver.childFitnessIsBetter(*candidate, bestEver) {
					if solver.PrintDiagnosticInfo {
						fmt.Print("+")
						solver.needNewlineBeforeDisplay = true
					}
					solver.printNewlineIfNecessary()
					if solver.PrintStrategyUsage {
						fmt.Print((*candidate).strategy.name)
					}
					display((*candidate).genes)
					bestEver = *candidate
					roundsSinceLastImprovement = 0
				}
			}
		}
	}()

	solver.initializePool(generationCount, numberOfGenesPerChromosome, geneSet, bestEver, getFitness, filteredDisplay)
	solver.initializeStrategies(numberOfGenesPerChromosome, getFitness)

	defer func() {
		solver.quit <- true
		solver.quit <- true
		solver.quit <- true
		<-solver.nextChromosome
		<-solver.nextGene
		for _, strategy := range solver.strategies {
			select {
			case <-strategy.results:
			default:
			}
		}
		solver.initialParent = ""
	}()

	maxLength := maxNumberOfChromosomes * numberOfGenesPerChromosome

	for len(bestEver.genes) <= maxLength &&
		roundsSinceLastImprovement < solver.MaxRoundsWithoutImprovement &&
		bestEver.fitness != bestPossibleFitness &&
		solver.pool.any() {

		result := solver.getBestWithInitialParent(getFitness,
			geneSet,
			len(bestEver.genes)/numberOfGenesPerChromosome,
			numberOfGenesPerChromosome)

		if solver.childFitnessIsBetter(result, bestEver) {
			roundsSinceLastImprovement = 0
			bestEver = result
			if bestEver.fitness == bestPossibleFitness {
				break
			}
		} else {
			roundsSinceLastImprovement++
			if roundsSinceLastImprovement >= solver.MaxRoundsWithoutImprovement {
				break
			}
		}

		generationCount++
		solver.printNewlineIfNecessary()

		if len(bestEver.genes) == maxLength {
			continue
		}

		solver.maxPoolSize = getMaxPoolSize(len(bestEver.genes)/numberOfGenesPerChromosome+1, numberOfGenesPerChromosome, len(geneSet))

		newPool := make([]sequenceInfo, 0, solver.maxPoolSize)
		distinctPool := make(map[string]bool, solver.maxPoolSize)

		improved := false
		climbStrategy := strategyInfo{name: "climb     "}

		for round := 0; round < 100 && !improved; round++ {
			for _, parent := range solver.pool.items {
				if len(parent.genes) >= maxLength {
					continue
				}
				childGenes := parent.genes + <-solver.nextChromosome
				if distinctPool[childGenes] {
					continue
				}
				distinctPool[childGenes] = true

				fitness := getFitness(childGenes)
				child := sequenceInfo{genes: childGenes, fitness: fitness, strategy: climbStrategy}
				if len(newPool) < solver.maxPoolSize {
					newPool = append(newPool, child)
				} else {
					newPool[len(newPool)-1] = child
				}
				insertionSort(newPool, solver.childFitnessIsSameOrBetter, len(newPool)-1)

				if solver.childFitnessIsBetter(child, bestEver) {
					improved = true
					filteredDisplay <- &child
				}
			}
		}

		solver.pool.truncateAndAddAll(newPool)
	}

	solver.printNewlineIfNecessary()
	solver.printStrategyUsage()

	return bestEver.genes
}

func (solver *Solver) With(initialParent string) *Solver {
	solver.initialParent = initialParent
	return solver
}

func (solver *Solver) GetBest(getFitness func(string) int,
	display func(string),
	geneSet string,
	numberOfChromosomes, numberOfGenesPerChromosome int) string {
	solver.isHillClimbing = false
	solver.geneSet = geneSet

	solver.initialize(geneSet, numberOfGenesPerChromosome, getFitness, -1)

	defer func() {
		solver.quit <- true
		<-solver.nextChromosome
		<-solver.nextGene
		for _, strategy := range solver.strategies {
			select {
			case <-strategy.results:
			default:
			}
		}
		solver.initialParent = ""
	}()

	initialParent := sequenceInfo{genes: solver.initialParent}
	if len(solver.initialParent) == 0 {
		initialParent = sequenceInfo{genes: generateParent(solver.nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome)}
	}
	initialParent.fitness = getFitness(initialParent.genes)

	displayCaptureBest := make(chan *sequenceInfo)

	best := *new(sequenceInfo)
	go func() {
		for {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case candidate := <-displayCaptureBest:
				if solver.PrintDiagnosticInfo {
					fmt.Print("+")
					solver.needNewlineBeforeDisplay = true
				}
				solver.printNewlineIfNecessary()
				if solver.PrintStrategyUsage {
					fmt.Print((*candidate).strategy.name)
				}
				display((*candidate).genes)
				best = *candidate
			}
		}
	}()

	solver.initializePool(numberOfChromosomes, numberOfGenesPerChromosome, geneSet, initialParent, getFitness, displayCaptureBest)
	solver.initializeStrategies(numberOfGenesPerChromosome, getFitness)

	solver.getBestWithInitialParent(getFitness,
		geneSet,
		numberOfChromosomes,
		numberOfGenesPerChromosome)

	solver.printNewlineIfNecessary()
	solver.printStrategyUsage()

	return best.genes
}

func (solver *Solver) getBestWithInitialParent(getFitness func(string) int,
	geneSet string,
	numberOfChromosomes, numberOfGenesPerChromosome int) sequenceInfo {

	start := time.Now()

	quit := make(chan bool)

	children := NewPool(solver.maxPoolSize,
		quit,
		solver.PrintDiagnosticInfo,
		solver.childFitnessIsSameOrBetter,
		func() { solver.needNewlineBeforeDisplay = true },
		solver.pool.addNewItem)
	poolBest := solver.pool.getBest()
	children.addNewItem <- &poolBest

	promoteChildrenIfFull := func() {
		if children.len() >= 20 || children.len() >= 10 && time.Since(start).Seconds() > solver.MaxSecondsToRunWithoutImprovement/2 {
			if solver.PrintDiagnosticInfo {
				fmt.Print(">")
				solver.needNewlineBeforeDisplay = true
			}

			solver.pool.truncateAndAddAll(children.items)

			bestParent := solver.pool.getBest()
			children.reset()
			children.addItem(bestParent)
		}
	}

	updatePools := func(child *sequenceInfo) bool {
		children.addItem(*child)

		poolBest := solver.pool.getBest()
		if solver.childFitnessIsBetter(*child, poolBest) {

			if poolBest.genes == (*(*child).parent).genes {
				solver.strategySuccessLock.Lock()
				solver.successParentIsBestParentCount++
				solver.strategySuccessLock.Unlock()
			}
			solver.strategySuccessLock.Lock()
			solver.numberOfImprovements++
			solver.strategySuccessLock.Unlock()

			solver.pool.addItem(*child)

			children.addItem(*(*child).parent)

			return true
		}

		return false
	}

	updateIfIsImprovement := func(child *sequenceInfo) {
		if solver.pool.contains(*child) {
			return
		}
		go func() {
			if solver.shouldAddChild(child, getFitness) {
				if updatePools(child) {
					solver.incrementStrategyUseCount((*child).strategy.index)
					start = time.Now()
				}
			}
		}()
	}

	timeout := make(chan bool, 1)
	go func() {
		for {
			time.Sleep(1 * time.Millisecond)
			select {
			case timeout <- true:
			case <-quit:
				quit <- true
				close(timeout)
				return
			}
		}
		close(timeout)
	}()

	defer func() {
		quit <- true
		solver.pool.addAll(children.items)
	}()

	for {
		// prefer successful strategies
		minStrategySuccess := solver.nextRand(solver.maxStrategySuccess)
		for index := 0; index < len(solver.strategies); index++ {
			if solver.strategies[index].successCount < minStrategySuccess {
				continue
			}
			select {
			case child := <-solver.strategies[index].results:
				updateIfIsImprovement(child)
			case <-timeout:
				if time.Since(start).Seconds() >= solver.MaxSecondsToRunWithoutImprovement {
					return solver.pool.getBest()
				}
				promoteChildrenIfFull()
			}
		}
	}
	return solver.pool.getBest()
}

func (solver *Solver) createFitnessComparisonFunctions(bestPossibleFitness int) {
	if !solver.isHillClimbing {
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
	} else {
		// checks distance from optimal
		// assumes negative fitnesses indicate invalid sequences

		checkIfEitherIsInvalid := func(childFitness, otherFitness int) (bool, bool) {
			if childFitness < 0 {
				if otherFitness < 0 {
					// both invalid, keep the newer one
					return true, true
				} else {
					// child is invalid but other is valid, keep it
					return true, false
				}
			} else if otherFitness < 0 {
				// child is valid but other is invalid, keep child
				return true, true
			}
			return false, false
		}

		solver.childFitnessIsBetter = func(child, other sequenceInfo) bool {
			eitherIsInvalid, toReturn := checkIfEitherIsInvalid(child.fitness, other.fitness)
			if eitherIsInvalid {
				return toReturn
			}

			childVsOptimalLower, childVsOptimalHigher := sort(child.fitness, bestPossibleFitness)
			otherVsOptimalLower, otherVsOptimalHigher := sort(other.fitness, bestPossibleFitness)
			if childVsOptimalHigher-childVsOptimalLower < otherVsOptimalHigher-otherVsOptimalLower {
				return solver.LowerFitnessesAreBetter && child.fitness >= bestPossibleFitness ||
					!solver.LowerFitnessesAreBetter && child.fitness <= bestPossibleFitness
			}
			return false
		}

		solver.childFitnessIsSameOrBetter = func(child, other sequenceInfo) bool {
			eitherIsInvalid, toReturn := checkIfEitherIsInvalid(child.fitness, other.fitness)
			if eitherIsInvalid {
				return toReturn
			}

			if child.fitness == bestPossibleFitness && other.fitness == bestPossibleFitness {
				// prefer the shorter optimal solution
				return len(child.genes) <= len(other.genes)
			}

			childVsOptimalLower, childVsOptimalHigher := sort(child.fitness, bestPossibleFitness)
			otherVsOptimalLower, otherVsOptimalHigher := sort(other.fitness, bestPossibleFitness)
			if childVsOptimalHigher-childVsOptimalLower <= otherVsOptimalHigher-otherVsOptimalLower {
				return solver.LowerFitnessesAreBetter && child.fitness >= bestPossibleFitness ||
					!solver.LowerFitnessesAreBetter && child.fitness <= bestPossibleFitness
			}
			return false
		}
	}
}

func (solver *Solver) ensureMaxSecondsToRunIsValid() {
	if solver.MaxSecondsToRunWithoutImprovement == 0 {
		solver.MaxSecondsToRunWithoutImprovement = 20
		fmt.Printf("\tSolver will run at most %v second(s) without improvement.\n", solver.MaxSecondsToRunWithoutImprovement)
	}
}

func (solver *Solver) incrementStrategyUseCount(strategyIndex int) {
	solver.strategySuccessLock.Lock()
	solver.strategies[strategyIndex].successCount++
	if solver.strategies[strategyIndex].successCount > solver.maxStrategySuccess {
		solver.maxStrategySuccess = solver.strategies[strategyIndex].successCount
	}
	solver.strategySuccessLock.Unlock()
}

func (solver *Solver) initialize(geneSet string, numberOfGenesPerChromosome int, getFitness func(string) int, optimalFitness int) {
	if solver.MaxRoundsWithoutImprovement == 0 {
		solver.MaxRoundsWithoutImprovement = 2
	}
	solver.random = createRandomNumberGenerator()
	solver.ensureMaxSecondsToRunIsValid()
	solver.createFitnessComparisonFunctions(optimalFitness)
	solver.initializeChannels(geneSet, numberOfGenesPerChromosome)
	solver.needNewlineBeforeDisplay = false
}

func (solver *Solver) initializeChannels(geneSet string, numberOfGenesPerChromosome int) {
	solver.quit = make(chan bool)
	solver.nextGene = make(chan string, 1+numberOfGenesPerChromosome)
	go generateGene(solver.nextGene, geneSet, solver.quit)

	solver.nextChromosome = make(chan string, 1)
	go generateChromosome(solver.nextChromosome, solver.nextGene, geneSet, numberOfGenesPerChromosome, solver.quit)
}

func (solver *Solver) initializePool(numberOfChromosomes, numberOfGenesPerChromosome int, geneSet string, initialParent sequenceInfo, getFitness func(string) int, display chan *sequenceInfo) {
	solver.maxPoolSize = getMaxPoolSize(numberOfChromosomes, numberOfGenesPerChromosome, len(geneSet))

	solver.pool = NewPool(solver.maxPoolSize,
		solver.quit,
		solver.PrintDiagnosticInfo,
		solver.childFitnessIsSameOrBetter,
		func() { solver.needNewlineBeforeDisplay = true },
		display)

	solver.pool.populatePool(solver.nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome, solver.childFitnessIsBetter, getFitness, initialParent)

	solver.numberOfImprovements = 1
	solver.randomParent = make(chan *sequenceInfo, 10)
	go func() {
		for {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			default:
				sendBestParent := solver.random.Intn(solver.numberOfImprovements) <= solver.successParentIsBestParentCount
				if sendBestParent {
					parent := solver.pool.getBest()
					select {
					case <-solver.quit:
						solver.quit <- true
						return
					case solver.randomParent <- &parent:
					}
				}
				parent := solver.pool.getRandomItem()
				select {
				case <-solver.quit:
					solver.quit <- true
					return
				case solver.randomParent <- &parent:
				}
			}
		}
	}()
}

func (solver *Solver) nextRand(limit int) int {
	return solver.random.Intn(limit)
}

func (solver *Solver) printNewlineIfNecessary() {
	if solver.needNewlineBeforeDisplay {
		solver.needNewlineBeforeDisplay = false
		fmt.Println()
	}
}

func (solver *Solver) printStrategyUsage() {
	if !solver.PrintStrategyUsage {
		return
	}

	fmt.Println("\nsuccessful strategy usage:")
	for _, strategy := range solver.strategies {
		successCount := strategy.successCount - initialStrategySuccess
		fmt.Println(
			strategy.name, "\t",
			successCount, "\t",
			100*successCount/solver.numberOfImprovements, "%")
	}
	fmt.Println()

	fmt.Println("\nNew champions were children of the reigning champion",
		100*solver.successParentIsBestParentCount/solver.numberOfImprovements,
		"% of the time.")
}

func (solver *Solver) shouldAddChild(child *sequenceInfo, getFitness func(string) int) bool {
	(*child).fitness = getFitness((*child).genes)

	if !solver.pool.any() {
		return false // already returned final result
	}
	poolWorst := solver.pool.getWorst()
	if !solver.childFitnessIsSameOrBetter(*child, poolWorst) {
		return false
	}

	if (*child).fitness == poolWorst.fitness {
		solver.pool.addItem(*child)
		return false
	}

	return true
}

func getMaxPoolSize(numberOfChromosomes, numberOfGenesPerChromosome, numberOfGenes int) int {
	max := numberOfGenes
	for i := 1; i < numberOfChromosomes*numberOfGenesPerChromosome && max < 500; i++ {
		max *= numberOfGenes
	}
	if max > 500 {
		return 500
	}
	return max
}

func populateDistinctPoolFitnessesMap(pool *pool) map[int]bool {
	distinctChildrenFitnesses := make(map[int]bool, pool.maxPoolSize)
	distinctChildrenFitnesses[pool.items[0].fitness] = true
	distinctChildrenFitnesses[pool.items[len(pool.items)/2].fitness] = true
	distinctChildrenFitnesses[pool.items[len(pool.items)-1].fitness] = true
	return distinctChildrenFitnesses
}
