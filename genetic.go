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

	solver.initializePool(generationCount, numberOfGenesPerChromosome, geneSet, bestEver, getFitness)
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

	filteredDisplay := make(chan *sequenceInfo)
	go func() {
		for {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case candidate := <-filteredDisplay:
				if solver.childFitnessIsBetter(*candidate, bestEver) {
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

	bestEver = solver.pool[1]
	{
		tempBest := solver.pool[0]
		filteredDisplay <- &tempBest
	}

	maxLength := maxNumberOfChromosomes * numberOfGenesPerChromosome

	for len(bestEver.genes) <= maxLength &&
		roundsSinceLastImprovement < solver.MaxRoundsWithoutImprovement &&
		bestEver.fitness != bestPossibleFitness &&
		len(solver.pool) > 0 {

		result := solver.getBestWithInitialParent(getFitness,
			filteredDisplay,
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
			for _, parent := range solver.pool {
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

		solver.pool = newPool
		solver.distinctPool = distinctPool
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

	solver.initializePool(numberOfChromosomes, numberOfGenesPerChromosome, geneSet, initialParent, getFitness)
	solver.initializeStrategies(numberOfGenesPerChromosome, getFitness)

	best := *new(sequenceInfo)
	displayCaptureBest := make(chan *sequenceInfo)
	go func() {
		for {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case candidate := <-displayCaptureBest:
				solver.printNewlineIfNecessary()
				if solver.PrintStrategyUsage {
					fmt.Print((*candidate).strategy.name)
				}
				display((*candidate).genes)
				best = *candidate
			}
		}
	}()

	solver.getBestWithInitialParent(getFitness,
		displayCaptureBest,
		geneSet,
		numberOfChromosomes,
		numberOfGenesPerChromosome)

	solver.printNewlineIfNecessary()
	solver.printStrategyUsage()

	return best.genes
}

func (solver *Solver) getBestWithInitialParent(getFitness func(string) int,
	display chan *sequenceInfo,
	geneSet string,
	numberOfChromosomes, numberOfGenesPerChromosome int) sequenceInfo {

	start := time.Now()

	children := make([]sequenceInfo, 1, solver.maxPoolSize)
	children[0] = solver.pool[0]

	distinctChildren := make(map[string]bool, len(solver.pool))
	distinctChildrenFitnesses := populateDistinctPoolFitnessesMap(solver.pool)

	quit := make(chan bool)

	promoteChildrenIfFull := func() {
		if len(children) >= 20 || len(children) >= 10 && time.Since(start).Seconds() > solver.MaxSecondsToRunWithoutImprovement/2 {
			if solver.PrintDiagnosticInfo {
				fmt.Print(">")
				solver.needNewlineBeforeDisplay = true
			}

			solver.poolLock.Lock()
			solver.distinctPoolLock.Lock()
			solver.pool = children
			solver.distinctPool = distinctChildren
			solver.distinctPoolLock.Unlock()
			solver.poolLock.Unlock()

			bestParent := solver.pool[0]
			children = make([]sequenceInfo, 1, solver.maxPoolSize)
			children[0] = bestParent

			distinctChildren = make(map[string]bool, len(children))
			distinctChildren[bestParent.genes] = true

			distinctChildrenFitnesses = make(map[int]bool, len(solver.pool))
			distinctChildrenFitnesses[bestParent.fitness] = true
		}
	}

	updatePools := func(child *sequenceInfo) bool {
		addToChildren := func(item *sequenceInfo) {
			if len(children) < solver.maxPoolSize &&
				(len(distinctChildrenFitnesses) < 4 ||
					(*item).fitness == children[len(children)-1].fitness) {

				children = append(children, *item)

				if solver.PrintDiagnosticInfo {
					fmt.Print(".")
					solver.needNewlineBeforeDisplay = true
				}
				insertionSort(children, solver.childFitnessIsSameOrBetter, len(children)-1)
			} else if solver.childFitnessIsSameOrBetter(*item, children[len(children)-1]) {
				children[len(children)-1] = *item
				insertionSort(children, solver.childFitnessIsSameOrBetter, len(children)-1)
			}

			distinctChildren[(*item).genes] = true
			distinctChildrenFitnesses[(*item).fitness] = true
		}
		addToChildren(child)

		if solver.childFitnessIsBetter(*child, solver.pool[0]) {

			go func() {
				display <- child
			}()

			if solver.pool[0].genes == (*(*child).parent).genes {
				solver.strategySuccessLock.Lock()
				solver.successParentIsBestParentCount++
				solver.strategySuccessLock.Unlock()
			}
			solver.strategySuccessLock.Lock()
			solver.numberOfImprovements++
			solver.strategySuccessLock.Unlock()

			if solver.PrintDiagnosticInfo {
				fmt.Print("+")
				solver.needNewlineBeforeDisplay = true
			}

			solver.poolLock.Lock()
			solver.pool[len(solver.pool)-1] = *child
			insertionSort(solver.pool, solver.childFitnessIsSameOrBetter, len(solver.pool)-1)
			solver.poolLock.Unlock()

			if !distinctChildren[(*(*child).parent).genes] {
				addToChildren(child.parent)
			}

			return true
		}

		return false
	}

	updateIfIsImprovement := func(child *sequenceInfo) {
		if solver.inPool((*child).genes) {
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
		for _, child := range children {
			solver.pool = append(solver.pool, child)
		}
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
					return solver.pool[0]
				}
				promoteChildrenIfFull()
			}
		}
	}
	return solver.pool[0]
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

func (solver *Solver) initializePool(numberOfChromosomes, numberOfGenesPerChromosome int, geneSet string, initialParent sequenceInfo, getFitness func(string) int) {
	solver.maxPoolSize = getMaxPoolSize(numberOfChromosomes, numberOfGenesPerChromosome, len(geneSet))
	solver.pool = make([]sequenceInfo, solver.maxPoolSize, solver.maxPoolSize)
	solver.distinctPool = populatePool(solver.pool, solver.nextChromosome, geneSet, numberOfChromosomes, numberOfGenesPerChromosome, solver.childFitnessIsBetter, getFitness, initialParent)

	solver.numberOfImprovements = 1
	solver.randomParent = make(chan *sequenceInfo, 10)
	go func() {
		for {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			default:
				useBestParent := solver.random.Intn(solver.numberOfImprovements) <= solver.successParentIsBestParentCount
				if useBestParent {
					parent := solver.pool[0]
					select {
					case <-solver.quit:
						solver.quit <- true
						return
					case solver.randomParent <- &parent:
					}
				}
				parent := solver.pool[solver.random.Intn(len(solver.pool))]
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
	if len(solver.pool) == 0 {
		return false
	}
	if !solver.childFitnessIsSameOrBetter(*child, solver.pool[len(solver.pool)-1]) {
		return false
	}

	if (*child).fitness == solver.pool[len(solver.pool)-1].fitness {
		if len(solver.pool) < solver.maxPoolSize {
			solver.poolLock.Lock()
			solver.pool = append(solver.pool, *child)
			solver.poolLock.Unlock()
		} else {
			solver.poolLock.Lock()
			solver.pool[len(solver.pool)-1] = *child
			insertionSort(solver.pool, solver.childFitnessIsSameOrBetter, len(solver.pool)-1)
			solver.poolLock.Unlock()
		}

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

func populateDistinctPoolFitnessesMap(pool []sequenceInfo) map[int]bool {
	distinctChildrenFitnesses := make(map[int]bool, len(pool))
	distinctChildrenFitnesses[pool[0].fitness] = true
	distinctChildrenFitnesses[pool[len(pool)/2].fitness] = true
	distinctChildrenFitnesses[pool[len(pool)-1].fitness] = true
	return distinctChildrenFitnesses
}
