package genetic

import (
	"fmt"
	"math"
)

type Solver struct {
	MaxSecondsToRunWithoutImprovement float64
	MaxRoundsWithoutImprovement       int
	LowerFitnessesAreBetter           bool
	PrintStrategyUsage                bool
	PrintDiagnosticInfo               bool

	initialParentGenes             string
	initialParent                  sequenceInfo
	strategies                     map[string]*strategyInfo
	successParentIsBestParentCount int
	numberOfImprovements           int

	childFitnessIsBetter, childFitnessIsSameOrBetter func(child, other sequenceInfo) bool
}

func (solver *Solver) GetBest(getFitness func(string) int,
	display func(string),
	geneSet string,
	numberOfChromosomes, numberOfGenesPerChromosome int) string {

	quit := make(chan bool)
	solver.initialize(getFitness, -1, false)

	defer func() {
		quit <- true
		solver.initialParentGenes = ""
	}()

	bestEver := solver.initialParent
	displayCaptureBest := make(chan *sequenceInfo)

	go func() {
		for {
			select {
			case <-quit:
				return
			case candidate := <-displayCaptureBest:
				if !solver.childFitnessIsBetter(*candidate, bestEver) {
					continue
				}
				if solver.PrintDiagnosticInfo {
					fmt.Print((*candidate).strategy.name)
				}
				display((*candidate).genes)

				solver.incrementStrategyUseCount(candidate, &bestEver)

				bestEver = *candidate
			}
		}
	}()

	e := evolver{
		maxSecondsToRunWithoutImprovement: solver.MaxSecondsToRunWithoutImprovement,
		maxRoundsWithoutImprovement:       solver.MaxRoundsWithoutImprovement,
		lowerFitnessesAreBetter:           solver.LowerFitnessesAreBetter,
		childFitnessIsBetter:              solver.childFitnessIsBetter,
		childFitnessIsSameOrBetter:        solver.childFitnessIsSameOrBetter,
		geneSet:                           geneSet,
		numberOfGenesPerChromosome:        numberOfGenesPerChromosome,
		initialParent:                     solver.initialParent,
		display:                           displayCaptureBest,
		getFitness:                        getFitness,
	}
	e.getBest(numberOfChromosomes)

	solver.printStrategyUsage()

	return bestEver.genes
}

func (solver *Solver) GetBestUsingHillClimbing(getFitness func(string) int,
	display func(string),
	geneSet string,
	maxNumberOfChromosomes, numberOfGenesPerChromosome int,
	bestPossibleFitness int) string {

	quit := make(chan bool)
	solver.initialize(getFitness, bestPossibleFitness, true)

	defer func() {
		quit <- true
		solver.initialParentGenes = ""
	}()

	bestEver := solver.initialParent
	displayCaptureBest := make(chan *sequenceInfo)

	go func() {
		for {
			select {
			case <-quit:
				return
			case candidate := <-displayCaptureBest:
				if !solver.childFitnessIsBetter(*candidate, bestEver) {
					continue
				}
				if solver.PrintDiagnosticInfo {
					fmt.Print((*candidate).strategy.name)
				}
				display((*candidate).genes)

				solver.incrementStrategyUseCount(candidate, &bestEver)

				bestEver = *candidate
			}
		}
	}()

	e := evolver{
		maxSecondsToRunWithoutImprovement: solver.MaxSecondsToRunWithoutImprovement,
		maxRoundsWithoutImprovement:       solver.MaxRoundsWithoutImprovement,
		lowerFitnessesAreBetter:           solver.LowerFitnessesAreBetter,
		childFitnessIsBetter:              solver.childFitnessIsBetter,
		childFitnessIsSameOrBetter:        solver.childFitnessIsSameOrBetter,
		geneSet:                           geneSet,
		numberOfGenesPerChromosome:        numberOfGenesPerChromosome,
		initialParent:                     solver.initialParent,
		display:                           displayCaptureBest,
		getFitness:                        getFitness,
	}

	e.getBestUsingHillClimbing(maxNumberOfChromosomes, bestPossibleFitness)

	solver.printStrategyUsage()

	return bestEver.genes
}

func (solver *Solver) With(initialParentGenes string) *Solver {
	solver.initialParentGenes = initialParentGenes
	return solver
}

func (solver *Solver) createFitnessComparisonFunctions(bestPossibleFitness int, isHillClimbing bool) {
	if !isHillClimbing {
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

		if solver.LowerFitnessesAreBetter {
			solver.childFitnessIsBetter = func(child, other sequenceInfo) bool {
				eitherIsInvalid, toReturn := checkIfEitherIsInvalid(child.fitness, other.fitness)
				if eitherIsInvalid {
					return toReturn
				}

				childVsOptimalLower, childVsOptimalHigher := sort(child.fitness, bestPossibleFitness)
				otherVsOptimalLower, otherVsOptimalHigher := sort(other.fitness, bestPossibleFitness)
				if childVsOptimalHigher-childVsOptimalLower < otherVsOptimalHigher-otherVsOptimalLower {
					return child.fitness >= bestPossibleFitness
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
					return child.fitness >= bestPossibleFitness
				}
				return false
			}
		} else {
			solver.childFitnessIsBetter = func(child, other sequenceInfo) bool {
				eitherIsInvalid, toReturn := checkIfEitherIsInvalid(child.fitness, other.fitness)
				if eitherIsInvalid {
					return toReturn
				}

				childVsOptimalLower, childVsOptimalHigher := sort(child.fitness, bestPossibleFitness)
				otherVsOptimalLower, otherVsOptimalHigher := sort(other.fitness, bestPossibleFitness)
				if childVsOptimalHigher-childVsOptimalLower < otherVsOptimalHigher-otherVsOptimalLower {
					return child.fitness <= bestPossibleFitness
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
					return child.fitness <= bestPossibleFitness
				}
				return false
			}
		}
	}
}

func (solver *Solver) ensureMaxSecondsToRunIsValid() {
	if solver.MaxSecondsToRunWithoutImprovement == 0 {
		solver.MaxSecondsToRunWithoutImprovement = 20
		fmt.Printf("\tSolver will run at most %v second(s) without improvement.\n", solver.MaxSecondsToRunWithoutImprovement)
	}
}

func (solver *Solver) incrementStrategyUseCount(candidate, bestEver *sequenceInfo) {
	if bestEver.genes == (*(*candidate).parent).genes {
		solver.successParentIsBestParentCount++
	}
	solver.numberOfImprovements++

	strategyName := (*candidate).strategy.name
	strategy, exists := solver.strategies[strategyName]
	if !exists {
		strategy = &strategyInfo{name: strategyName}
		solver.strategies[strategyName] = strategy
	}
	(*strategy).successCount++
}

func (solver *Solver) initialize(getFitness func(string) int, optimalFitness int, isHillClimbing bool) {
	if solver.MaxRoundsWithoutImprovement == 0 {
		solver.MaxRoundsWithoutImprovement = 2
	}
	solver.ensureMaxSecondsToRunIsValid()
	solver.createFitnessComparisonFunctions(optimalFitness, isHillClimbing)

	solver.strategies = make(map[string]*strategyInfo, 10)

	initialParent := sequenceInfo{genes: solver.initialParentGenes}
	if len(initialParent.genes) == 0 {
		if solver.LowerFitnessesAreBetter {
			initialParent.fitness = math.MaxInt32
		} else {
			initialParent.fitness = math.MinInt32
		}
	} else {
		initialParent.fitness = getFitness(solver.initialParent.genes)
	}
	initialParent.parent = &solver.initialParent
	solver.initialParent = initialParent

}

func (solver *Solver) printStrategyUsage() {
	if !solver.PrintStrategyUsage {
		return
	}

	fmt.Println("\nsuccessful strategy usage:")
	for _, strategy := range solver.strategies {
		fmt.Println(
			strategy.name, "\t",
			strategy.successCount, "\t",
			100*strategy.successCount/solver.numberOfImprovements, "%")
	}
	fmt.Println()

	fmt.Println("\nNew champions were children of the reigning champion",
		100*solver.successParentIsBestParentCount/solver.numberOfImprovements,
		"% of the time.")
}
